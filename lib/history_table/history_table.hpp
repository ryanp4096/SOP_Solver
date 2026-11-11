/**
 * The history table stores prefixes and subpaths that have previously been visited by branch and bound
 */

#ifndef HISTORY_TABLE_H
#define HISTORY_TABLE_H

#include <immintrin.h>

#include <array>
#include <cstdint>
#include <cstddef>

#include "hash.hpp"
#include "memory.hpp"
#include "entries.hpp"

/* Change number of entries stored in each bucket */
constexpr size_t slots_per_bucket = 4;

constexpr size_t history_entry_size = std::max(sizeof(PrefixEntry), sizeof(SubpathEntry));
constexpr size_t history_entry_align = std::max(alignof(PrefixEntry), alignof(SubpathEntry));


/* Generic history entry, can store either prefix or subpath */
struct alignas(history_entry_align) HistoryEntryTemplate {
    char padding[history_entry_size];
};


/* Stores lower 96 bits of 128-bit hash to check if key matches entry */
using HashCheck = std::array<uint32_t, 3>;


/**
 * Represents a bucket in the history table. A bucket stores 4 history entries by default.
 * The upper 32 bits of the 128-bit hash of the history key are used to determine which
 * bucket a key will be found in.
 */
struct alignas(64) HistoryBucket {    
    uint8_t lock_byte;                                  // 0 = unlocked, 1 = locked
    uint16_t occupied_slot_mask;                         // 1 bit for each slot, 1 = occupied, 0 = empty
    HashCheck hash_checks[slots_per_bucket];            // the lower 96 bits of the hash of the entry stored in each slot
    HistoryEntryMetadata metadata[slots_per_bucket];    // extra data to choose which entry to replace when bucket is full
    HistoryEntryTemplate entries[slots_per_bucket];     // the entries in the bucket
};


/* Each bucket is aligned to fit on a certain number of 64-byte cache lines */
constexpr size_t cache_lines_per_bucket = sizeof(HistoryBucket) / 64;


/* The outcome of the retrieve operation */
enum HistoryOutcome {
    Found,          // Key was found in history table
    Inserted,       // Key was not found, inserted
    NotInserted     // Key was not found, could not be inserted
};


/* A reference to a previously retrieved history entry, to more efficiently re-access */
struct HistoryReference {
    size_t bucket = 0;  // bucket 
    int slot = -1;      // slot of the bucket
    HashCheck check;    // used to check if the slot still contains the same hash

    /* Whether this reference actually references a history entry */
    inline bool exists() const { return slot >= 0; }
};


/**
 * The history table stores prefixes and subpaths that have previously been visited
 * in branch and bound.
 * 
 * A prefix is a partial path starting at node 0 and ending at any node.
 * - If we already visited a prefix containing the same nodes with the same
 *   last node as the current prefix, and the old prefix had a better cost,
 *   then the current path can be pruned as it is known to be inferior.
 * - If the current prefix is better than the prefix in the history table,
 *   then we send a thread stop request to the thread that inserted the prefix,
 *   so if it is still running, it will stop as it is known to be inferior.
 * 
 * A subpath is a slice of a prefix that can start and end at any node.
 * - If we already visited a subpath containing the same nodes with the same
 *   first and last node as the current subpath, and the existing subpath had
 *   a better cost, the current path can be pruned as it is inferior.
 * 
 */
class HistoryTable {
private:
    static constexpr unsigned all_slots_mask = (1u << slots_per_bucket) - 1;

    HistoryTableMemory memory;      // Memory block for history table
    size_t bucket_count;            // Number of buckets in table
    HistoryBucket *buckets;         // Pointer to memory storing table

    /**
     * Lock for history table buckets.
     * Locks when constructed, releases when destroyed.
     */
    class BucketLock {
    private:
        uint8_t* lock_byte;
    public:
        explicit BucketLock(HistoryBucket& bucket) : lock_byte(&bucket.lock_byte) {
            for (;;) {
                // Try to take the lock with one atomic exchange
                if (__atomic_exchange_n(lock_byte, static_cast<uint8_t>(1), __ATOMIC_ACQUIRE) == 0)
                    return;
                
                // Someone else holds it. Spin with plain reads, which keep the cache
                // line shared instead of fighting over ownership, until it looks free
                while (__atomic_load_n(lock_byte, __ATOMIC_RELAXED) != 0)
                    _mm_pause();
            }
        }
        ~BucketLock() {
            __atomic_store_n(lock_byte, static_cast<uint8_t>(0), __ATOMIC_RELEASE);
        }
    };


public:
    /**
     * Create a new history table.
     * @param memory_bytes The size of the history table, in bytes
     */
    explicit HistoryTable(HistoryTableMemoryOptions options = {})
          : memory(options),
            bucket_count(memory.get_memory_bytes() / sizeof(HistoryBucket)),
            buckets(reinterpret_cast<HistoryBucket *>(memory.get_memory_start()))
            {}

    /**
     * Send request to start loading bucket into cache.
     * Can load in background, then be accessed quickly later when needed
     * @param hash The hash to retrieve
     */
    void prefetch(Hash128 hash) const {
        const char* bucket_start = reinterpret_cast<const char*>(&buckets[hash_to_bucket(hash)]);
        for (size_t i = 0; i < cache_lines_per_bucket; i++)
            __builtin_prefetch(bucket_start + i * 64, 1, 3);
    }

    /**
     * Send request to start loading bucket into cache.
     * Can load in background, then be accessed quickly later when needed
     * @param reference A HistoryReference to a previously retrieved history entry
     */
    void prefetch(const HistoryReference& reference) const {
        if (!reference.exists()) return;
        const char* bucket_start = reinterpret_cast<const char*>(&buckets[reference.bucket]);
        for (size_t i = 0; i < cache_lines_per_bucket; i++)
            __builtin_prefetch(bucket_start + i * 64, 1, 3);
    }

    /**
     * Look up a key in the table. If absent, insert the key into the table.
     * @param hash The hash to search for
     * @param callback A function to run immediately after retrieving/inserting,
     *                 while lock is still held, for quick operations only
     * @param reference Optional return variable. A pointer to a HistoryReference to store a reference to the
     *                  location of the retrieved/inserted entry for future access
     * @return Whether the entry was found or inserted
     */
    template <class EntryType, class Callback>
    HistoryOutcome retrieve_or_insert(Hash128 hash, size_t depth, Callback&& callback, HistoryReference *reference = nullptr) {
        size_t bucket_index = hash_to_bucket(hash);
        HistoryBucket &bucket = buckets[bucket_index];
        HashCheck check = hash_to_check(hash);
        BucketLock lock(bucket); // Locks bucket, automatically unlocks on return

        int slot = search_bucket(bucket, check);
        if (slot >= 0) {
            // Entry found in bucket
            callback(HistoryOutcome::Found, reinterpret_cast<EntryType *>(&bucket.entries[slot]));
            if (reference) *reference = {bucket_index, slot, check};
            return HistoryOutcome::Found;
        }

        const unsigned free_slot_mask = ~static_cast<unsigned>(bucket.occupied_slot_mask) & all_slots_mask;
        if (free_slot_mask == 0) {
            // No free slots
            callback(HistoryOutcome::NotInserted, nullptr);
            return HistoryOutcome::NotInserted;
        }

        // Insert in first free slot
        slot = __builtin_ctz(free_slot_mask); // First bit equal to 1 in free mask
        bucket.hash_checks[slot] = check;
        bucket.occupied_slot_mask |= static_cast<uint16_t>(1u << slot);
        bucket.metadata[slot].depth = depth;
        callback(HistoryOutcome::Inserted, reinterpret_cast<EntryType *>(&bucket.entries[slot]));
        if (reference) *reference = {bucket_index, slot, check};

        return HistoryOutcome::Inserted;
    }

    /**
     * Update a previously found or inserted entry from retrieve_or_insert.
     * @param reference A HistoryReference returned from retrieve_or_insert
     * @param callback A function to run immediately after finding entry, 
     *                 while lock is still held, for quick operations only
     * @return true if entry found and updated, false if entry is no longer present (replaced)
     */
    template <class EntryType, class Callback>
    bool update_entry(const HistoryReference &reference, Callback&& callback) {
        if (!reference.exists()) return false;

        HistoryBucket &bucket = buckets[reference.bucket];
        BucketLock lock(bucket); // Locks bucket until return

        if (bucket.hash_checks[reference.slot] == reference.check) {
            callback(true, reinterpret_cast<EntryType *>(&bucket.entries[reference.slot]));
            return true;
        } else {
            callback(false, nullptr);
            return false;
        }
    }

private:
    /**
     * Searches a bucket to see if it contains a hash.
     * @param bucket Bucket to search
     * @param check Hash to search for
     * @return Slot index containing the hash, or -1 if not found
     */
    int search_bucket(const HistoryBucket &bucket, const HashCheck &check) {
        unsigned slot_mask = bucket.occupied_slot_mask;
        while (slot_mask != 0) {
            int slot_index = __builtin_ctz(slot_mask);  // Get next occupied slot (first bit equal to 1)
            
            if (bucket.hash_checks[slot_index] == check)
                return slot_index;
            
            slot_mask &= slot_mask - 1; // Move to next slot (clear first bit equal to 1)
        }
        return -1;
    }

    /**
     * Get the bucket of a 128-bit hash, based on the upper 32 bits
     * @param hash A hash
     * @return The bucket that hash would be contained in
     */
    size_t hash_to_bucket(Hash128 hash) const {
        // converts high 32 bits of hash from range [0, 2^32) to range [0, bucket_count)
        // this code is equivalent to high_32_bits * (bucket_count / 2^32)
        return static_cast<size_t>(((hash.high >> 32) * bucket_count) >> 32);
    }

    /**
     * Get the lower 96 bits of a 128-bit hash, separated into three 32-bit segments.
     * Used to check for a match within a bucket.
     * @param hash A hash
     * @return The lower 96 bits of the hash, separated into 32-bit pieces
     */
    HashCheck hash_to_check(Hash128 hash) {
        return {
            static_cast<uint32_t>(hash.high),
            static_cast<uint32_t>(hash.low),
            static_cast<uint32_t>(hash.low >> 32)
        };
    }
};


#endif