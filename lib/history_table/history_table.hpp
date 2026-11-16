/**
 * The history table stores prefixes and subpaths that have previously been visited by branch and bound
 */

#ifndef HISTORY_TABLE_H
#define HISTORY_TABLE_H

#include <array>
#include <cstdint>
#include <cstddef>
#include <utility>

#include "hash.hpp"
#include "memory.hpp"
#include "entries.hpp"
#include "lock.hpp"
#include "buckets.hpp"
#include "expansion.hpp"


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
    static constexpr unsigned all_expansion_slots_mask = (1u << slots_per_expansion) - 1;

    HistoryTableMemory memory;      // Memory block for history table
    size_t bucket_count;            // Number of buckets in table
    HistoryBucket *buckets;         // Pointer to memory storing table
    ExpansionManager expansion_manager;     // Manager for expanding buckets

public:
    /**
     * Create a new history table.
     * @param memory_bytes The size of the history table, in bytes
     */
    explicit HistoryTable(HistoryTableMemoryOptions options = {})
          : memory(options),
            bucket_count(static_cast<size_t>(memory.get_memory_bytes() * (1 - expansion_space)) / sizeof(HistoryBucket)),
            buckets(reinterpret_cast<HistoryBucket *>(memory.get_memory_start())),
            expansion_manager(
                reinterpret_cast<HistoryExpansion *>(&buckets[bucket_count]),
                static_cast<size_t>(memory.get_memory_bytes() * expansion_space) / sizeof(HistoryExpansion),
                options.thread_count
            )
            {}

    /**
     * Get bucket count
     */
    size_t get_bucket_count() const { return bucket_count; }

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
    HistoryOutcome retrieve_or_insert_one_bucket(Hash128 hash, size_t depth, bool is_subpath, Callback&& callback, HistoryReference *reference = nullptr) {
        size_t bucket_index = hash_to_bucket(hash);
        HistoryBucket &bucket = buckets[bucket_index];
        HashCheck check = hash_to_check(hash);
        BucketLock lock(&bucket.lock_byte); // Locks bucket, automatically unlocks on return

        int slot = search_bucket(bucket, check);
        if (slot >= 0) {
            // Entry found in bucket
            bucket.metadata[slot].hits++;
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
        bucket.entry_count++;
        bucket.metadata[slot].type = HistoryEntryType(is_subpath, depth);
        bucket.metadata[slot].hits = 1;
        callback(HistoryOutcome::Inserted, reinterpret_cast<EntryType *>(&bucket.entries[slot]));
        if (reference) *reference = {bucket_index, slot, check};

        return HistoryOutcome::Inserted;
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
    HistoryOutcome retrieve_or_insert_two_bucket(Hash128 hash, size_t depth, bool is_subpath, Callback&& callback, HistoryReference *reference = nullptr) {
        size_t bucket_index1 = hash_to_bucket(hash);
        size_t bucket_index2 = other_bucket(bucket_index1, hash);
        HistoryBucket &bucket1 = buckets[bucket_index1];
        HistoryBucket &bucket2 = buckets[bucket_index2];

        size_t lower_index = bucket_index1; size_t higher_index = bucket_index2;
        if (lower_index > higher_index) std::swap(lower_index, higher_index);
        HistoryBucket &lower_bucket = buckets[lower_index];
        HistoryBucket &higher_bucket = buckets[higher_index];

        HashCheck check = hash_to_check(hash);
        BucketDoubleLock lock(&lower_bucket.lock_byte, &higher_bucket.lock_byte); // Locks lower bucket, locks higher bucket later. Both locks release on return

        int slot = search_bucket(lower_bucket, check);
        if (slot >= 0) {
            // Entry found in lower bucket
            lower_bucket.metadata[slot].hits++;
            callback(HistoryOutcome::Found, reinterpret_cast<EntryType *>(&lower_bucket.entries[slot]));
            if (reference) *reference = {lower_index, slot, check};
            return HistoryOutcome::Found;
        }

        if (lower_index != higher_index) {
            lock.lock_second();
            slot = search_bucket(higher_bucket, check);
            if (slot >= 0) {
                // Entry found in bucket 2
                higher_bucket.metadata[slot].hits++;
                callback(HistoryOutcome::Found, reinterpret_cast<EntryType *>(&higher_bucket.entries[slot]));
                if (reference) *reference = {higher_index, slot, check};
                return HistoryOutcome::Found;
            }
        }

        size_t bucket_index = bucket1.entry_count < bucket2.entry_count ? bucket_index1 : bucket_index2;
        HistoryBucket &bucket = buckets[bucket_index];
        
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
        bucket.entry_count++;
        bucket.metadata[slot].type = HistoryEntryType(is_subpath, depth);
        bucket.metadata[slot].hits = 1;
        callback(HistoryOutcome::Inserted, reinterpret_cast<EntryType *>(&bucket.entries[slot]));
        if (reference) *reference = {bucket_index, slot, check};

        return HistoryOutcome::Inserted;
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
    HistoryOutcome retrieve_or_insert_expand(Hash128 hash, size_t depth, bool is_subpath, size_t thread_id, Callback&& callback, HistoryReference *reference = nullptr) {
        size_t bucket_index = hash_to_bucket(hash);
        HistoryBucket &bucket = buckets[bucket_index];
        HashCheck check = hash_to_check(hash);
        BucketLock lock(&bucket.lock_byte); // Locks bucket, automatically unlocks on return

        int slot = search_bucket(bucket, check);
        if (slot >= 0) {
            // Entry found in bucket
            bucket.metadata[slot].hits++;
            callback(HistoryOutcome::Found, reinterpret_cast<EntryType *>(&bucket.entries[slot]));
            if (reference) *reference = {bucket_index, slot, check};
            return HistoryOutcome::Found;
        }
        HistoryExpansion *expansion = bucket.expansion;
        int slot_shift = slots_per_bucket;
        while (expansion != nullptr) {
            int slot = search_expansion(*expansion, check);
            if (slot >= 0) {
                expansion->metadata[slot].hits++;
                callback(HistoryOutcome::Found, reinterpret_cast<EntryType *>(&expansion->entries[slot]));
                if (reference) *reference = {bucket_index, slot_shift + slot, check};
                return HistoryOutcome::Found;
            }
            expansion = expansion->expansion;
            slot_shift += slots_per_expansion;
        }

        unsigned entry_count = bucket.entry_count;
        if (entry_count < slots_per_bucket) {
            const unsigned free_slot_mask = ~static_cast<unsigned>(bucket.occupied_slot_mask) & all_slots_mask;
            slot = __builtin_ctz(free_slot_mask); // First bit equal to 1 in free mask
            bucket.hash_checks[slot] = check;
            bucket.occupied_slot_mask |= static_cast<uint16_t>(1u << slot);
            bucket.entry_count++;
            bucket.metadata[slot].type = HistoryEntryType(is_subpath, depth);
            bucket.metadata[slot].hits = 1;
            callback(HistoryOutcome::Inserted, reinterpret_cast<EntryType *>(&bucket.entries[slot]));
            if (reference) *reference = {bucket_index, slot, check};
            return HistoryOutcome::Inserted;
        }
        if (bucket.expansion == nullptr) {
            bucket.expansion = expansion_manager.allocate_expansion(thread_id);
            if (bucket.expansion == nullptr) {
                // No free slots
                callback(HistoryOutcome::NotInserted, nullptr);
                return HistoryOutcome::NotInserted;
            }
        }
        slot_shift = slots_per_bucket;
        entry_count -= slots_per_bucket;
        expansion = bucket.expansion;
        while (entry_count >= slots_per_expansion) {
            if (expansion->expansion == nullptr) {
                expansion->expansion = expansion_manager.allocate_expansion(thread_id);
                if (expansion->expansion == nullptr) {
                    // No free slots
                    callback(HistoryOutcome::NotInserted, nullptr);
                    return HistoryOutcome::NotInserted;
                }
            }
            expansion = expansion->expansion;
            entry_count -= slots_per_expansion;
            slot_shift += slots_per_expansion;
        }
        const unsigned free_slot_mask = ~static_cast<unsigned>(expansion->occupied_slot_mask) & all_expansion_slots_mask;
        slot = __builtin_ctz(free_slot_mask); // First bit equal to 1 in free mask
        expansion->hash_checks[slot] = check;
        expansion->occupied_slot_mask |= static_cast<uint16_t>(1u << slot);
        bucket.entry_count++;
        expansion->metadata[slot].type = HistoryEntryType(is_subpath, depth);
        expansion->metadata[slot].hits = 1;
        callback(HistoryOutcome::Inserted, reinterpret_cast<EntryType *>(&expansion->entries[slot]));
        if (reference) *reference = {bucket_index, slot_shift + slot, check};
        return HistoryOutcome::Inserted;
    }

    template <class Callback>
    HistoryOutcome retrieve_or_insert_prefix(Hash128 hash, size_t depth, size_t thread_id, Callback&& callback, HistoryReference *reference = nullptr) {
        return retrieve_or_insert_expand<PrefixEntry>(hash, depth, false, thread_id, callback, reference);
    }

    template <class Callback>
    HistoryOutcome retrieve_or_insert_subpath(Hash128 hash, size_t depth, size_t thread_id, Callback&& callback, HistoryReference *reference = nullptr) {
        return retrieve_or_insert_expand<SubpathEntry>(hash, depth, true, thread_id, callback, reference);
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
        BucketLock lock(&bucket.lock_byte); // Locks bucket until return

        if (reference.slot < static_cast<int>(slots_per_bucket)) {
            if (bucket.hash_checks[reference.slot] == reference.check) {
                callback(true, reinterpret_cast<EntryType *>(&bucket.entries[reference.slot]));
                return true;
            } else {
                callback(false, nullptr);
                return false;
            }
        } else {
            int slot = reference.slot - slots_per_bucket;
            HistoryExpansion *expansion = bucket.expansion;
            while (expansion != nullptr && slot >= static_cast<int>(slots_per_expansion)) {
                expansion = expansion->expansion;
                slot -= slots_per_expansion;
            }
            if (expansion == nullptr || expansion->hash_checks[slot] != reference.check) {
                callback(false, nullptr);
                return false;
            } else {
                callback(true, reinterpret_cast<EntryType *>(&expansion->entries[slot]));
                return true;
            }
        }
    }

    /**
     * Diagnostic: distribution of entries in buckets
     * @return Vector containing number of buckets containing i elements for each index i
     */
    std::vector<size_t> bucket_distribution() const {
        std::vector<size_t> distribution(slots_per_bucket + 1, 0);
        for (size_t i = 0; i < bucket_count; i++) {
            size_t entries = buckets[i].entry_count;
            if (entries >= distribution.size()) {
                distribution.resize(entries + 1);
            }
            distribution[entries]++;
        }
        return distribution;
    }

    /**
     * Diagnostic: distribution of how many times each entry is hit
     */
    std::pair<std::vector<size_t>, std::vector<size_t>> hit_distribution() const {
        std::vector<size_t> prefix_dist(16, 0);
        std::vector<size_t> subpath_dist(16, 0);
        for (size_t i = 0; i < bucket_count; i++) {
            unsigned slot_mask = buckets[i].occupied_slot_mask;
            while (slot_mask != 0) {
                int j = __builtin_ctz(slot_mask);

                unsigned hits = buckets[i].metadata[j].hits;
                int block = hits ? (31 - __builtin_clz(hits)) : 0;

                if (buckets[i].metadata[j].type.is_subpath())
                    subpath_dist[block]++;
                else
                    prefix_dist[block]++;

                slot_mask &= slot_mask - 1;
            }
        }
        return std::make_pair(prefix_dist, subpath_dist);
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

    int search_expansion(const HistoryExpansion &expansion, const HashCheck &check) {
        unsigned slot_mask = expansion.occupied_slot_mask;
        while (slot_mask != 0) {
            int slot_index = __builtin_ctz(slot_mask);  // Get next occupied slot (first bit equal to 1)
            
            if (expansion.hash_checks[slot_index] == check)
                return slot_index;
            
            slot_mask &= slot_mask - 1; // Move to next slot (clear first bit equal to 1)
        }
        return -1;
    }

    /**
     * Get the bucket of a 128-bit hash, based on the upper ~32 bits
     * @param hash A hash
     * @return The bucket that hash would be contained in
     */
    size_t hash_to_bucket(Hash128 hash) const {
        // converts high bits of hash from range [0, 2^64) to range [0, bucket_count)
        // this code is equivalent to high_bits * (bucket_count / 2^64)
        return static_cast<uint64_t>((static_cast<unsigned __int128>(hash.high) * bucket_count) >> 64);
    }

    /**
     * Get other bucket from one bucket
     */
    size_t other_bucket(size_t bucket, Hash128 hash) const {
        size_t hash_bits =  static_cast<uint64_t>((static_cast<unsigned __int128>(hash.low) * bucket_count) >> 64);
        if (hash_bits < bucket) hash_bits += bucket_count;
        return hash_bits - bucket;
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