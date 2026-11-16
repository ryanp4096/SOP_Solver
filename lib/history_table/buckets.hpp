/**
 * Defines how each bucket in the history table is stored
 */

#ifndef HISTORY_TABLE_BUCKET_H
#define HISTORY_TABLE_BUCKET_H

#include <immintrin.h>

#include <array>
#include <cstdint>
#include <cstddef>
#include <vector>

#include "hash.hpp"
#include "memory.hpp"
#include "entries.hpp"
#include "lock.hpp"
#include <mutex>

/* Storage configuration */
constexpr size_t slots_per_bucket = 8;          // Number of entries stored in each bucket
constexpr size_t slots_per_expansion = 4;       // Number of entries each bucket is expanded by when full
constexpr size_t expansion_pool_size = 10000;   // Number of expansion spaces allocated at once
constexpr double expansion_space = 0.25;        // Space left in memory to be used for expanding buckets


constexpr size_t history_entry_size = std::max(sizeof(PrefixEntry), sizeof(SubpathEntry));
constexpr size_t history_entry_align = std::max(alignof(PrefixEntry), alignof(SubpathEntry));

/* Generic history entry, can store either prefix or subpath */
struct alignas(history_entry_align) HistoryEntryTemplate {
    char padding[history_entry_size];
};

/* Stores lower 96 bits of 128-bit hash to check if key matches entry */
using HashCheck = std::array<uint32_t, 3>;


/**
 * Expansion of history bucket
 */
struct alignas(64) HistoryExpansion {
    uint16_t occupied_slot_mask;                           // 1 bit for each slot, 1 = occupied, 0 = empty
    HistoryExpansion *expansion;                           // if bucket is expanded, points to expansion
    HashCheck hash_checks[slots_per_expansion];            // the lower 96 bits of the hash of the entry stored in each slot
    HistoryEntryMetadata metadata[slots_per_expansion];    // extra data to choose which entry to replace when bucket is full
    HistoryEntryTemplate entries[slots_per_expansion];     // the entries in the bucket
};


/**
 * Represents a bucket in the history table. A bucket stores 4 history entries by default.
 * The upper 32 bits of the 128-bit hash of the history key are used to determine which
 * bucket a key will be found in.
 */
struct alignas(64) HistoryBucket {    
    uint8_t lock_byte;                                  // 0 = unlocked, 1 = locked
    uint8_t entry_count;
    uint16_t occupied_slot_mask;                         // 1 bit for each slot, 1 = occupied, 0 = empty
    HistoryExpansion *expansion;                        // if bucket is expanded, points to expansion
    HashCheck hash_checks[slots_per_bucket];            // the lower 96 bits of the hash of the entry stored in each slot
    HistoryEntryMetadata metadata[slots_per_bucket];    // extra data to choose which entry to replace when bucket is full
    HistoryEntryTemplate entries[slots_per_bucket];     // the entries in the bucket
};


/* Each bucket is aligned to fit on a certain number of 64-byte cache lines */
constexpr size_t cache_lines_per_bucket = sizeof(HistoryBucket) / 64;
constexpr size_t cache_lines_per_expansion = sizeof(HistoryExpansion) / 64;


#endif