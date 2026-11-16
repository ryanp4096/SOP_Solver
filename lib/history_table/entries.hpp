/**
 * Defines what is stored for each entry in the history table
 */

#ifndef HISTORY_ENTRIES_H
#define HISTORY_ENTRIES_H

#include <cstdint>
#include <boost/dynamic_bitset.hpp>


/* Represents a key in the history table*/
struct HistoryKey {
    boost::dynamic_bitset<> bit_vector;
    int first_node = -1;
    int last_node = -1;
};


/* State of thread exploring a prefix */
enum HistoryState : uint8_t { UNEXPLORED, EXPLORING, EXPLORED };


/* Data stored in history table for a prefix */
struct PrefixEntry {
    int32_t cost;               // the cost of this prefix
    int32_t lower_bound;        // the lower bound on the cost of the entire path that starts with this prefix
    HistoryState state;    // status of exploring the subtree of the prefix
    int8_t active_thread;       // the thread currently exploring
};


/* Data stored in history table for a subpath */
struct SubpathEntry {
    int32_t cost;               // the cost of the prefix
    uint32_t active_threads;    // bitmask of threads exploring a path that contains this subpath
};

/* Compact format to store entry type (prefix/subpath and depth) */
struct HistoryEntryType {
private:
    int16_t data;
public:
    HistoryEntryType(bool is_subpath, size_t depth)
        // : data{static_cast<uint16_t>((depth << 1) & static_cast<size_t>(is_subpath))} {}
        : data{static_cast<int16_t>(depth * (is_subpath ? -1 : 1))} {}

    // inline bool is_subpath() { return static_cast<bool>(static_cast<unsigned>(data) & 1u); }
    // inline size_t depth() { return static_cast<size_t>(data) >> 1; }
    inline bool is_subpath() { return data < 0; }
    inline size_t depth() { return std::abs(data); }
};

/* Extra data on each entry used to choose which entry to replace when bucket is full */
struct HistoryEntryMetadata {
    HistoryEntryType type;      // contains is_subpath and depth in compact format
    uint16_t hits;              // times entry has been requested
};


#endif