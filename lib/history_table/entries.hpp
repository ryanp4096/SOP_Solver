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
    int32_t ordering;
};


/* Extra data on each entry used to choose which entry to replace when bucket is full */
struct HistoryEntryMetadata {
    int16_t depth;              // depth of prefix/subpath
    inline bool is_subpath() { return depth < 0; } // negative depth is used to indicate entry is a subpath
};


#endif