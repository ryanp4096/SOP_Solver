#ifndef TRACE_H
#define TRACE_H

#include <fstream>
#include <string>
#include <iostream>

using namespace std;

/**
 * @brief Codes representing actions in a trace.
 */
enum TraceCode {
    TRACE_PRUNE_COST, // prune by path cost >= best cost
    TRACE_PRUNE_LEAF, // prune by being at the end of a path
    TRACE_HISTORY_NO_MATCH, // no matching entry in history table
    TRACE_HISTORY_PRUNE_COST, // matching history entry found with better prefix cost
    TRACE_HISTORY_PRUNE_LB, // matching history entry found, updated lower bound is worse than best cost
    TRACE_HISTORY_NO_PRUNE, // matching history entry found, but not pruned
    TRACE_PRUNE_LB, // prune by lower bound >= best cost
    TRACE_NO_PRUNE, // not pruned and added to ready list
    TRACE_CANCEL_THREAD_STOP, // cancelled by thread stopping
    TRACE_CANCEL_PRECHECK, // cancelled by enumeration precheck (lower bound check with potentially updated best cost)
    TRACE_ENUMERATE // recursively call enumerate from this node
};

/**
 * @brief Code representing the end of a list of items in a trace.
 */
#define TRACE_END_LIST UINT16_MAX

/**
 * @brief Version number included in trace file to differentiate between trace formats
 */
#define TRACE_VERSION_NUMBER 1

/**
 * @brief For writing a binary file representing a trace of the algorithm.
 */
class Trace {
private:
    bool opened;
    std::string path;
    std::ofstream file;
public:
    /**
     * @brief Create empty trace file.
     */
    Trace();

    /**
     * @brief Open trace file.
     */
    void open(std::string path);

    /**
     * @brief Close trace file.
     */
    void close();

    /**
     * @brief Check if trace file is opened
     */
    bool is_open();

    /**
     * @brief Write bytes to the trace file.
     * @param data The data to write.
     * @param bytes The number of bytes to write.
     */
    void write(unsigned long long data, size_t bytes);
};

#endif