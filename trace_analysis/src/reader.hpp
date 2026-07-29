#ifndef READER_H
#define READER_H

#include <fstream>
#include <filesystem>
#include <string>
#include <vector>
#include <cstdint>
#include <iostream>

#define TRACE_END_LIST UINT16_MAX
#define MAX_ITERATIONS UINT16_MAX

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

enum TraceMatchCode {
    MATCH_NOT_FOUND, // didn't match prefix or subpath
    MATCH_PREFIX, // matched lkh prefix
    MATCH_NOT_AVAILABLE, // lkh solution not processed
    MATCH_SUBPATH // matched lkh subpath
};

enum TraceDetailLevel {
    DETAIL_NORMAL,
    DETAIL_COMPACT
};

struct NodeCheck {
    int node;
    int depth;
    int path_cost;
    int best_cost;
    TraceMatchCode lkh_match_type;
    int lkh_match_cost = -1;
    TraceCode action;
    bool history_match;
    int history_lower_bound = -1;
    int history_path_cost = -1;
    bool history_is_best_suffix;
    int lower_bound = -1;
};

struct NodeCall {
    int node;
    TraceCode action;
    int depth;
    bool thread_stop_prefix_key_matched;
    int enumeration_precheck_best_cost = -1;
    int enumeration_precheck_lower_bound = -1;
};

class TraceReader {
public:
    std::string path;
    std::ifstream file;
    std::vector<int> current_path;
    unsigned long long version_number;
    int thread_id;
    TraceDetailLevel detail_level;
    int instance_size_shift;
    int instance_size;

    unsigned long long enumerated_nodes = 0;
    unsigned long long ready_nodes = 0;
    unsigned long long recursive_nodes = 0;
    
    TraceReader(const std::string& path);
    unsigned long long next(size_t bytes);
    unsigned long long next_detail(size_t bytes);
    int read_node();
    void read_header();
    void parse();
    void parse_initial_state();
    unsigned long long parse_node();
    bool parse_node_check(NodeCheck *node);
    bool parse_node_call(NodeCall *node);
};

#endif