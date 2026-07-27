#include "reader.hpp"

TraceReader::TraceReader(const std::string& path)
    : path(path), file(path, std::ios::binary)
{

}

unsigned long long TraceReader::next(size_t bytes) {
    switch (bytes) {
        case 1:
            uint8_t x8;
            file.read(reinterpret_cast<char *>(&x8), 1);
            return static_cast<unsigned long long>(x8);
        
        case 2:
            uint16_t x16;
            file.read(reinterpret_cast<char *>(&x16), 2);
            return static_cast<unsigned long long>(x16);

        case 4:
            uint32_t x32;
            file.read(reinterpret_cast<char *>(&x32), 4);
            return static_cast<unsigned long long>(x32);
        
        case 8:
            uint64_t x64;
            file.read(reinterpret_cast<char *>(&x64), 8);
            return static_cast<unsigned long long>(x64);

        default:
            return 0;
    }
}

void TraceReader::read_header() {
    version_number = next(4);
    thread_id = next(1);
}

void TraceReader::parse() {
    read_header();
    while (!file.eof()) {
        parse_initial_state();
        parse_node();
        file.peek();
    }
}

void TraceReader::parse_initial_state() {
    current_path.clear();
    for (int i = 0; i < MAX_ITERATIONS; i++) {
        int n = next(2);
        if (n == TRACE_END_LIST) return;
        current_path.push_back(n);
    }
}

unsigned long long TraceReader::parse_node() {
    unsigned long long en = 0;
    recursive_nodes++;
    for (int i = 0; i < MAX_ITERATIONS; i++) {
        NodeCheck node;
        bool x = parse_node_check(&node);
        if (!x) break;
        enumerated_nodes++;
        en++;
    }
    for (int i = 0; i < MAX_ITERATIONS; i++) {
        NodeCall node;
        bool x = parse_node_call(&node);
        if (!x) break;
        ready_nodes++;
        if (node.action == TRACE_ENUMERATE) {
            current_path.push_back(node.node);
            en += parse_node();
            current_path.pop_back();
        }
    }
    return en;
}

bool TraceReader::parse_node_check(NodeCheck *node) {
    NodeCheck temp;
    if (node == NULL) node = &temp;

    int n = next(2);
    if (n == TRACE_END_LIST) return false;

    node->node = n;
    node->depth = current_path.size() + 1;
    node->path_cost = next(4);
    node->best_cost = next(4);
    node->lkh_match = next(1) == 1;
    if (node->lkh_match) {
        node->lkh_match_cost = next(4);
    }

    TraceCode action = static_cast<TraceCode>(next(1));
    if (action != TRACE_PRUNE_COST && action != TRACE_PRUNE_LEAF) {
        if (action == TRACE_HISTORY_NO_MATCH) {
            node->history_match = false;
            node->lower_bound = next(4);
            action = static_cast<TraceCode>(next(1));

        } else {
            node->history_match = true;
            node->history_lower_bound = next(4);
            node->history_path_cost = next(4);
            node->history_is_best_suffix = next(1);

            if (node->history_is_best_suffix && action != TRACE_HISTORY_PRUNE_COST)
                node->lower_bound = node->history_lower_bound - (node->history_path_cost - node->path_cost);

            if (!node->history_is_best_suffix) {
                if (action != TRACE_HISTORY_PRUNE_COST)
                    node->lower_bound = node->history_lower_bound;
                node->history_lower_bound = -1;
            }

            if (action == TRACE_HISTORY_NO_PRUNE)
                action = static_cast<TraceCode>(next(1));
        }
    }
    node->action = action;

    return true;
}

bool TraceReader::parse_node_call(NodeCall *node) {
    NodeCall temp;
    if (node == NULL) node = &temp;

    int n = next(2);
    if (n == TRACE_END_LIST) return false;

    node->node = n;
    node->depth = current_path.size() + 1;

    TraceCode action = static_cast<TraceCode>(next(1));

    if (action == TRACE_CANCEL_THREAD_STOP) {
        node->thread_stop_prefix_key_matched = next(1);
    }

    if (action == TRACE_CANCEL_PRECHECK) {
        node->enumeration_precheck_best_cost = next(4);
        node->enumeration_precheck_lower_bound = next(4);
    }

    node->action = action;

    return true;
}