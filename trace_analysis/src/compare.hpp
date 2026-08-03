#ifndef COMPARE_H
#define COMPARE_H

#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>

#include "reader.hpp"

struct ResultsTreeNode {
    int node;
    bool different; // whether this node is different between R1 and R2
    std::string r1; // if different: r1's behavior
    std::string r2; // if different: r2's behavior
    std::vector<ResultsTreeNode> children; // if same: any child paths with differences
};

class TraceCompare {
public:
    std::string path1;
    std::string path2;
    TraceReader r1;
    TraceReader r2;

    std::vector<ResultsTreeNode> node_results;
    int print_indent_count;

    unsigned long long shared_nodes = 0;
    unsigned long long r1_only_times = 0;
    unsigned long long r1_only_nodes = 0;
    unsigned long long r2_only_times = 0;
    unsigned long long r2_only_nodes = 0;
    std::vector<unsigned long long> r1_only_times_by_type = std::vector<unsigned long long>(11);
    std::vector<unsigned long long> r1_only_nodes_by_type = std::vector<unsigned long long>(11);
    std::vector<unsigned long long> r2_only_times_by_type = std::vector<unsigned long long>(11);
    std::vector<unsigned long long> r2_only_nodes_by_type = std::vector<unsigned long long>(11);
    unsigned long long r1_only_error_times = 0;
    unsigned long long r1_only_error_nodes = 0;
    unsigned long long r2_only_error_times = 0;
    unsigned long long r2_only_error_nodes = 0;

    TraceCompare(const std::string& path1, const std::string& path2);
    void parse();
    void compare_node();
    void print_results();
    void print_tree();
    void print_tree_node(ResultsTreeNode& node, int indent);
};

#endif