#ifndef COMPARE_H
#define COMPARE_H

#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>

#include "reader.hpp"

class TraceCompare {
public:
    std::string path1;
    std::string path2;
    TraceReader r1;
    TraceReader r2;

    unsigned long long shared_nodes = 0;
    unsigned long long r1_only_times = 0;
    unsigned long long r1_only_nodes = 0;
    unsigned long long r2_only_times = 0;
    unsigned long long r2_only_nodes = 0;
    std::vector<unsigned long long> r1_only_times_by_type = std::vector<unsigned long long>(11);
    std::vector<unsigned long long> r1_only_nodes_by_type = std::vector<unsigned long long>(11);
    std::vector<unsigned long long> r2_only_times_by_type = std::vector<unsigned long long>(11);
    std::vector<unsigned long long> r2_only_nodes_by_type = std::vector<unsigned long long>(11);

    TraceCompare(const std::string& path1, const std::string& path2);
    void parse();
    void compare_node();
    void print_results();
};

#endif