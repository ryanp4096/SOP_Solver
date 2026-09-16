#ifndef SUBPATH_DIAGNOSTICS_H
#define SUBPATH_DIAGNOSTICS_H

#include <vector>
#include <iostream>

struct subpath_depth_data {
    unsigned long long checks{0}; // total times a subpath at this depth was checked
    unsigned long long checks_pruned{0}; // there was an existing superior matching subpath, leading to pruning
    unsigned long long checks_equal{0}; // there was an existing matching subpath with equal cost (can't prune)
    unsigned long long checks_improved{0}; // there was an existing matching subpath, but it was inferior to the current subpath. this could be a thread stop once that is implemented
    unsigned long long checks_no_match{0}; // there was not an existing matching subpath
};

struct subpath_thread_data {
    unsigned long long nodes{0}; // total nodes that reached subpath checking phase
    unsigned long long nodes_pruned{0}; // total nodes that were pruned by subpath checking
    unsigned long long nodes_not_pruned{0}; // total nodes that were not pruned by subpath checking
    
    unsigned long long checks{0};
    unsigned long long checks_pruned{0};
    unsigned long long checks_equal{0};
    unsigned long long checks_improved{0};
    unsigned long long checks_no_match{0};

    std::vector<subpath_depth_data> by_depth;

    subpath_thread_data(unsigned instance_size) : by_depth(instance_size + 1) {}
};

struct subpath_data {
    std::vector<subpath_thread_data> threads{};
    unsigned thread_count;
    unsigned instance_size;

    subpath_data(unsigned thread_count, unsigned instance_size)
        : thread_count{thread_count}, instance_size{instance_size}
    {
        threads.reserve(thread_count);
        for (unsigned i = 0; i < thread_count; i++) {
            threads.push_back(subpath_thread_data(instance_size));
        }
    }

    void print_results() {
        subpath_thread_data totals{instance_size};

        for (unsigned i = 0; i < thread_count; i++) {
            totals.nodes += threads[i].nodes;
            totals.nodes_pruned += threads[i].nodes_pruned;
            totals.nodes_not_pruned += threads[i].nodes_not_pruned;
            totals.checks += threads[i].checks;
            totals.checks_pruned += threads[i].checks_pruned;
            totals.checks_equal += threads[i].checks_equal;
            totals.checks_improved += threads[i].checks_improved;
            totals.checks_no_match += threads[i].checks_no_match;
            
            for (unsigned j = 0; j < instance_size + 1; j++) {
                totals.by_depth[j].checks += threads[i].by_depth[j].checks;
                totals.by_depth[j].checks_pruned += threads[i].by_depth[j].checks_pruned;
                totals.by_depth[j].checks_equal += threads[i].by_depth[j].checks_equal;
                totals.by_depth[j].checks_improved += threads[i].by_depth[j].checks_improved;
                totals.by_depth[j].checks_no_match += threads[i].by_depth[j].checks_no_match;
            }
        }

        std::cout << "======== SUBPATH DIAGNOSTICS ========" << std::endl;
        std::cout << "Total Nodes Processed: " << totals.nodes << std::endl;
        std::cout << "   Pruned: " << totals.nodes_pruned << std::endl;
        std::cout << "   Not Pruned: " << totals.nodes_not_pruned << std::endl << std::endl;

        std::cout << "Total Checks: " << totals.checks << std::endl;
        std::cout << "   Pruned: " << totals.checks_pruned << std::endl;
        std::cout << "   Not Pruned (Equal Cost): " << totals.checks_equal << std::endl;
        std::cout << "   Not Pruned (Improvement): " << totals.checks_improved << std::endl;
        std::cout << "   Not Pruned (No Match): " << totals.checks_no_match << std::endl << std::endl;

        std::cout << "Checks By Depth: " << std::endl;

        // Readable
        // for (unsigned j = 0; j < instance_size + 1; j++) {
        //     std::cout << "   Depth " << j << ": ";
        //     std::cout << totals.by_depth[j].checks << " (";
        //     std::cout << "P: " << totals.by_depth[j].checks_pruned;
        //     std::cout << " Eq: " << totals.by_depth[j].checks_equal;
        //     std::cout << " Imp: " << totals.by_depth[j].checks_improved;
        //     std::cout << " NM: " << totals.by_depth[j].checks_no_match;
        //     std::cout << ")" << std::endl;
        // }

        // Table (for spreadsheets)
        std::cout << "Depth\tChecks\tPruned\tEqual\tImproved\tNoMatch" << std::endl;
        for (unsigned j = 0; j < instance_size + 1; j++) {
            std::cout << j << "\t" << totals.by_depth[j].checks << "\t";
            std::cout << totals.by_depth[j].checks_pruned << "\t" << totals.by_depth[j].checks_equal << "\t";
            std::cout << totals.by_depth[j].checks_improved << "\t" << totals.by_depth[j].checks_no_match << std::endl;
        }

        std::cout << std::endl << "=====================================" << std::endl;
    }
};

#endif