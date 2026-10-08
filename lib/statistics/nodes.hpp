#ifndef STATS_NODES_H
#define STATS_NODES_H

#include <vector>
#include <iostream>

/* Enumerated nodes statistics for a specific depth */
struct stats_nodes_depth {
    unsigned long long enumerated{0};   // total enumerated nodes at this depth
    unsigned long long ready{0};        // total ready nodes at this depth
    unsigned long long recursive{0};    // total recursive nodes at this depth

    stats_nodes_depth &operator+=(const stats_nodes_depth &d) {
        enumerated += d.enumerated;
        ready += d.ready;
        recursive += d.recursive;
        return *this;
    }
};

/* Enumerated nodes statistics for a thread */
struct stats_nodes_thread {
    unsigned long long enumerated{0};               // total nodes that were checked
    unsigned long long ready{0};                    // total nodes that were added to local pool (not pruned)
    unsigned long long popped{0};                   // total nodes that were removed from local pool
    unsigned long long recursive{0};                // total nodes that recursively called (not thread stopped or pruned by precheck)

    unsigned long long prune_cost{0};               // pruned by current cost >= best cost
    unsigned long long prune_leaf{0};               // pruned by reaching a leaf node
    unsigned long long prune_prefix_history{0};     // pruned by being inferior to another prefix
    unsigned long long prune_lower_bound{0};        // pruned by lower bound >= best cost
    unsigned long long prune_subpath_history{0};    // pruned by a subpath being inferior to another subpath

    unsigned long long prune_thread_stop{0};        // pruned by thread stopping
    unsigned long long prune_precheck{0};           // pruned by re-checking lower bound before recursively calling
    std::vector<stats_nodes_depth> by_depth;

    stats_nodes_thread(unsigned instance_size) : by_depth{instance_size + 1} {}

    stats_nodes_thread &operator+=(const stats_nodes_thread &t) {
        enumerated += t.enumerated;
        ready += t.ready;
        popped += t.popped;
        recursive += t.recursive;
        prune_cost += t.prune_cost;
        prune_leaf += t.prune_leaf;
        prune_prefix_history += t.prune_prefix_history;
        prune_lower_bound += t.prune_lower_bound;
        prune_subpath_history += t.prune_subpath_history;
        prune_thread_stop += t.prune_thread_stop;
        prune_precheck += t.prune_precheck;
        for (std::size_t i = 0; i < by_depth.size(); i++)
            by_depth[i] += t.by_depth[i];
        return *this;
    }
};

/* Enumerated nodes statistics for a run */
struct stats_nodes {
    std::vector<stats_nodes_thread> threads{};
    unsigned thread_count;
    unsigned instance_size;

    stats_nodes(unsigned thread_count, unsigned instance_size)
        : thread_count{thread_count}, instance_size{instance_size}
    {
        threads.reserve(thread_count);
        for (unsigned i = 0; i < thread_count; i++) {
            threads.push_back(stats_nodes_thread(instance_size));
        }
    }

    void print_results() {
        stats_nodes_thread totals{instance_size};
        for (unsigned i = 0; i < thread_count; i++)
            totals += threads[i];
        
        std::cout << "Enumerated Nodes: " << totals.enumerated << "\n";
        unsigned long long remaining = totals.enumerated;

        std::cout << "    Pruned by cost:             " << totals.prune_cost << " / " << remaining << "\n";
        remaining -= totals.prune_cost;

        std::cout << "    Pruned by leaf:             " << totals.prune_leaf << " / " << remaining << "\n";
        remaining -= totals.prune_leaf;

        std::cout << "    Pruned by prefix history:   " << totals.prune_prefix_history << " / " << remaining << "\n";
        remaining -= totals.prune_prefix_history;

        std::cout << "    Pruned by lower bound:      " << totals.prune_lower_bound << " / " << remaining << "\n";
        remaining -= totals.prune_lower_bound;

        std::cout << "    Pruned by subpath history:  " << totals.prune_subpath_history << " / " << remaining << "\n";
        remaining -= totals.prune_subpath_history;

        std::cout << "Ready Nodes: " << totals.ready << "\n";
        if (remaining != totals.ready)
            std::cout << "!! ERROR: Miscount: " << remaining << " remaining nodes != " << totals.ready << " ready nodes\n";

        if (totals.ready != totals.popped)
            std::cout << "Popped Nodes: " << totals.popped << "\n";
        remaining = totals.popped;

        std::cout << "    Pruned by thread stop:      " << totals.prune_thread_stop << " / " << remaining << "\n";
        remaining -= totals.prune_thread_stop;

        std::cout << "    Pruned by precheck:         " << totals.prune_precheck << " / " << remaining << "\n";
        remaining -= totals.prune_precheck;

        if (remaining != totals.recursive)
            std::cout << "!! ERROR: Miscount: " << remaining << " remaining nodes != " << totals.recursive << " recursive nodes\n";

        std::cout << "Recursive Nodes: " << totals.recursive << std::endl;
    }
};

#endif