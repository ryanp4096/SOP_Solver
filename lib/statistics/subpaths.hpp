#ifndef STATS_SUBPATHS_H
#define STATS_SUBPATHS_H

#include <vector>
#include <iostream>

struct stats_subpaths_depth {
    unsigned long long checks{0};           // total times a subpath at this depth was checked
    unsigned long long checks_pruned{0};    // there was an existing superior matching subpath, leading to pruning
    unsigned long long checks_equal{0};     // there was an existing matching subpath with equal cost (can't prune)
    unsigned long long checks_improved{0};  // there was an existing matching subpath, but it was inferior to the current subpath. this could be a thread stop once that is implemented
    unsigned long long checks_no_match{0};  // there was not an existing matching subpath

    stats_subpaths_depth &operator+=(const stats_subpaths_depth &d) {
        checks += d.checks;
        checks_pruned += d.checks_pruned;
        checks_equal += d.checks_equal;
        checks_improved += d.checks_improved;
        checks_no_match += d.checks_no_match;
        return *this;
    }
};

struct stats_subpaths_thread {
    unsigned long long checks{0};           // total times an individual subpath was checked

    unsigned long long found{0};            // matching subpath found in history table
    unsigned long long pruned{0};           // there was an existing superior matching subpath, leading to pruning
    unsigned long long equal{0};            // there was an existing matching subpath with equal cost (can't prune)
    unsigned long long updated{0};          // this subpath was better than the existing subpath. this could be a thread stop once that is implemented

    unsigned long long not_found{0};        // subpath not found in history table
    unsigned long long inserted{0};         // not in history table, inserted
    unsigned long long not_inserted{0};     // not in history table, not inserted (out of memory or lkh subpaths only)

    std::vector<stats_subpaths_depth> by_depth;

    stats_subpaths_thread(unsigned instance_size) : by_depth(instance_size + 1) {}

    stats_subpaths_thread &operator+=(const stats_subpaths_thread &t) {
        checks += t.checks;
        found += t.found;
        pruned += t.pruned;
        equal += t.equal;
        updated += t.updated;
        not_found += t.not_found;
        inserted += t.inserted;
        not_inserted += t.not_inserted;
        for (std::size_t i = 0; i < by_depth.size(); i++)
            by_depth[i] += t.by_depth[i];
        return *this;
    }
};

struct stats_subpaths {
    std::vector<stats_subpaths_thread> threads{};
    unsigned thread_count;
    unsigned instance_size;

    stats_subpaths(unsigned thread_count, unsigned instance_size)
        : thread_count{thread_count}, instance_size{instance_size}
    {
        threads.reserve(thread_count);
        for (unsigned i = 0; i < thread_count; i++) {
            threads.push_back(stats_subpaths_thread(instance_size));
        }
    }

    void print_results() {
        stats_subpaths_thread totals{instance_size};
        for (unsigned i = 0; i < thread_count; i++)
            totals += threads[i];

        std::cout << "SUBPATH HISTORY\n";
        std::cout << "Checks:            " << totals.checks << "\n";
        std::cout << "Found:             " << totals.found << "\n";
        std::cout << "    Pruned:        " << totals.pruned << "\n";
        std::cout << "    Equal Cost:    " << totals.equal << "\n";
        std::cout << "    Updated:       " << totals.updated << "\n";
        std::cout << "Not Found:         " << totals.not_found << "\n";
        std::cout << "    Inserted:      " << totals.inserted << "\n";
        std::cout << "    Not Inserted:  " << totals.not_inserted << std::endl;

        // std::cout << "Checks By Depth: " << std::endl;

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
        // std::cout << "Depth\tChecks\tPruned\tEqual\tImproved\tNoMatch" << std::endl;
        // for (unsigned j = 0; j < instance_size + 1; j++) {
        //     std::cout << j << "\t" << totals.by_depth[j].checks << "\t";
        //     std::cout << totals.by_depth[j].checks_pruned << "\t" << totals.by_depth[j].checks_equal << "\t";
        //     std::cout << totals.by_depth[j].checks_improved << "\t" << totals.by_depth[j].checks_no_match << std::endl;
        // }
    }
};

#endif