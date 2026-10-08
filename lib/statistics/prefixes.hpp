#ifndef STATS_PREFIXES_H
#define STATS_PREFIXES_H

#include <vector>
#include <iostream>

struct stats_prefixes_thread {
    unsigned long long checks{0};               // total prefixes checked

    unsigned long long found{0};                // prefix found in history table
    unsigned long long pruned{0};               // inferior to the existing prefix
    unsigned long long updated{0};              // better than the existing prefix, update entry, stop inferior thread

    unsigned long long not_found{0};            // prefix not found in history table
    unsigned long long inserted{0};             // not already in history table, inserted
    unsigned long long not_inserted{0};         // not already in history table, did not insert (out of memory)

    stats_prefixes_thread &operator+=(const stats_prefixes_thread &t) {
        checks += t.checks;
        found += t.found;
        pruned += t.pruned;
        updated += t.updated;
        not_found += t.not_found;
        inserted += t.inserted;
        not_inserted += t.not_inserted;
        return *this;
    }
};

struct stats_prefixes {
    std::vector<stats_prefixes_thread> threads{};
    unsigned thread_count;

    stats_prefixes(unsigned thread_count)
        : thread_count{thread_count}
    {
        threads.reserve(thread_count);
        for (unsigned i = 0; i < thread_count; i++) {
            threads.push_back(stats_prefixes_thread());
        }
    }

    void print_results() {
        stats_prefixes_thread totals{};
        for (unsigned i = 0; i < thread_count; i++)
            totals += threads[i];
        
        std::cout << "PREFIX HISTORY\n";
        std::cout << "Checks:            " << totals.checks << "\n";
        std::cout << "Found:             " << totals.found << "\n";
        std::cout << "    Pruned:        " << totals.pruned << "\n";
        std::cout << "    Updated:       " << totals.updated << "\n";
        std::cout << "Not Found:         " << totals.not_found << "\n";
        std::cout << "    Inserted:      " << totals.inserted << "\n";
        std::cout << "    Not Inserted:  " << totals.not_inserted << std::endl;
    }
};

#endif