#ifndef STATS_HISTORY_TABLE_H
#define STATS_HISTORY_TABLE_H

#include <vector>
#include <iostream>
#include <iomanip>
#include "../history_table/history_table.hpp"

struct stats_history_table_thread {
    unsigned long long found{0};        // request was found in history table
    unsigned long long inserted{0};     // request was inserted into history table
    unsigned long long replaced{0};     // request was inserted by replacing another entry
    unsigned long long not_inserted{0}; // request was not inserted into history table

    stats_history_table_thread &operator+=(const stats_history_table_thread &t) {
        found += t.found;
        inserted += t.inserted;
        replaced += t.replaced;
        not_inserted += t.not_inserted;
        return *this;
    }
};

struct stats_history_table {
    std::vector<const stats_history_table_thread *> threads;
    unsigned thread_count;
    const HistoryTable &history_table;

    stats_history_table(unsigned thread_count, const HistoryTable &history_table)
        : threads{thread_count}, thread_count{thread_count}, history_table{history_table} {}

    void register_thread(unsigned thread_id, const stats_history_table_thread &thread) {
        threads[thread_id] = &thread;
    }

    void print_results() {
        stats_history_table_thread totals{};
        for (unsigned i = 0; i < thread_count; i++)
            if (threads[i] != nullptr) totals += *threads[i];
        
        auto percent = [](unsigned long long x, unsigned long long y) {
            return (y == 0) ? 0.0 : ((100.0 * x) / y);
        };

        unsigned long long insert_attempts = totals.inserted + totals.not_inserted + totals.replaced;
        unsigned long long requests = insert_attempts + totals.found;
        size_t buckets = history_table.get_bucket_count();
        size_t entries = buckets * slots_per_bucket;

        std::cout << "HISTORY TABLE\n";
        std::cout << "Requests:           " << std::setw(12) << std::right << requests << "\n";
        std::cout << "    Found:          " << std::setw(12) << std::right << totals.found << "  (" << percent(totals.found, requests) << "%)\n";
        std::cout << "    Not Found:      " << std::setw(12) << std::right << insert_attempts << "  (" << percent(insert_attempts, requests) << "%)\n";
        std::cout << "Insertions:         " << std::setw(12) << std::right << insert_attempts << "\n";
        std::cout << "    Inserted:       " << std::setw(12) << std::right << totals.inserted << "  (" << percent(totals.inserted, insert_attempts) << "%)\n";
        std::cout << "    Replaced:       " << std::setw(12) << std::right << totals.replaced << "  (" << percent(totals.replaced, insert_attempts) << "%)\n";
        std::cout << "    Not Inserted:   " << std::setw(12) << std::right << totals.not_inserted << "  (" << percent(totals.not_inserted, insert_attempts) << "%)\n";
        std::cout << "Capacity:\n";
        std::cout << "    Buckets:        " << buckets << "\n";
        std::cout << "    Entries:        " << totals.inserted << " / " << entries
                                            << "  (" << percent(totals.inserted, entries) << "%)\n";
        std::cout << "    Avg Per Bucket: " << (static_cast<double>(totals.inserted) / buckets) << " / " << slots_per_bucket << std::endl;
        std::cout << "    Memory Usage:   " << getPeakMemoryUsage() / 1048576 << " MB" << std::endl;

        std::vector<size_t> dist = history_table.bucket_distribution();
        std::cout << "Bucket Distribution:\n";
        for (size_t i = 0; i < dist.size(); i++) {
            size_t b = dist[i];
            std::cout << "    Size " << i << ": " << b << " buckets (" << percent(b, buckets) << "%)" << std::endl;
        }

        // auto dists = history_table.hit_distribution();
        // auto &prefix_dist = dists.first;
        // auto &subpath_dist = dists.second;
        // std::cout << "Hit Distribution:\n";
        // for (size_t block = 0; block < 16; block++) {
        //     std::cout << "    " << (1u << block) << " - " << ((1u << (block + 1)) - 1) << " hits: ";
        //     std::cout << prefix_dist[block] << " prefixes, " << subpath_dist[block] << " subpaths" << std::endl;
        // }
    }
};


#endif