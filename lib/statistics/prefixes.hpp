#ifndef STATS_PREFIXES_H
#define STATS_PREFIXES_H

#include <vector>
#include <cstddef>
#include <iostream>
#include <iomanip>

struct stats_prefixes_depth {
    unsigned long long pruned{0};       // inferior to the existing prefix
    unsigned long long updated{0};      // better than the existing prefix
    unsigned long long inserted{0};     // not already in history table, inserted
    unsigned long long not_inserted{0}; // not already in history table, did not insert

    stats_prefixes_depth &operator+=(const stats_prefixes_depth &d) {
        pruned += d.pruned;
        updated += d.updated;
        inserted += d.inserted;
        not_inserted += d.not_inserted;
        return *this;
    }
};

struct stats_prefixes_thread {
    unsigned long long checks{0};               // total prefixes checked

    unsigned long long found{0};                // prefix found in history table
    unsigned long long pruned{0};               // inferior to the existing prefix
    unsigned long long updated{0};              // better than the existing prefix, update entry, stop inferior thread

    unsigned long long not_found{0};            // prefix not found in history table
    unsigned long long inserted{0};             // not already in history table, inserted
    unsigned long long not_inserted{0};         // not already in history table, did not insert (out of memory)

    std::vector<stats_prefixes_depth> by_depth;

    stats_prefixes_thread(unsigned instance_size) : by_depth(instance_size + 1) {}

    stats_prefixes_thread &operator+=(const stats_prefixes_thread &t) {
        checks += t.checks;
        found += t.found;
        pruned += t.pruned;
        updated += t.updated;
        not_found += t.not_found;
        inserted += t.inserted;
        not_inserted += t.not_inserted;
        for (size_t i = 0; i < by_depth.size(); i++)
            by_depth[i] += t.by_depth[i];
        return *this;
    }
};

struct stats_prefixes {
    std::vector<const stats_prefixes_thread *> threads;
    unsigned thread_count;
    unsigned instance_size;

    stats_prefixes(unsigned thread_count, unsigned instance_size)
        : threads{thread_count}, thread_count{thread_count}, instance_size{instance_size} {}

    void register_thread(unsigned thread_id, const stats_prefixes_thread &thread) {
        threads[thread_id] = &thread;
    }

    void print_results() {
        stats_prefixes_thread totals{instance_size};
        for (unsigned i = 0; i < thread_count; i++)
            if (threads[i] != nullptr) totals += *threads[i];

        auto percent = [](unsigned long long x, unsigned long long y) {
            return (y == 0) ? 0.0 : ((100.0 * x) / y);
        };

        std::cout << std::fixed << std::setprecision(3);

        std::cout << "PREFIX HISTORY\n";
        std::cout << "Checks:           " << std::setw(12) << std::right << totals.checks << "\n";
        std::cout << "Found:            " << std::setw(12) << std::right << totals.found << "  (" << percent(totals.found, totals.checks) << "%)\n";
        std::cout << "    Pruned:       " << std::setw(12) << std::right << totals.pruned << "  (" << percent(totals.pruned, totals.checks) << "%)\n";
        std::cout << "    Updated:      " << std::setw(12) << std::right << totals.updated << "  (" << percent(totals.updated, totals.checks) << "%)\n";
        std::cout << "Not Found:        " << std::setw(12) << std::right << totals.not_found << "  (" << percent(totals.not_found, totals.checks) << "%)\n";
        std::cout << "    Inserted:     " << std::setw(12) << std::right << totals.inserted << "  (" << percent(totals.inserted, totals.checks) << "%)\n";
        std::cout << "    Not Inserted: " << std::setw(12) << std::right << totals.not_inserted << "  (" << percent(totals.not_inserted, totals.checks) << "%)\n";

        unsigned first;
        for (first = 0; first <= instance_size; first++) {
            const auto &d = totals.by_depth[first];
            unsigned long long checks = d.pruned + d.updated + d.inserted + d.not_inserted;
            if (checks != 0) break;
        }
        unsigned last;
        for (last = instance_size; last >= 0; last--) {
            const auto &d = totals.by_depth[last];
            unsigned long long checks = d.pruned + d.updated + d.inserted + d.not_inserted;
            if (checks != 0) break;
        }

        std::cout << std::setw(9) << std::left << "By Depth:" << " ";
        std::cout << std::setw(12) << std::right << "Checks" << " ";
        std::cout << std::setw(8) << std::right << "Pruned" << " ";
        std::cout << std::setw(8) << std::right << "Updated" << " ";
        std::cout << std::setw(8) << std::right << "NotFound" << " ";
        std::cout << std::setw(12) << std::right << "Entries" << " ";
        std::cout << std::setw(8) << std::right << " Prunes/Entry" << "\n";

        for (unsigned i = first; i <= last; i++) {
            const auto &d = totals.by_depth[i];
            unsigned long long found = d.pruned + d.updated;
            unsigned long long not_found = d.inserted + d.not_inserted;
            unsigned long long checks = found + not_found;
            std::cout << std::setw(9) << std::left << "    [" + std::to_string(i) + "]" << " ";
            std::cout << std::setw(12) << std::right << checks << " ";
            std::cout << std::setw(7) << std::right << percent(d.pruned, checks) << "% ";
            std::cout << std::setw(7) << std::right << percent(d.updated, checks) << "% ";
            std::cout << std::setw(7) << std::right << percent(not_found, checks) << "% ";
            std::cout << std::setw(12) << std::right << d.inserted << " ";
            std::cout << std::setw(7) << std::right << ((checks == 0) ? 0.0 : (static_cast<double>(d.pruned) / d.inserted)) << std::endl;
        }
    }
};

#endif