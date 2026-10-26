/**
 * CPU TIMER
 * The CPU timer measures the number of clock cycles spent executing various code segments.
 * The measurement adds overhead, so it is disabled by default.
 * To enable, you must compile with `make ENABLE_CPU_TIMER=1`
 */

#ifndef STATS_CPU_TIMER_H
#define STATS_CPU_TIMER_H

#include <vector>
#include <iostream>
#include <iomanip>

#ifdef ENABLE_CPU_TIMER
#include <x86intrin.h>
#endif

struct stats_cpu_timer_segment {
    std::string name;
    const bool enabled;
    unsigned long long total_count{0}; // total number of times this code segment has been timed
    unsigned long long total_time{0};  // total cpu cycles this segment has executed for
    unsigned long long start_time{0};  // cpu cycle count when start was called for this segment

    stats_cpu_timer_segment(const std::string &name, bool enabled = true) : name{name}, enabled{enabled} {}

    inline void start_timer() {
        #ifdef ENABLE_CPU_TIMER
        if (!enabled) return;
        unsigned temp;
        _mm_mfence();
        start_time = __rdtscp(&temp);
        _mm_mfence();
        #endif
    }

    inline void stop_timer() {
        #ifdef ENABLE_CPU_TIMER
        if (!enabled) return;
        unsigned temp;
        _mm_mfence();
        unsigned long long end_time = __rdtscp(&temp);
        _mm_mfence();
        total_count++;
        total_time += end_time - start_time;
        #endif
    }

    stats_cpu_timer_segment &operator+=(const stats_cpu_timer_segment &s) {
        total_count += s.total_count;
        total_time += s.total_time;
        return *this;
    }
};

struct stats_cpu_timer_thread {
    stats_cpu_timer_segment enumerate_node_start{"Node Start", false};
    stats_cpu_timer_segment enumerate_prefix_history{"Prefix History", false};
    stats_cpu_timer_segment enumerate_subpath_history{"Subpath History"};
    stats_cpu_timer_segment enumerate_node_end{"Node End", false};
    stats_cpu_timer_segment enumerate_sort_pool{"Sort Pool", false};
    stats_cpu_timer_segment enumerate_thread_stop{"Thread Stopping", false};
    stats_cpu_timer_segment enumerate_recursive_start{"Recursive Start", false};
    stats_cpu_timer_segment enumerate_recursive_end{"Recursive End", false};

    stats_cpu_timer_segment lower_bound{"Lower Bound", false};

    stats_cpu_timer_segment thread_stop_request{"Thread Stop Request", false};
    stats_cpu_timer_segment thread_stop_check{"Thread Stop Check", false};

    stats_cpu_timer_segment prefix_history_lookup{"Prefix History Lookup", false};
    stats_cpu_timer_segment prefix_history_entry{"Prefix History Entry", false};
    
    stats_cpu_timer_segment subpath_history_lookup{"Subpath History Lookup"};
    stats_cpu_timer_segment subpath_history_entry{"Subpath History Entry"};
    stats_cpu_timer_segment subpath_history_hash{"Subpath History Hash"};
    stats_cpu_timer_segment subpath_history_search{"Subpath History Search"};
    stats_cpu_timer_segment subpath_history_insert{"Subpath History Insert"};

    std::vector<stats_cpu_timer_segment *> segments() {
        return {
            &enumerate_node_start, &enumerate_prefix_history, &enumerate_subpath_history,
            &enumerate_node_end, &enumerate_sort_pool, &enumerate_thread_stop,
            &enumerate_recursive_start, &enumerate_recursive_end,
            &lower_bound, &thread_stop_request, &thread_stop_check,
            &prefix_history_lookup, &prefix_history_entry, &subpath_history_lookup,
            &subpath_history_entry, &subpath_history_hash, &subpath_history_search, &subpath_history_insert
        };
    }

    stats_cpu_timer_thread &operator+=(const stats_cpu_timer_thread &t) {
        enumerate_node_start += t.enumerate_node_start;
        enumerate_prefix_history += t.enumerate_prefix_history;
        enumerate_subpath_history += t.enumerate_subpath_history;
        enumerate_node_end += t.enumerate_node_end;
        enumerate_sort_pool += t.enumerate_sort_pool;
        enumerate_thread_stop += t.enumerate_thread_stop;
        enumerate_recursive_start += t.enumerate_recursive_start;
        enumerate_recursive_end += t.enumerate_recursive_end;

        lower_bound += t.lower_bound;

        thread_stop_request += t.thread_stop_request;
        thread_stop_check += t.thread_stop_check;

        prefix_history_lookup += t.prefix_history_lookup;
        prefix_history_entry += t.prefix_history_entry;
    
        subpath_history_lookup += t.subpath_history_lookup;
        subpath_history_entry += t.subpath_history_entry;
        subpath_history_hash += t.subpath_history_hash;
        subpath_history_search += t.subpath_history_search;
        subpath_history_insert += t.subpath_history_insert;

        return *this;
    }
};

struct stats_cpu_timer {
    std::vector<const stats_cpu_timer_thread *> threads;
    unsigned thread_count;

    stats_cpu_timer(unsigned thread_count)
        : threads{thread_count}, thread_count{thread_count} {}

    void register_thread(unsigned thread_id, const stats_cpu_timer_thread &thread) {
        threads[thread_id] = &thread;
    }

    void print_results() {
        #ifdef ENABLE_CPU_TIMER
        stats_cpu_timer_thread totals{};

        for (unsigned i = 0; i < thread_count; i++)
            if (threads[i] != nullptr) totals += *threads[i];
        
        std::cout << "CPU TIMER\n";

        std::cout << "    " << std::setw(30) << std::left << "CODE SEGMENT" << "  ";
        std::cout << std::setw(18) << std::right << "AVERAGE TIME" << "  ";
        std::cout << std::setw(18) << std::right << "COUNT" << "  ";
        std::cout << std::setw(18) << std::right << "TOTAL TIME" << std::endl;

        for (stats_cpu_timer_segment *segment : totals.segments()) {
            if (!segment->enabled) continue;
            unsigned long long average = 0;
            if (segment->total_count > 0) average = segment->total_time / segment->total_count;

            std::cout << "    " << std::setw(30) << std::left << segment->name << "  ";
            std::cout << std::setw(18) << std::right << average << "  ";
            std::cout << std::setw(18) << std::right << segment->total_count << "  ";
            std::cout << std::setw(18) << std::right << segment->total_time << std::endl;
        }
        #endif
    }
};

#endif