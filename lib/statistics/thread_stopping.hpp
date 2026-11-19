#ifndef STATS_THREAD_STOPPING_H
#define STATS_THREAD_STOPPING_H

#include <vector>
#include <iostream>

struct stats_thread_stopping_thread {
    unsigned long long reqs_attempted{0};       // times a thread wanted to push a stop request to the queue
    unsigned long long reqs_pushed{0};          // stop requests successfully pushed to the queue

    unsigned long long queue_checks{0};
    unsigned long long reqs_processed{0};       // stop requests processed by the target thread
    unsigned long long reqs_compared{0};        // stop requests that matched the current path and could be compared
    unsigned long long reqs_pruned{0};          // stop requests that were successfully pruned

    unsigned long long rechecks_processed{0};   // times a thread re-checked the history table to see if path was dominated 
    unsigned long long rechecks_pruned{0};      // times pruned by re-checking history table to see if path was dominated

    stats_thread_stopping_thread &operator+=(const stats_thread_stopping_thread &t) {
        reqs_attempted += t.reqs_attempted;
        reqs_pushed += t.reqs_pushed;
        queue_checks += t.queue_checks;
        reqs_processed += t.reqs_processed;
        reqs_compared += t.reqs_compared;
        reqs_pruned += t.reqs_pruned;
        rechecks_processed += t.rechecks_processed;
        rechecks_pruned += t.rechecks_pruned;
        return *this;
    }
};

struct stats_thread_stopping {
    std::vector<const stats_thread_stopping_thread *> threads;
    unsigned thread_count;

    stats_thread_stopping(unsigned thread_count)
        : threads{thread_count}, thread_count{thread_count} {}

    void register_thread(unsigned thread_id, const stats_thread_stopping_thread &thread) {
        threads[thread_id] = &thread;
    }

    void print_results() {
        stats_thread_stopping_thread totals{};
        for (unsigned i = 0; i < thread_count; i++)
            if (threads[i] != nullptr) totals += *threads[i];
        
        std::cout << "THREAD STOPPING\n";
        std::cout << "Requests:\n";
        std::cout << "    Attempted:         " << totals.reqs_attempted << "\n";
        std::cout << "    Pushed:            " << totals.reqs_pushed << "\n";
        std::cout << "Queue Checks:          " << totals.queue_checks << "\n";
        std::cout << "    Processed:         " << totals.reqs_processed << "\n";
        std::cout << "    Compared:          " << totals.reqs_compared << "\n";
        std::cout << "    Pruned:            " << totals.reqs_pruned << "\n";
        std::cout << "Rechecks:              " << totals.rechecks_processed << "\n";
        std::cout << "    Pruned:            " << totals.rechecks_pruned << std::endl;
    }
};

#endif