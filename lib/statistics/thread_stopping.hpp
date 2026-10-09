#ifndef STATS_THREAD_STOPPING_H
#define STATS_THREAD_STOPPING_H

#include <vector>
#include <iostream>

struct stats_thread_stopping_thread {
    unsigned long long requests{0};            // requested a thread stop to another thread
    unsigned long long checks{0};               // handled a thread stop request from another thread
    unsigned long long compared{0};           // compared current bit vector to request bit vector
    unsigned long long success{0};              // thread successfully stopped

    stats_thread_stopping_thread &operator+=(const stats_thread_stopping_thread &t) {
        requests += t.requests;
        checks += t.checks;
        compared += t.compared;
        success += t.success;
        return *this;
    }
};

struct stats_thread_stopping {
    std::vector<stats_thread_stopping_thread> threads{};
    unsigned thread_count;

    stats_thread_stopping(unsigned thread_count)
        : thread_count{thread_count}
    {
        threads.reserve(thread_count);
        for (unsigned i = 0; i < thread_count; i++) {
            threads.push_back(stats_thread_stopping_thread());
        }
    }

    void print_results() {
        stats_thread_stopping_thread totals{};
        for (unsigned i = 0; i < thread_count; i++)
            totals += threads[i];
        
        std::cout << "THREAD STOPPING\n";
        std::cout << "Stops Requested:       " << totals.requests << "\n";
        std::cout << "Requests Checked:      " << totals.checks << "\n";
        std::cout << "Requests Compared:     " << totals.compared << "\n";
        std::cout << "Successfully Stopped:  " << totals.success << std::endl;
    }
};

#endif