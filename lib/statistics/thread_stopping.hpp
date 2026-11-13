#ifndef STATS_THREAD_STOPPING_H
#define STATS_THREAD_STOPPING_H

#include <vector>
#include <iostream>

struct stats_thread_stopping_thread {
    unsigned long long requests{0};            // requested a thread stop to another thread
    unsigned long long checks{0};               // handled a thread stop request from another thread
    unsigned long long compared{0};           // compared current bit vector to request bit vector
    unsigned long long success{0};              // thread successfully stopped
    unsigned long long cancel_schedule{0};      // a scheduled node in local pool was dominated before being recursively called

    stats_thread_stopping_thread &operator+=(const stats_thread_stopping_thread &t) {
        requests += t.requests;
        checks += t.checks;
        compared += t.compared;
        success += t.success;
        cancel_schedule += t.cancel_schedule;
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
        std::cout << "Stops Requested:       " << totals.requests << "\n";
        std::cout << "Requests Checked:      " << totals.checks << "\n";
        std::cout << "Requests Compared:     " << totals.compared << "\n";
        std::cout << "Successfully Stopped:  " << totals.success << "\n";
        std::cout << "Schedule Cancelled:    " << totals.cancel_schedule << std::endl;
    }
};

#endif