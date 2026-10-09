#ifndef STATS_WORK_STEALING_H
#define STATS_WORK_STEALING_H

#include <vector>
#include <iostream>

struct stats_work_stealing_target {
    unsigned long long attempts{0};     // tried to steal from this thread
    unsigned long long success{0};      // successfully stole from this thread

    stats_work_stealing_target &operator+=(const stats_work_stealing_target &t) {
        attempts += t.attempts;
        success += t.success;
        return *this;
    }
};

struct stats_work_stealing_thread {
    unsigned long long attempts{0};     // tried to steal from another thread
    unsigned long long success{0};      // successfully stole work from another thread
    double time_stealing{0.0};          // time spent work stealing

    std::vector<stats_work_stealing_target> by_target;  // stats by thread targeted

    stats_work_stealing_thread(unsigned thread_count) : by_target(thread_count) {}

    stats_work_stealing_thread &operator+=(const stats_work_stealing_thread &t) {
        attempts += t.attempts;
        success += t.success;
        time_stealing += t.time_stealing;
        for (std::size_t i = 0; i < by_target.size(); i++)
            by_target[i] += t.by_target[i];
        return *this;
    }
};

struct stats_work_stealing {
    std::vector<stats_work_stealing_thread> threads{};
    unsigned thread_count;

    stats_work_stealing(unsigned thread_count)
        : thread_count{thread_count}
    {
        threads.reserve(thread_count);
        for (unsigned i = 0; i < thread_count; i++) {
            threads.push_back(stats_work_stealing_thread(thread_count));
        }
    }

    double total_time_stealing() {
        double total = 0.0;
        for (unsigned i = 0; i < thread_count; i++)
            total += threads[i].time_stealing;
        return total;
    }

    void print_results() {
        stats_work_stealing_thread totals{thread_count};
        for (unsigned i = 0; i < thread_count; i++)
            totals += threads[i];
        
        std::cout << "WORK STEALING\n";
        std::cout << "Attempts:         " << totals.attempts << "\n";
        std::cout << "Successes:        " << totals.success << "\n";
        std::cout << "By Thread:\n";
        std::cout << "    Stole:        ";
        for (unsigned i = 0; i < thread_count - 1; i++)
            std::cout << threads[i].success << ", ";
        std::cout << threads[thread_count - 1].success << "\n";
        std::cout << "    Stolen From:  ";
        for (unsigned i = 0; i < thread_count - 1; i++)
            std::cout << totals.by_target[i].success << ", ";
        std::cout << totals.by_target[thread_count - 1].success << std::endl;
    }
};

#endif