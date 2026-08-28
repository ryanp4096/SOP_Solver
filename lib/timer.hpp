#ifndef TIMER_H
#define TIMER_H
#include <unordered_map>
#include <chrono>
#include <vector>
#include <array>
#include <cstdint>
#include <x86intrin.h>
#include <iostream>

class timer
{
private:
    std::unordered_map<std::string, std::chrono::time_point<std::chrono::system_clock>> markers;
    std::chrono::time_point<std::chrono::system_clock> epoch;

public:
    timer();                                     // contructor - starts timer
    void add_new_marker(std::string marker);     // adds a marker to the timer
    void restart();                              // restarts the main timer
    double get_time_seconds();                   // returns the main time in seconds
    double get_time_seconds(std::string marker); // returns the specified marker time in seconds
    double get_time_millis();                    // returns the main time in milliseconds
    double get_time_millis(std::string marker);  // returns the specified marker time in milliseconds
};

class cpu_timer
{
public:
    enum marker {
        NODE_SETUP, HISTORY_UTILIZATION, SUBPATH_HISTORY, NODE_END,
        POOL_SORT,
        RECURSIVE_THREAD_STOP, RECURSIVE_START, RECURSIVE_END,
        HISTORY_LOOKUP, LOWER_BOUND, ENTRY_EDIT, THREAD_STOP
    };
private:
    static constexpr int marker_count = 12;

    struct marker_data {
        uint64_t time{0};
        uint64_t count{0};
        uint64_t start{0};
    };

    int thread_count;
    std::vector<std::array<marker_data, marker_count>> data;

public:
    void initialize(int thread_count)
    {
        this->thread_count = thread_count;
        data.assign(thread_count, std::array<marker_data, marker_count>{});
    }

    inline void start(marker m, int thread_id)
    {
        uint64_t *d = &data[thread_id][m].start;
        _mm_lfence();
        *d = __rdtsc();
        _mm_lfence();
    }

    inline void stop(marker m, int thread_id)
    {
        unsigned temp;
        uint64_t end = __rdtscp(&temp);
        _mm_lfence();
        marker_data &d = data[thread_id][m];
        d.count++;
        d.time += end - d.start;
    }

    void print_results()
    {
        std::cout << "------ CPU TIMER RESULTS ------" << std::endl;
        for (int marker = 0; marker < marker_count; marker++) {
            uint64_t time = 0;
            uint64_t count = 0;
            for (int thread = 0; thread < thread_count; thread++) {
                time += data[thread][marker].time;
                count += data[thread][marker].count;
            }
            uint64_t average = count == 0 ? 0 : time / count;
            
            std::cout << "Marker: " << marker << "   ";
            std::cout << "Average: " << average << "   ";
            std::cout << "Count: " << count << "   ";
            std::cout << "Total: " << time << std::endl;
        }
        std::cout << "-------------------------------" << std::endl;
    }
};

#endif