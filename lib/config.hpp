#ifndef CONFIG_H
#define CONFIG_H

#include <fstream>
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

/**
 * @brief Stores parameters from config file.
 */
struct Config {
    /* Stops running the instance if it exceeds this time limit (seconds) */
    int time_limit;

    /* Size of initial global workload pool */
    int global_pool_size;

    /* Percent of available system memory to use at maximum */
    float memory_percent;

    /* Enable stealing work from other threads */
    bool enable_work_stealing;

    /* Enable stopping an inferior thread when a better cost is found */
    bool enable_thread_stopping;

    /* Enable using one thread to run LKH */
    bool enable_lkh;

    /* Number of buckets in history table */
    int number_of_buckets;

    /* Size of each bucket in history table (set to 0 to split the entries in each bucket equally) */
    int bucket_size;


    /* UNUSED */

    int assign_workload_level;
    int history_depth;
    float exploitation_percent;
    int group_sample_time;
    int group_thread_count;
    bool enable_progress_estimation;
    bool enable_heuristic;
};

/**
 * @brief Parses config file into struct.
 * 
 * Parses a config file from the specified path and returns a config_t struct
 * containing values set in the file.
 * 
 * @return A Config struct containing the parsed values.
 */
Config parse_config(string path);

#endif