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
    int time_limit = 3600;

    /* Size of initial global workload pool */
    int global_pool_size = 32;

    /* Percent of available system memory to use at maximum */
    float memory_percent = 0.9;

    /* Enable stealing work from other threads */
    bool enable_work_stealing = true;

    /* Enable stopping an inferior thread when a better cost is found */
    bool enable_thread_stopping = true;

    /* Enable using one thread to run LKH */
    bool enable_lkh = true;

    /* Number of buckets in history table */
    int number_of_buckets = 1;

    /* Size of each bucket in history table (set to 0 to split the entries in each bucket equally) */
    int bucket_size = 0;

    /* Time to run lkh for
     * If reuse thread is enabled, lkh thread will be reused for branch and bound after this time
     * Set to -1 for no time limit
     */
    int end_lkh_time = 20;

    /* Enable processing the best lkh tour into the history table */
    bool process_lkh_best_tour = true;

    /* Enable reusing the lkh thread for branch and bound after end_lkh_time */
    bool reuse_lkh_thread = true;

    /* Finish lkh before starting branch and bound, instead of running in parallel (for debugging) */
    bool finish_lkh_before_bb = false;

    /* Enable Heuristic: Treat 3 subtable as 1 table if the global pool is empty before hitting the first threshold */
    bool enable_heuristic = false;

    /* Process each subpath of the lkh tour as a separate history table entry */
    bool process_lkh_subpaths = true;

    /* If trace is enabled, set the detail level of the trace */
    int trace_detail_level = 0;

    /* UNUSED */

    int assign_workload_level = 150;
    int history_depth = -1;
    int exploitation_percent = 50;
    int group_sample_time = 3600;
    int group_thread_count = 4;
    bool enable_progress_estimation = true;
    int stable_lkh_entry_duration = 10;
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