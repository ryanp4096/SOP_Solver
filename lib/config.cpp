#include "config.hpp"


/** 
 * @brief Parses key and value strings out of a single line
 */
void parse_key_value(string line, string *key, string *value);

/**
 * @brief Trims leading and trailing whitespace of string.
 */
void trim(string& str);

/**
 * @brief Sets a specified key in the config struct to a given value.
 * 
 * Modify this function when adding new config options.
 */
void set_key(Config *config, string key, string value);


Config parse_config(string path) {
    ifstream infile(path);
    Config config;
    string line;
    int enable_count = 0;

    while (getline(infile, line))
    {
        // Parse key and value
        string key;
        string value;
        parse_key_value(line, &key, &value);
        if (key.empty() || value.empty()) continue;

        // Legacy: supports multiple ordered parameters named "enable"
        if (key == "enable") {
            key += to_string(enable_count);
            enable_count++;
        }
        set_key(&config, key, value);
        line.clear();
    }

    infile.close();
    return config;
}

void parse_key_value(string line, string *key, string *value) {
    // Remove comment
    const size_t comment_pos = line.find("//");
    const string content = line.substr(0, comment_pos);

    // Extract key and value
    const size_t equals_pos = content.find('=');
    if (equals_pos == string::npos) return;
    string k = content.substr(0, equals_pos);
    string v = content.substr(equals_pos + 1);

    // Trim whitespace
    trim(k);
    trim(v);

    // Convert key to lowercase
    std::transform(k.begin(), k.end(), k.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    *key = k;
    *value = v;
}

void trim(string& str) {
    const size_t start = str.find_first_not_of(" \t\n\r");
    const size_t end = str.find_last_not_of(" \t\n\r");

    if (start == string::npos)
        str.clear();
    else
        str = str.substr(start, end - start + 1);
}

void set_key(Config *config, string key, string value) {
    const char *cvalue = value.c_str();
    if (key == "time_limit")
        config->time_limit = atoi(cvalue);

    else if (key == "global_pool_size")
        config->global_pool_size = atoi(cvalue);

    else if (key == "level")
        config->assign_workload_level = atoi(cvalue);

    else if (key == "restrict_per")
        config->memory_percent = atof(cvalue);

    else if (key == "history_depth")
        config->history_depth = atoi(cvalue);

    else if (key == "ratio")
        config->exploitation_percent = atoi(cvalue);

    else if (key == "cycle_time")
        config->group_sample_time = atoi(cvalue);

    else if (key == "group_thread_count")
        config->group_thread_count = atoi(cvalue);

    else if (key == "enable_work_stealing" || key == "enable0")
        config->enable_work_stealing = atoi(cvalue);

    else if (key == "enable_thread_stopping" || key == "enable1")
        config->enable_thread_stopping = atoi(cvalue);

    else if (key == "enable_lkh" || key == "enable2")
        config->enable_lkh = atoi(cvalue);

    else if (key == "enable_progress_estimation" || key == "enable3")
        config->enable_progress_estimation = atoi(cvalue);

    else if (key == "number_of_buckets")
        config->number_of_buckets = atoi(cvalue);

    else if (key == "bucket_size")
        config->bucket_size = atoi(cvalue);

    else if (key == "enable_heuristic" || key == "enable4")
        config->enable_heuristic = atoi(cvalue);

    else if (key == "end_lkh_time")
        config->end_lkh_time = atoi(cvalue);
    
    else if (key == "process_lkh_best_tour")
        config->process_lkh_best_tour = atoi(cvalue);

    else if (key == "reuse_lkh_thread")
        config->reuse_lkh_thread = atoi(cvalue);
    
    else if (key == "finish_lkh_before_bb")
        config->finish_lkh_before_bb = atoi(cvalue);

    else if (key == "stable_lkh_entry_duration")
        config->stable_lkh_entry_duration = atoi(cvalue);

    else if (key == "process_lkh_subpaths")
        config->process_lkh_subpaths = atoi(cvalue);

    else if (key == "trace_detail_level")
        config->trace_detail_level = atoi(cvalue);

    else if (key == "enable_manual_match_check")
        config->enable_manual_match_check = atoi(cvalue);
    
    else {
        cout << "[Config] Unknown key " << key << endl;
        return;
    }
    
    cout << "[Config] " << key << ": " << value << endl;
}