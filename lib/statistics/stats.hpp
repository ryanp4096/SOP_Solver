#ifndef STATS_H
#define STATS_H

#include "../history_table/history_table.hpp"

#include "nodes.hpp"
#include "subpaths.hpp"
#include "prefixes.hpp"
#include "thread_stopping.hpp"
#include "work_stealing.hpp"
#include "cpu_timer.hpp"
#include "history_table.hpp"

/* Stats for a specific thread, pointing to multiple thread stats structs stored in stats_global*/
struct stats_thread {
    stats_nodes_thread nodes;
    stats_prefixes_thread prefixes;
    stats_subpaths_thread subpaths;
    stats_thread_stopping_thread thread_stopping;
    stats_work_stealing_thread work_stealing;
    stats_cpu_timer_thread cpu_timer;
    stats_history_table_thread history_table;

    stats_thread(unsigned instance_size)
        : nodes{instance_size},
          prefixes{},
          subpaths{instance_size},
          thread_stopping{},
          work_stealing{instance_size},
          cpu_timer{},
          history_table{}
          {}
};

/* Contains all stats for all threads */
struct stats_global {
    unsigned thread_count;
    unsigned instance_size;

    stats_nodes nodes;
    stats_prefixes prefixes;
    stats_subpaths subpaths;
    stats_thread_stopping thread_stopping;
    stats_work_stealing work_stealing;
    stats_cpu_timer cpu_timer;
    stats_history_table history_table;

    stats_global(unsigned thread_count, unsigned instance_size, const HistoryTable &history_table)
        : thread_count{thread_count},
          instance_size{instance_size},
          nodes{thread_count, instance_size},
          prefixes{thread_count},
          subpaths{thread_count, instance_size},
          thread_stopping{thread_count},
          work_stealing{thread_count},
          cpu_timer{thread_count},
          history_table{thread_count, history_table}
          {}
    
    void register_thread(unsigned thread_id, const stats_thread &thread) {
        nodes.register_thread(thread_id, thread.nodes);
        prefixes.register_thread(thread_id, thread.prefixes);
        subpaths.register_thread(thread_id, thread.subpaths);
        thread_stopping.register_thread(thread_id, thread.thread_stopping);
        work_stealing.register_thread(thread_id, thread.work_stealing);
        cpu_timer.register_thread(thread_id, thread.cpu_timer);
        history_table.register_thread(thread_id, thread.history_table);
    }
};

#endif