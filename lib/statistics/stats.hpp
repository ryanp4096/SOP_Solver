#ifndef STATS_H
#define STATS_H

#include "nodes.hpp"
#include "subpaths.hpp"
#include "prefixes.hpp"
#include "thread_stopping.hpp"
#include "work_stealing.hpp"
#include "cpu_timer.hpp"

/* Stats for a specific thread, pointing to multiple thread stats structs stored in stats_global*/
struct stats_thread {
    stats_nodes_thread *const nodes{};
    stats_prefixes_thread *const prefixes{};
    stats_subpaths_thread *const subpaths{};
    stats_thread_stopping_thread *const thread_stopping{};
    stats_work_stealing_thread *const work_stealing{};
    stats_cpu_timer_thread *const cpu_timer{};
};

/* Contains all stats for all threads */
struct stats_global {
    stats_nodes nodes;
    stats_prefixes prefixes;
    stats_subpaths subpaths;
    stats_thread_stopping thread_stopping;
    stats_work_stealing work_stealing;
    stats_cpu_timer cpu_timer;

    stats_global(unsigned thread_count, unsigned instance_size)
        : nodes{thread_count, instance_size},
          prefixes{thread_count},
          subpaths{thread_count, instance_size},
          thread_stopping{thread_count},
          work_stealing{thread_count},
          cpu_timer{thread_count}
          {}
    
    stats_thread thread(unsigned thread_id) {
        return {
            .nodes = &nodes.threads[thread_id],
            .prefixes = &prefixes.threads[thread_id],
            .subpaths = &subpaths.threads[thread_id],
            .thread_stopping = &thread_stopping.threads[thread_id],
            .work_stealing = &work_stealing.threads[thread_id],
            .cpu_timer = &cpu_timer.threads[thread_id]
        };
    }
};

#endif