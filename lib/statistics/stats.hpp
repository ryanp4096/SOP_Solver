#ifndef STATS_H
#define STATS_H

#include "nodes.hpp"
#include "subpaths.hpp"
#include "prefixes.hpp"
#include "thread_stopping.hpp"

/* Stats for a specific thread, pointing to multiple thread stats structs stored in stats_global*/
struct stats_thread {
    stats_nodes_thread *nodes{};
    stats_prefixes_thread *prefixes{};
    stats_subpaths_thread *subpaths{};
    stats_thread_stopping_thread *thread_stopping{};
};

/* Contains all stats for all threads */
struct stats_global {
    stats_nodes nodes;
    stats_prefixes prefixes;
    stats_subpaths subpaths;
    stats_thread_stopping thread_stopping;

    stats_global(unsigned thread_count, unsigned instance_size)
        : nodes{thread_count, instance_size},
          prefixes{thread_count},
          subpaths{thread_count, instance_size},
          thread_stopping{thread_count}
          {}
    
    stats_thread thread(unsigned thread_id) {
        return {
            .nodes = &nodes.threads[thread_id],
            .prefixes = &prefixes.threads[thread_id],
            .subpaths = &subpaths.threads[thread_id],
            .thread_stopping = &thread_stopping.threads[thread_id]
        };
    }
};

#endif