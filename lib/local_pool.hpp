#ifndef LOCAL_H
    #define LOCAL_H

    #include <deque> //for the local pool structure
    #include <queue>
    #include <iostream>
    #include "synchronization.hpp"
    #include "graph.hpp"

    /* The local pool construct, consisting of pools of nodes for each thread. In 
        each thread, the pool is organized by depth, so that nodes are stolen from 
        only the shallowest part of the pool, and added at the deepest level. */
    
    class local_pool_thread {
    private:
        int instance_size;
        spin_lock lock{};
        std::vector<std::deque<path_node>> pool;
        int zero_depth{0};
        int depth{0};

    public:
        local_pool_thread(int instance_size)
            : instance_size{instance_size}, pool(instance_size) {}

        local_pool_thread(const local_pool_thread &l)
            : instance_size{l.instance_size}, lock{}, pool{l.pool}, zero_depth{l.zero_depth}, depth{l.depth} {}

        local_pool_thread(local_pool_thread &&l)
            : instance_size{l.instance_size}, lock{}, pool{std::move(l.pool)}, zero_depth{l.zero_depth}, depth{l.depth} {}

        int level() { return depth - zero_depth; }

        void initial_depth(int init_depth);

        /*Grabs a node from the shallowest / zero pool*/
        bool pop_from_zero_list(path_node &result_node);
        
        /*Grabs a node from the deepest / active pool*/
        bool pop_from_active_list(path_node &result_node);

        /*Pushes new list to the back of the local pool*/
        void push_list(const std::deque<path_node> &list);
        
        /*Removes active list once empty*/
        void pop_active_list();

        // value is compared when choosing which thread to steal from
        unsigned long long node_value();
    };
    
    class local_pool {
        private:
            int thread_count;
            int instance_size;
            std::vector<local_pool_thread> threads;

        public:
            local_pool(int thread_count, int instance_size)
                : thread_count{thread_count}, instance_size{instance_size}, threads{}
                {
                    threads.reserve(thread_count);
                    for (int i = 0; i < thread_count; i++)
                        threads.push_back(local_pool_thread(instance_size));
                }

            /* Returns a specific thread's local pool */
            local_pool_thread &thread(int thread_number) { return threads[thread_number]; }

            /* Returns a thread number of the best victim, other than you, for workstealing. 
                thread_number - this thread's number, to ensure you aren't recommended to steal from yourself
                Return - the thread number of the thread to steal from */
            int choose_victim(int thread_number,std::vector<std::atomic<unsigned long long>>& work_remaining, int stolen_from);

            //diagnostic
            void print();
    };

#endif