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
        spin_lock lock{};
        std::deque<std::deque<path_node>> pool{};
        int depth{0};

    public:
        /*Grabs a node from the shallowest / zero pool*/
        bool pop_from_zero_list(path_node &result_node);
        
        /*Grabs a node from the deepest / active pool*/
        bool pop_from_active_list(path_node &result_node);

        /*Pushes new list to the back of the local pool*/
        void push_list(const std::deque<path_node> &list);
        
        /*Removes active list once empty*/
        void pop_active_list();
        
        /* Determines if a specific thread's local pool is completely empty. */
        bool out_of_work() { return pool.size() == 0; }
        
        //sets the relative depth of the pool
        void set_pool_depth(int depth) { this->depth = depth; }

        // diagnostic
        int pool_size() { return pool.size(); }

        // diagnostic
        int active_pool_size() { return pool.back().size(); }

        // value is compared when choosing which thread to steal from
        unsigned long long node_value();
    };
    
    class local_pool {
        private:
            int thread_count;
            std::vector<local_pool_thread> threads;

        public:
            local_pool(int thread_count);

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