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
    
    class local_pool_list {
    private:
        int instance_size;
        std::vector<path_node> nodes;
        std::vector<int> queue{};
        int last_popped{-1};
        boost::dynamic_bitset<> present;
    
    public:
        local_pool_list(int instance_size)
            : instance_size{instance_size}, nodes(instance_size), present(instance_size, false)
            { queue.reserve(instance_size); }
        
        local_pool_list(const local_pool_list &l)
            : instance_size{l.instance_size}, nodes{l.nodes}, queue{l.queue}, last_popped{l.last_popped}, present{l.present}
            { queue.reserve(instance_size); }

        local_pool_list(local_pool_list &&l)
            : instance_size{l.instance_size}, nodes{std::move(l.nodes)}, queue{std::move(l.queue)}, last_popped{l.last_popped}, present{std::move(l.present)}
            { queue.reserve(instance_size); }
        
        /* Empty the list */
        void clear();

        /* Check if the queue is empty */
        bool empty() { return queue.empty(); }

        /* Get the next node in the queue */
        path_node &back() { return nodes[queue.back()]; }

        /* Pop a node from the queue */
        void pop_back();

        /* Push a node to the queue */
        void push_back(const path_node &node);

        /* Push a node to the queue */
        void push_back(path_node &&node);

        /* Update work remaining values after all nodes have been added */
        void set_node_value(unsigned long long next_work_above);

        /* Sort queue by lower bound after all nodes have been added */
        void sort();

        /* Find a specific path */
        path_node *get(int last_node);
    };

    class local_pool_thread {
    private:
        int instance_size;
        spin_lock lock{};
        std::vector<local_pool_list> pool;
        int zero_depth{0};
        int depth{0};

    public:
        local_pool_thread(int instance_size)
            : instance_size{instance_size}, pool{}
            {
                pool.reserve(instance_size);
                for (int i = 0; i < instance_size; i++)
                    pool.push_back(local_pool_list(instance_size));
            }

        local_pool_thread(const local_pool_thread &l)
            : instance_size{l.instance_size}, lock{}, pool{l.pool}, zero_depth{l.zero_depth}, depth{l.depth} {}

        local_pool_thread(local_pool_thread &&l)
            : instance_size{l.instance_size}, lock{}, pool{std::move(l.pool)}, zero_depth{l.zero_depth}, depth{l.depth} {}

        /* The number of lists in the pool that contain */
        int level() { return depth - zero_depth; }

        /* The next list to be added */
        local_pool_list &ready_list() { return pool[depth]; }

        /* Establishes the depth of the problem state before enumeration */
        void initial_depth(int init_depth);

        /*Grabs a node from the shallowest / zero pool*/
        bool pop_from_zero_list(path_node &result_node);
        
        /*Grabs a node from the deepest / active pool*/
        bool pop_from_active_list(path_node &result_node);

        /*Initializes new list that will later be pushed to the back of the local pool*/
        void start_ready_list();

        /*Pushes new list to the back of the local pool*/
        void push_ready_list();
        
        /*Removes active list once empty*/
        void pop_active_list();

        // value is compared when choosing which thread to steal from
        unsigned long long node_value();

        /* Find a specific path*/
        path_node *get(int depth, int last_node);
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