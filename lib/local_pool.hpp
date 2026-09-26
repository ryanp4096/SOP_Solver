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

    class local_pool_thread;

    struct local_pool_node {
        int lower_bound;
        unsigned long long current_node_value;
        HistoryNode *history_node;
        int cost;
    };

    struct local_pool_node_ref {
        local_pool_node *node;
        local_pool_thread *thread;
        int taken_node;
        int depth;

        inline int lower_bound() const { return node->lower_bound; }
        inline unsigned long long current_node_value() const { return node->current_node_value; }
        inline HistoryNode *history_node() const { return node->history_node; }
        inline int cost() const { return node->cost; }
        std::vector<int> sequence() const;
        boost::dynamic_bitset<> bit_vector() const;
        void to_path_node(path_node &node) const;
    };

    class local_pool_list {
    private:
        int instance_size;
        std::vector<local_pool_node> nodes;
        std::vector<int> queue{};
        boost::dynamic_bitset<> present;
    
    public:
        local_pool_list(int instance_size)
            : instance_size{instance_size}, nodes(instance_size), present(instance_size, false)
            { queue.reserve(instance_size); }
        
        local_pool_list(const local_pool_list &l)
            : instance_size{l.instance_size}, nodes{l.nodes}, queue{l.queue}, present{l.present}
            { queue.reserve(instance_size); }

        local_pool_list(local_pool_list &&l)
            : instance_size{l.instance_size}, nodes{std::move(l.nodes)}, queue{std::move(l.queue)}, present{std::move(l.present)}
            { queue.reserve(instance_size); }
        
        /* Empty the list */
        void clear();

        /* Check if the queue is empty */
        bool empty() { return queue.empty(); }

        int back() { return queue.back(); }

        /* Get the next node in the queue */
        local_pool_node &back_node() { return nodes[queue.back()]; }

        /* Pop a node from the queue */
        void pop_back() { queue.pop_back(); }

        /* Push a node to the queue */
        void push_back(int taken_node, const local_pool_node &node);

        /* Push a node to the queue */
        void push_back(int taken_node, local_pool_node &&node);

        /* Update work remaining values after all nodes have been added */
        void set_node_value(unsigned long long next_work_above);

        /* Sort queue by lower bound after all nodes have been added */
        void sort();

        /* Find a specific path */
        local_pool_node *get(int last_node);
    };

    class local_pool_thread {
    private:
        int instance_size;
        spin_lock lock{};
        std::vector<local_pool_list> pool;
    public:
        std::vector<int> current_path{};
        boost::dynamic_bitset<> current_key;
    private:
        int zero_depth{0};
        int depth{0};

    public:
        local_pool_thread(int instance_size)
            : instance_size{instance_size}, pool{}, current_key(instance_size, false)
            {
                current_path.reserve(instance_size);
                pool.reserve(instance_size);
                for (int i = 0; i < instance_size; i++)
                    pool.push_back(local_pool_list(instance_size));
            }

        local_pool_thread(const local_pool_thread &l)
            : instance_size{l.instance_size}, lock{}, pool{l.pool}, current_path{l.current_path}, current_key{l.current_key}, zero_depth{l.zero_depth}, depth{l.depth} {}

        local_pool_thread(local_pool_thread &&l)
            : instance_size{l.instance_size}, lock{}, pool{std::move(l.pool)}, current_path{std::move(l.current_path)}, current_key{std::move(l.current_key)}, zero_depth{l.zero_depth}, depth{l.depth} {}

        /* The number of lists in the pool that contain */
        int level() { return depth - zero_depth; }

        /* Establishes the depth of the problem state before enumeration */
        void initial_state(const std::vector<int> &sequence);

        /*Grabs a node from the shallowest / zero pool*/
        bool pop_from_zero_list(path_node &result_node);
        
        /*Grabs a node from the deepest / active pool*/
        bool pop_from_active_list(local_pool_node_ref &result_node);

        /*Initializes new list that will later be pushed to the back of the local pool*/
        void start_ready_list();

        /* Pushes a new node to the ready list */
        void push_to_ready_list(int taken_node, int lower_bound, HistoryNode *history_node, int cost);

        /* Updates work remaining values of nodes in ready list */
        void update_ready_list_node_value(unsigned long long node_value) { pool[depth].set_node_value(node_value); }

        /* Sorts the ready list by lower bound */
        void sort_ready_list() { pool[depth].sort(); }

        /*Pushes new list to the back of the local pool*/
        void push_ready_list();
        
        /*Removes active list once empty*/
        void pop_active_list();

        // value is compared when choosing which thread to steal from
        unsigned long long node_value();

        /* Find a specific path*/
        bool get(int depth, int last_node, local_pool_node_ref &result_node);

        /* Add to the current path */
        void activate(int taken_node);

        /* Remove from the current path */
        void deactivate();
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