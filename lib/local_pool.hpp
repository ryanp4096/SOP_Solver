#ifndef LOCAL_H
    #define LOCAL_H

    #include <deque> //for the local pool structure
    #include <queue>
    #include <iostream>
    #include "synchronization.hpp"
    #include "graph.hpp"

    struct local_pool_node {
        NodeState state{UNEXPLORED};
        int lower_bound{-1};
        unsigned long long current_node_value = -1;
        HistoryNode *prefix_node{};
    };

    // struct local_pool_node_ref {
    //     local_pool_thread *thread;
    //     int depth;
    //     int taken_node;
    //     local_pool_node *content;

    //     local_pool_node_ref()
    //         : thread{}, depth{-1}, taken_node{-1}, content{} {}

    //     local_pool_node_ref(local_pool_thread *thread, int depth, int taken_node, local_pool_node *content)
    //         : thread(thread), depth(depth), taken_node(taken_node), content(content) {}

    //     bool valid() { return content != nullptr; }
    // };

    class local_pool_list {
    private:
        std::vector<local_pool_node> nodes;
        boost::dynamic_bitset<> available;
        std::deque<int> queue{};
        int active_node{-1};

    public:
        local_pool_list(unsigned instance_size)
            : nodes(instance_size), available(instance_size, false) {}
        
        local_pool_list(const local_pool_list &l)
            : nodes(l.nodes), available(l.available), queue(l.queue), active_node(l.active_node) {}

        void clear() {
            available.reset();
            queue.clear();
            active_node = -1;
        }

        bool empty() { return queue.empty(); }

        void push(int node, int lower_bound = -1, unsigned long long current_node_value = -1) {
            nodes[node] = {
                .lower_bound = lower_bound,
                .current_node_value = current_node_value
            };
            available[node] = true;
            queue.push_back(node);
        }

        void sort() {
            if (queue.empty()) return;
            std::sort(
                queue.begin(), queue.end(),
                [&](int src, int dst){ return nodes[src].lower_bound > nodes[dst].lower_bound; }
            );
        }

        int pop() {
            int node = queue.back();
            available[node] = false;
            queue.pop_back();
            return node;
        }

        int back() {
            return queue.back();
        }

        local_pool_node *get_node(int node) {
            // if (!available[node]) return nullptr;
            return &nodes[node];
        }

        void activate(int node) {
            active_node = node;
        }

        void deactivate() {
            active_node = -1;
        }

        void set_work(unsigned long long work) {
            for (int n : queue) {
                nodes[n].current_node_value = work;
            }
        }
    };

    class local_pool_thread {
    private:
        int instance_size;
        std::vector<local_pool_list> lists;
        spin_lock lock{};
        int base_depth{-1};
        int depth{-1};
    public:
        std::vector<int> current_path{};
        boost::dynamic_bitset<> current_key;

    public:
        local_pool_thread(unsigned instance_size)
            : instance_size(instance_size), lists{}, current_key(instance_size, false)
        {
            lists.reserve(instance_size);
            for (int i = 0; i < instance_size; i++) lists.push_back(local_pool_list(instance_size));
            current_path.reserve(instance_size);
        }

        local_pool_thread(const local_pool_thread &l)
            : instance_size(l.instance_size), lists(l.lists), base_depth{l.base_depth}, depth{l.depth}, current_path{l.current_path}, current_key(l.current_key) {}

        local_pool_list &active_list() { return lists[depth]; }
        local_pool_list &base_list() { return lists[base_depth]; }

        int level() { return depth - base_depth; }

        void initial_state(const std::vector<int> &path, const boost::dynamic_bitset<> &key) {
            assert(path.size() >= 1);

            lock.lock();
            current_path = path;
            current_key = key;
            base_depth = path.size();
            depth = path.size();
            lock.unlock();
        }

        void activate(int node) {
            lock.lock();
            current_path.push_back(node);
            current_key[node] = true;
            active_list().activate(node);
            depth++;
            lock.unlock();
        }

        void deactivate() {
            lock.lock();
            if (depth < instance_size) active_list().clear();
            depth--;
            active_list().deactivate();
            int node = current_path.back();
            current_path.pop_back();
            current_key[node] = false;
            lock.unlock();
        }

        void push(int node, int lower_bound = -1, unsigned long long current_node_value = -1) {
            if (depth < instance_size) active_list().push(node, lower_bound, current_node_value);
        }

        void sort() {
            if (depth < instance_size) active_list().sort();
        }
        void set_work(unsigned long long work) {
            if (depth < instance_size) active_list().set_work(work);
        }

        bool pop(path_node &result_node) {
            if (active_list().empty()) return false;
            int taken_node = active_list().pop();
            local_pool_node *content = active_list().get_node(taken_node);

            result_node.current_node_value = content->current_node_value;
            result_node.lower_bound = content->lower_bound;
            result_node.sequence = current_path;
            result_node.sequence.push_back(taken_node);
            result_node.history_key.last_node = taken_node;
            result_node.history_key.bit_vector = current_key;
            result_node.history_key.bit_vector[taken_node] = true;
            result_node.origin_node = result_node.sequence[1];

            return true;
        }

        bool steal(path_node &result_node) {
            if (depth - base_depth < 1 || base_list().empty()) return false;
            lock.lock();
            if (depth - base_depth < 1 || base_list().empty()) {
                lock.unlock();
                return false;
            }
            int taken_node = base_list().pop();
            local_pool_node *content = base_list().get_node(taken_node);

            std::vector<int> sequence{};
            boost::dynamic_bitset<> bit_vector(instance_size, false);
            for (int i = 0; i < base_depth; i++) {
                sequence.push_back(current_path[i]);
                bit_vector[current_path[i]] = true;
            }
            sequence.push_back(taken_node);
            bit_vector[taken_node] = true;

            result_node.current_node_value = content->current_node_value;
            result_node.lower_bound = content->lower_bound;
            result_node.sequence = std::move(sequence);
            result_node.history_key.last_node = taken_node;
            result_node.history_key.bit_vector = std::move(bit_vector);
            result_node.origin_node = result_node.sequence[1];

            lock.unlock();
            return true;
        }

        unsigned long long steal_value() {
            if (depth - base_depth < 1 || base_list().empty()) return 0;
            lock.lock();
            if (depth - base_depth < 1 || base_list().empty()) {
                lock.unlock();
                return 0;
            }
            unsigned long long value = base_list().get_node(base_list().back())->current_node_value;
            lock.unlock();
            return value;
        }
    };

    /* The local pool construct, consisting of pools of nodes for each thread. In 
        each thread, the pool is organized by depth, so that nodes are stolen from 
        only the shallowest part of the pool, and added at the deepest level. */
    class local_pool {
        // private:
        //     struct local_pool_node {
        //         NodeState state{UNEXPLORED};
        //         int lower_bound{-1};
        //         unsigned long long current_node_value = -1;
        //         HistoryNode *prefix_node{};
        //     };

        //     struct local_pool_list {
        //         std::vector<local_pool_node> nodes;
        //         std::deque<int> queue{};
        //         int current_node{-1};

        //         local_pool_list(unsigned instance_size) : nodes(instance_size) {}
        //         local_pool_list(const local_pool_list &l) : nodes(l.nodes) {}
        //     };

        //     struct local_pool_thread {
        //         std::vector<local_pool_list> lists;
        //         spin_lock lock{};
        //         int zero_depth{-1};
        //         int depth{-1};
        //         std::vector<int> current_path{};
        //         boost::dynamic_bitset<> current_key;
        //         int current_node{-1};

        //         local_pool_thread(unsigned instance_size) : current_key(instance_size, false) {
        //             current_path.reserve(instance_size);
        //             lists = std::vector<local_pool_list>(instance_size, local_pool_list(instance_size));
        //         }
        //         local_pool_thread(const local_pool_thread &l) : current_key{l.current_key}, lists{l.lists} { current_path.reserve(l.current_key.size()); }
        //         local_pool_list &get_list(int d) { return lists[d - 1]; }
        //         local_pool_list &ready_list() { return lists[depth]; }
        //         local_pool_list &active_list() { return lists[depth - 1]; }
        //         local_pool_list &zero_list() { return lists[zero_depth]; }
        //         int level() { return depth - zero_depth; }
        //         void push() {
        //             assert(current_node != -1);
        //             current_path.push_back(current_node);
        //             current_key[current_node] = true;
        //             current_node = -1;
        //             depth++;
        //         }
        //         void pop() {
        //             int node = current_path.back();
        //             current_path.pop_back();
        //             current_key[node] = false;
        //             current_node = node;
        //             depth--;
        //         }
        //     };

            int thread_count;
            unsigned instance_size;
            std::vector<local_pool_thread> threads;

            // int current_target = 0;
            // std::vector<spin_lock> locks;
            // std::vector<std::deque<std::deque<path_node>>> pools;
            // std::vector<int> depths;
            // int level(int thread_count) {
            //     if (threads[thread_count].zero_depth == -1) return -1;
            //     return threads[thread_count].depth - threads[thread_count].zero_depth;
            // }
        public:
            local_pool(int thread_count, unsigned instance_size)
                : thread_count{thread_count}, instance_size{instance_size}, threads{}
            {
                threads.reserve(thread_count);
                for (int i = 0; i < thread_count; i++) threads.push_back(local_pool_thread(instance_size));
            }
            
            local_pool_thread &thread(int i) { return threads[i]; }

            // void print_top_sequence_sizes(int thread_total);
            // void print_top_sequence_sizes_end(int thread_total);
            // void add_to_depth_queue(int thread);
            // void initial_state(int thread_number, const std::vector<int> &path, const boost::dynamic_bitset<> &key);
            // /*Grabs a node from the shallowest / zero pool*/
            // bool pop_from_zero_list(int thread_number, path_node &result_node, int stealing_thread);
            // /*Grabs a node from the deepest / active pool*/
            // bool pop_from_active_list(int thread_number, path_node &result_node);
            // /*Pushes new list to the back of the local pool*/
            // // void push_list(int thread_number, std::deque<path_node> list);
            // /*removes */
            // void pop_active_list(int thread_number);

            // void push_to_ready_list(int thread_number, path_node &path_node);

            // void push_ready_list(int thread_number, unsigned long long next_work_above);

            // /* Determines if a specific thread's local pool is completely empty. */
            // bool out_of_work(int thread_number);
            // /* Returns a thread number of the best victim, other than you, for workstealing. 
            //     thread_number - this thread's number, to ensure you aren't recommended to steal from yourself
            //     Return - the thread number of the thread to steal from */
            int choose_victim(int thread_number,std::vector<std::atomic<unsigned long long>>& work_remaining, int stolen_from);
            // //diagnostic
            // int active_pool_size(int thread_number); 
            // //diagnostic
            void print();
            // //sets the relative depth of the pool
            // void set_pool_depth(int thread_id, int depth);
    };

#endif