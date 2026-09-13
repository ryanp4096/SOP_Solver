#include "local_pool.hpp"
// // static bool isPrinted = false;
// // static bool isPrintedLast = false;
// // local_pool::local_pool(int thread_count)
// // {
// //     locks = std::vector<spin_lock>(thread_count);
// //     pools = std::vector<std::deque<std::deque<path_node>>>(thread_count);
// //     depths = std::vector<int>(thread_count);
// //     this->thread_count = thread_count;
// // }

// // void local_pool::print_top_sequence_sizes_end(int thread_total)
// // {
// //     std::vector<int> results;
// //     results.reserve(thread_total);

// //     if (!isPrintedLast)
// //     {
// //         isPrintedLast = true;
// //         for (int i = 0; i < thread_total; i++)
// //         {
// //             locks[i].lock(); // Lock the pool to ensure thread safety

// //             if (pools[i].size() > 0)
// //             {
// //                 while (pools[i].front().empty() && pools[i].size() > 1)
// //                 {
// //                     pools[i].pop_front();
// //                 }
// //                 if (pools[i].size() <= 1 && pools[i].front().empty())
// //                 {
// //                     results.push_back(-2); // Top deque is empty
// //                 }
// //                 else
// //                 {
// //                     // Access the top deque and get the first element in it
// //                     const path_node &top_node = pools[i].front().front();
// //                     // Add the size of the sequence to results
// //                     results.push_back(top_node.sequence.size());
// //                 }
// //             }
// //             else
// //             {
// //                 results.push_back(-1); // Pool is empty
// //             }

// //             locks[i].unlock(); // Unlock the pool after accessing
// //         }

// //         // Print all collected results as a comma-separated string
// //         for (size_t i = 0; i < results.size(); ++i)
// //         {
// //             std::cout << results[i];
// //             if (i < results.size() - 1)
// //                 std::cout << ",";
// //         }
// //         std::cout << std::endl;
// //     }
// // }

// // void local_pool::print_top_sequence_sizes(int thread_total)
// // {
// //     std::vector<int> results;
// //     results.reserve(thread_total);

// //     if (!isPrinted)
// //     {
// //         isPrinted = true;
// //         for (int i = 0; i < thread_total; i++)
// //         {
// //             locks[i].lock(); // Lock the pool to ensure thread safety

// //             if (pools[i].size() > 0)
// //             {
// //                 while (pools[i].front().empty() && pools[i].size() > 1)
// //                 {
// //                     pools[i].pop_front();
// //                 }
// //                 if (pools[i].size() <= 1 && pools[i].front().empty())
// //                 {
// //                     results.push_back(-2); // Top deque is empty
// //                 }
// //                 else
// //                 {
// //                     // Access the top deque and get the first element in it
// //                     const path_node &top_node = pools[i].front().front();
// //                     // Add the size of the sequence to results
// //                     results.push_back(top_node.sequence.size());
// //                 }
// //             }
// //             else
// //             {
// //                 results.push_back(-1); // Pool is empty
// //             }

// //             locks[i].unlock(); // Unlock the pool after accessing
// //         }

// //         // Print all collected results as a comma-separated string
// //         for (size_t i = 0; i < results.size(); ++i)
// //         {
// //             std::cout << results[i];
// //             if (i < results.size() - 1)
// //                 std::cout << ",";
// //         }
// //         std::cout << std::endl;
// //     }
// // }

// void local_pool::initial_state(int thread_number, const std::vector<int> &path, const boost::dynamic_bitset<> &key)
// {
//     assert(path.size() >= 1);
//     threads[thread_number].lock.lock();
//     threads[thread_number].current_path = path;
//     threads[thread_number].current_node = path.back();
//     threads[thread_number].current_path.pop_back();
//     threads[thread_number].current_key = key;
//     threads[thread_number].current_key[threads[thread_number].current_node] = false;
//     threads[thread_number].zero_depth = path.size();
//     threads[thread_number].depth = path.size();
//     threads[thread_number].lock.unlock();
// }

// bool local_pool::pop_from_zero_list(int thread_number, path_node &result_node, int stealing_thread)
// {
//     local_pool_thread &thread = threads[thread_number];
//     if (level(thread_number) <= 1 || thread.zero_list().queue.empty())
//     {
//         return false;
//     }
//     thread.lock.lock();

//     if (level(thread_number) <= 1 || thread.zero_list().queue.empty())
//     {
//         thread.lock.unlock();
//         return false;
//     }

//     int taken_node = thread.zero_list().queue.back();
//     local_pool_node &node = thread.zero_list().nodes[taken_node];

//     result_node.current_node_value = node.current_node_value;
//     result_node.lower_bound = node.lower_bound;
//     result_node.origin_node = thread.zero_depth <= 0 ? -1 : (thread.zero_depth == 1 ? taken_node : thread.current_path.at(1));
//     result_node.history_key.bit_vector.reset();
//     result_node.history_key.bit_vector.resize(instance_size, false);
//     result_node.sequence.clear();
//     result_node.sequence.reserve(thread.zero_depth + 1);
//     for (int d = 0; d < thread.zero_depth; d++) {
//         int tn = thread.current_path[d];
//         result_node.sequence.push_back(tn);
//         result_node.history_key.bit_vector[tn] = true;
//     }
//     result_node.sequence.push_back(taken_node);
//     result_node.history_key.bit_vector[taken_node] = true;
//     result_node.history_key.last_node = taken_node;

//     thread.zero_list().queue.pop_back();

//     thread.lock.unlock();
//     return true;
// };

// bool local_pool::pop_from_active_list(int thread_number, path_node &result_node)
// {

//     if (level(thread_number) <= 0)
//         return false;

//     deque<int> &queue = threads[thread_number].active_list().queue;
//     if (queue.empty())
//     {
//         return false;
//     }

//     int taken_node = queue.back();
//     local_pool_thread &thread = threads[thread_number];
//     thread.current_node = taken_node;
//     local_pool_node &node = thread.zero_list().nodes[taken_node];

//     result_node.current_node_value = node.current_node_value;
//     result_node.lower_bound = node.lower_bound;
//     result_node.origin_node = thread.current_path.size() == 0 ? -1 : (thread.current_path.size() == 1 ? thread.current_node : thread.current_path.at(1));
//     result_node.sequence.reserve(thread.current_path.size() + 1);
//     result_node.sequence = thread.current_path;
//     result_node.sequence.push_back(thread.current_node);
//     result_node.history_key.bit_vector = thread.current_key;
//     result_node.history_key.bit_vector[thread.current_node] = true;
//     result_node.history_key.last_node = thread.current_node;

//     queue.pop_back();

//     threads[thread_number].current_path.back() = taken_node;

//     return true;
// };

// // void local_pool::push_list(int thread_number, std::deque<path_node> list)
// // {
// //     locks[thread_number].lock();

// //     pools[thread_number].push_back(list);

// //     locks[thread_number].unlock();
// // };

// void local_pool::pop_active_list(int thread_number)
// {
//     threads[thread_number].lock.lock();
//     if (level(thread_number) > 0) {
//         assert(threads[thread_number].active_list().queue.empty());
//         threads[thread_number].pop();
//     }
//     threads[thread_number].lock.unlock();
// };

// void local_pool::push_to_ready_list(int thread_number, path_node &path_node)
// {
//     int taken_node = path_node.sequence.back();
//     threads[thread_number].ready_list().queue.push_back(taken_node);
//     local_pool_node &node = threads[thread_number].ready_list().nodes[taken_node];
//     node.current_node_value = path_node.current_node_value;
//     node.lower_bound = path_node.lower_bound;
// }

// void local_pool::push_ready_list(int thread_number, unsigned long long next_work_above)
// {
//     threads[thread_number].lock.lock();
//     deque<int> &queue = threads[thread_number].ready_list().queue;
//     std::vector<local_pool_node> &nodes = threads[thread_number].ready_list().nodes;
//     for (auto it = queue.begin(); it != queue.end(); ++it) {
//         nodes[*it].current_node_value = next_work_above;
//     }
//     if (!queue.empty())
//         std::sort(queue.begin(), queue.end(), [nodes](int src, int dst){ return nodes[src].lower_bound > nodes[dst].lower_bound; });
//     threads[thread_number].push();
//     threads[thread_number].lock.unlock();
// }


// bool local_pool::out_of_work(int thread_number)
// {
//     return level(thread_number) <= 0;
// };

int local_pool::choose_victim(int thread_number, std::vector<std::atomic<unsigned long long>> &work_remaining, int stolen_from)
{
    unsigned long long max_value = 0;
    int max_id = -1;
    bool flag = false;
    for (int i = 0; i < thread_count; i++)
    {
        if ((stolen_from & (1 << i)) != 0 || i == thread_number)
            continue;
        unsigned long long node_value = threads[i].steal_value();
        if (node_value > max_value)
        {
            max_value = node_value;
            max_id = i;
            continue;
        }
        if (max_value == 0 && (!flag || work_remaining[i] > work_remaining[max_id]))
        {
            max_value = node_value;
            max_id = i;
            flag = true;
        }
    }
    return max_id;
}

// // int local_pool::choose_victim(int thread_number, std::vector<std::atomic<unsigned long long>>& work_remaining, int stolen_from){
// //     double max_value = -1;
// //     int max_id = -1;
// //     //std::cout << stolen_from << std::endl;
// //     for(int i = 0; i < thread_count; i++){
// //         if((stolen_from & (1 << i)) != 0 )
// //             continue;
// //         if(i == thread_number)
// //             continue;
// //         if(work_remaining[i] / (depths[i] + 1) > max_value){
// //             max_value = work_remaining[i] / (depths[i] + 1);
// //             max_id = i;
// //             continue;
// //         }
// //         // if(work_remaining[i] == max_value && depths[i] < depths[max_id]){
// //         //     max_value = work_remaining[i];
// //         //     max_id = i;
// //         // }
// //     }
// //     return max_id;
// // }

// // int local_pool::choose_victim(int thread_number, std::vector<std::atomic<unsigned long long>>& work_remaining, int a){
// //     int target = rand() % 30;
// //     while(target == thread_number) target = rand() % 30;
// //     return target;
// // }

// int local_pool::active_pool_size(int thread_number)
// { // TODO: this is not strictly necessary
//     return threads[thread_number].active_list().queue.size();
// }

void local_pool::print()
{
    std::cout << "local_pool relative depth by thread: ";
    for (int i = 0; i < threads.size(); i++)
    {
        std::cout << threads[i].level() << ", ";
    }
    std::cout << std::endl;
}

// void local_pool::set_pool_depth(int thread_id, int depth)
// {
//     // depths[thread_id] = depth;
// }