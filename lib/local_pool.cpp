#include "local_pool.hpp"

local_pool::local_pool(int thread_count)
{
    threads = std::vector<local_pool_thread>(thread_count);
    this->thread_count = thread_count;
}

bool local_pool_thread::pop_from_zero_list(path_node &result_node)
{
    if (pool.size() <= 1)
    {
        return false;
    }
    lock.lock();

    while (pool.front().empty() && pool.size() > 1)
    {
        pool.pop_front();
        depth++;
    }

    if (pool.size() <= 1)
    {
        lock.unlock();
        return false;
    }

    result_node = pool.front().back();
    pool.front().pop_back();
    // depths[stealing_thread] = depth + 1;

    if (pool.front().empty())
    {
        pool.pop_front();
        depth++;
    }

    lock.unlock();
    return true;
};

bool local_pool_thread::pop_from_active_list(path_node &result_node)
{

    if (pool.size() == 0)
        return false;

    if (pool.size() == 0 || pool.back().empty())
    {
        if (pool.size() == 1)
        {
            lock.unlock();
        }
        return false;
    }

    result_node = pool.back().back();
    pool.back().pop_back();

    return true;
};

void local_pool_thread::push_list(const std::deque<path_node> &list)
{
    lock.lock();

    pool.push_back(list);

    lock.unlock();
};

void local_pool_thread::pop_active_list()
{
    lock.lock();
    if (pool.size() > 0)
        pool.pop_back();
    lock.unlock();
};

// bool local_pool_thread::out_of_work()
// {
//     return pool.size() == 0;
// };

unsigned long long local_pool_thread::node_value()
{
    lock.lock();
    unsigned long long node_value = 0;
    if (pool.size() > 1 && pool.front().size() != 0)
    {
        node_value = pool.front().back().current_node_value;
    }
    lock.unlock();
    return node_value;
}

int local_pool::choose_victim(int thread_number, std::vector<std::atomic<unsigned long long>> &work_remaining, int stolen_from)
{
    unsigned long long max_value = 0;
    int max_id = -1;
    bool flag = false;
    for (int i = 0; i < thread_count; i++)
    {
        if ((stolen_from & (1 << i)) != 0 || i == thread_number)
            continue;
        unsigned long long node_value = threads[i].node_value();
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

// int local_pool_thread::active_pool_size()
// {
//     return pool.back().size();
// }

void local_pool::print()
{
    for (int i = 0; i < threads.size(); i++)
    {
        std::cout << threads[i].pool_size() << ", ";
    }
    std::cout << std::endl;
}

// void local_pool_thread::set_pool_depth(int depth)
// {
//     this->depth = depth;
// }