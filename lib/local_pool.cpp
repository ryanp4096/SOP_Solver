#include "local_pool.hpp"

void local_pool_list::clear() {
    queue.clear();
    last_popped = -1;
    present.reset();
}

void local_pool_list::pop_back() {
    last_popped = queue.back();
    queue.pop_back();
}

void local_pool_list::push_back(const path_node &node) {
    int last_node = node.sequence.back();
    present[last_node] = true;
    nodes[last_node] = node;
    queue.push_back(last_node);
}

void local_pool_list::push_back(path_node &&node) {
    int last_node = node.sequence.back();
    present[last_node] = true;
    nodes[last_node] = node;
    queue.push_back(last_node);
}

void local_pool_list::set_node_value(unsigned long long next_work_above) {
    for (size_t i = 0; i < queue.size(); i++)
        nodes[queue[i]].current_node_value = next_work_above;
}

void local_pool_list::sort() {
    if (queue.empty()) return;
    std::sort(
        queue.begin(), queue.end(),
        [&](int a, int b){ return nodes[a].lower_bound > nodes[b].lower_bound; }
    );
}

path_node *local_pool_list::get(int last_node) {
    if (!present[last_node]) return nullptr;
    return &nodes[last_node];
}


void local_pool_thread::initial_depth(int init_depth) {
    lock.lock();
    zero_depth = init_depth;
    depth = init_depth;
    lock.unlock();
}

bool local_pool_thread::pop_from_zero_list(path_node &result_node)
{
    if (level() <= 1) return false;

    lock.lock();

    while (level() > 1 && pool[zero_depth].empty())
    {
        zero_depth++;
    }

    if (level() <= 1)
    {
        lock.unlock();
        return false;
    }

    result_node = pool[zero_depth].back();
    pool[zero_depth].pop_back();

    if (pool[zero_depth].empty())
    {
        zero_depth++;
    }

    lock.unlock();
    return true;
};

bool local_pool_thread::pop_from_active_list(path_node &result_node)
{
    if (level() <= 0 || pool[depth - 1].empty())
        return false;

    result_node = pool[depth - 1].back();
    pool[depth - 1].pop_back();

    return true;
};

void local_pool_thread::start_ready_list()
{
    ready_list().clear();
}

void local_pool_thread::push_ready_list()
{
    lock.lock();
    depth++;
    lock.unlock();
};

void local_pool_thread::pop_active_list()
{
    lock.lock();
    if (level() > 0) {
        depth--;
    }
    lock.unlock();
};

unsigned long long local_pool_thread::node_value()
{
    lock.lock();
    unsigned long long node_value = 0;
    if (level() > 1 && !pool[zero_depth].empty())
    {
        node_value = pool[zero_depth].back().current_node_value;
    }
    lock.unlock();
    return node_value;
}

path_node *local_pool_thread::get(int depth, int last_node)
{
    if (depth > this->depth) return nullptr;
    return pool[depth - 1].get(last_node);
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

void local_pool::print()
{
    for (int i = 0; i < threads.size(); i++)
    {
        std::cout << threads[i].level() << ", ";
    }
    std::cout << std::endl;
}