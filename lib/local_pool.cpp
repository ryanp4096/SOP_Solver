#include "local_pool.hpp"

std::vector<int> local_pool_node_ref::sequence() const {
    std::vector<int> path{};
    path.reserve(depth);
    for (int i = 0; i < depth - 1; i++) {
        path.push_back(thread->current_path[i]);
    }
    path.push_back(taken_node);
    return path;
}

boost::dynamic_bitset<> local_pool_node_ref::bit_vector() const {
    boost::dynamic_bitset<> bit_vector(thread->current_key.size(), false);
    for (int i = 0; i < depth - 1; i++) {
        bit_vector[thread->current_path[i]] = true;
    }
    bit_vector[taken_node] = true;
    return bit_vector;
}

void local_pool_node_ref::to_path_node(path_node &node) const {
    node.sequence = std::move(sequence());
    node.lower_bound = lower_bound();
    node.origin_node = depth > 1 ? node.sequence.at(1) : 0;
    node.current_node_value = current_node_value();
    node.history_key.last_node = taken_node;
    node.history_key.bit_vector = std::move(bit_vector());
}


void local_pool_list::clear() {
    queue.clear();
    present.reset();
}

void local_pool_list::push_back(int taken_node, const local_pool_node &node) {
    present[taken_node] = true;
    nodes[taken_node] = node;
    queue.push_back(taken_node);
}

void local_pool_list::push_back(int taken_node, local_pool_node &&node) {
    present[taken_node] = true;
    nodes[taken_node] = node;
    queue.push_back(taken_node);
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

bool local_pool_list::get(int last_node, local_pool_node_ref &node) {
    if (!present[last_node]) return false;
    node = local_pool_node_ref(&nodes[last_node], thread, depth, last_node);
    return true;
}


void local_pool_thread::initial_state(const std::vector<int> &sequence) {
    lock.lock();
    zero_depth = sequence.size();
    depth = sequence.size();
    current_path = sequence;
    current_key.reset();
    for (size_t i = 0; i < sequence.size(); i++) {
        pool[i].clear();
        current_key[sequence[i]] = true;
    }
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

    pool[zero_depth].back().to_path_node(result_node);
    pool[zero_depth].pop_back();

    if (pool[zero_depth].empty())
    {
        zero_depth++;
    }

    lock.unlock();
    return true;
};

bool local_pool_thread::pop_from_active_list(local_pool_node_ref &result_node)
{
    if (level() <= 0 || pool[depth - 1].empty())
        return false;

    result_node = pool[depth - 1].back();
    pool[depth - 1].pop_back();

    return true;
};

void local_pool_thread::start_ready_list()
{
    pool[depth].clear();
}

void local_pool_thread::push_to_ready_list(int taken_node, int lower_bound, HistoryNode *history_node, int cost)
{
    pool[depth].push_back(taken_node, local_pool_node{
        .lower_bound = lower_bound,
        .history_node = history_node,
        .cost = cost
    });
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
        node_value = pool[zero_depth].back().current_node_value();
    }
    lock.unlock();
    return node_value;
}

bool local_pool_thread::get(int depth, int last_node, local_pool_node_ref &result_node)
{
    if (depth > this->depth) return false;
    return pool[depth - 1].get(last_node, result_node);
}

void local_pool_thread::activate(int taken_node)
{
    current_path.push_back(taken_node);
    current_key[taken_node] = true;
}

void local_pool_thread::deactivate()
{
    current_key[current_path.back()] = false;
    current_path.pop_back();
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
    for (size_t i = 0; i < threads.size(); i++)
    {
        std::cout << threads[i].level() << ", ";
    }
    std::cout << std::endl;
}