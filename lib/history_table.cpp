#include "history_table.hpp"

#define MEMORY_RESTRIC 0.95
#define SAMPLE_FREQUENCY 60
#include <condition_variable>
#include <mutex>

static bool is_working = false;
int gp_depth = 0;
std::mutex mtx;
std::condition_variable cv;

template <typename T>
MemoryAllocator<T>::~MemoryAllocator() {
    free_all();
}
template <typename T>
T *MemoryAllocator<T>::allocate() {
    if (blocks.empty() || items_count >= items_per_block) {
        T *block = new T[items_per_block];
        blocks.push_back(block);
        items_count = 1;
        return block;
    } else {
        return &blocks.back()[items_count++];
    }
}

template<typename T>
void MemoryAllocator<T>::free_all() {
    for (T *b : blocks) {
        delete[] b;
    }
    blocks.clear();
}

// Memory_Module::Memory_Module()
// {
//     bucket_block = new Bucket[BUCKET_BLK_SIZE];
//     history_block = (HistoryNode *)malloc(HIS_BLK_SIZE * sizeof(HistoryNode));
//     bucket_counter = 0;
//     his_node_counter = 0;
// }
// Memory_Module::~Memory_Module()
// {
//     // cout << "destructor triggered\n";
//     delete[] bucket_block; // Use delete[] to free the array of Buckets
//     free(history_block);   // Use free to deallocate memory allocated with malloc
// }

// Bucket *Memory_Module::get_bucket()
// {
//     if (bucket_counter >= BUCKET_BLK_SIZE || bucket_block == NULL)
//     {
//         bucket_block = new Bucket[BUCKET_BLK_SIZE];
//         bucket_counter = 0;
//     }
//     Bucket *bucket = bucket_block + bucket_counter;
//     bucket_counter++;
//     return bucket;
// }

// HistoryNode *Memory_Module::retrieve_his_node()
// {
//     // HistoryNode* node = NULL;

//     if (his_node_counter >= HIS_BLK_SIZE || history_block == NULL)
//     {
//         history_block = (HistoryNode *)malloc(HIS_BLK_SIZE * sizeof(HistoryNode));
//         his_node_counter = 0;
//     }
//     HistoryNode *node = history_block + his_node_counter;
//     his_node_counter++;

//     return node;
// }

History_Table::History_Table(size_t size)
{
    // struct sysinfo info;
    // if (sysinfo(&info) != 0)
    // {
    //     cout << "can't retrieve system memory info\n";
    //     exit(EXIT_FAILURE);
    // }
    MemoryInfo info = getSystemMemory();

    num_buckets = size;
    total_ram = info.totalBytes;
    max_size = (double)info.availableBytes * MEMORY_RESTRIC - (size * 4);
    cout << "Total RAM (in bytes): " << total_ram / 1048576 << " MB" << std::endl;
    cout << "allowed RAM size (in bytes): " << max_size / 1048576 << " MB" << std::endl;
    current_size = 0;

    // std::cout << "Max bucket size is " << max_size / 1048576 << " MB" << std::endl;
    // std::cout << "Total Available Memory In OS is " << info.totalram / 1048576 << " MB" << std::endl;

    insert_count = 0;
}

void History_Table::initialize(int thread_num, size_t size, int number_of_groups, int group_size, timer *main_timer, bool enable_subpath_history_table)
{
    this->main_timer = main_timer;
    this->enable_subpath_history_table = enable_subpath_history_table;
    num_of_groups = number_of_groups;
    groups_size = group_size;
    block_count.resize(number_of_groups, 0);

    prefix_maps.resize(number_of_groups);

    if (enable_subpath_history_table) {
        subpath_maps.resize(number_of_groups);
    }

    blocked_groups.resize(number_of_groups, false);
    is_data_available.resize(number_of_groups, true);

    group_locks = vector<spin_lock>(number_of_groups);

    for (int i = 0; i < number_of_groups; i++)
    {
        prefix_maps[i].bucket_allocators.resize(thread_num, MemoryAllocator<PrefixEntry>(BUCKET_BLK_SIZE));
        prefix_maps[i].locks = vector<spin_lock>(size / COVER_AREA + 1);
        prefix_maps[i].buckets = vector<PrefixBucket>(size);

        if (enable_subpath_history_table) {
            subpath_maps[i].bucket_allocators.resize(thread_num, MemoryAllocator<SubpathEntry>(BUCKET_BLK_SIZE));
            subpath_maps[i].locks = vector<spin_lock>(size / COVER_AREA + 1);
            subpath_maps[i].buckets = vector<SubpathBucket>(size);
        }
    }
}

size_t History_Table::get_max_size() { return max_size; }
size_t History_Table::get_current_size() { return current_size; }

unsigned long History_Table::get_free_mem()
{
    // struct sysinfo info;
    // if (sysinfo(&info) != 0)
    // {
    //     cout << "can't retrieve sys mem info\n";
    //     exit(1);
    // }
    MemoryInfo info = getSystemMemory();
    return (double)info.availableBytes;
}

void History_Table::print_curmem()
{
    std::cout << "Current memory is " << current_size / 1048576 << " MB" << std::endl;
    return;
}

HistoryNode *History_Table::insert(PrefixKey &key, unsigned int depth, int prefix_cost, int lower_bound, HistoryNodeState state, unsigned int thread_id)
{
    int group_index = get_bucket_index(depth);

    if (blocked_groups[group_index])
        return NULL;
    
    PrefixMap &map = prefix_maps[group_index];

    size_t val = hash<boost::dynamic_bitset<>>{}(key.bit_vector);
    int bucket_index = (val + key.last_node * 1217) % num_buckets;

    spin_lock &lock = map.locks[bucket_index / COVER_AREA];
    lock.lock();

    PrefixEntry *entry = insert_prefix_entry(map, group_index, thread_id, bucket_index, key, prefix_cost, lower_bound, state);
    lock.unlock();
    return entry == NULL ? NULL : &entry->node;
}

HistoryNode *History_Table::retrieve_or_insert(PrefixKey &key, unsigned int depth, int prefix_cost, int lower_bound, HistoryNodeState state, unsigned thread_id, bool *inserted)
{
    int group_index = get_bucket_index(depth);
    *inserted = false;
    if (!is_data_available[group_index]) return NULL;
    if (blocked_groups[group_index]) return retrieve(key, depth);
    PrefixMap &map = prefix_maps[group_index];
    
    size_t val = hash<boost::dynamic_bitset<>>{}(key.bit_vector);
    int bucket_index = (val + key.last_node * 1217) % num_buckets;

    spin_lock &lock = map.locks[bucket_index / COVER_AREA];
    lock.lock();

    PrefixEntry *found_entry = search_prefix_bucket(map.buckets[bucket_index], key);
    if (found_entry != NULL) {
        lock.unlock();
        return &found_entry->node;
    }

    PrefixEntry *entry = insert_prefix_entry(map, group_index, thread_id, bucket_index, key, prefix_cost, lower_bound, state);
    lock.unlock();
    *inserted = true;
    return entry == NULL ? NULL : &entry->node;
}

HistoryNode *History_Table::retrieve(PrefixKey &key, unsigned int depth)
{
    if (depth < gp_depth) return NULL;
    int group_index = get_bucket_index(depth);
    if (!is_data_available[group_index]) return NULL;
    PrefixMap &map = prefix_maps[group_index];

    size_t val = hash<boost::dynamic_bitset<>>{}(key.bit_vector);
    int bucket_index = (val + key.last_node * 1217) % num_buckets;

    PrefixEntry *entry = search_prefix_bucket(map.buckets[bucket_index], key);
    if (entry == NULL) return NULL;
    return &entry->node;
}


SubpathHistoryNode *History_Table::insert_subpath(SubpathKey &key, unsigned int depth, int subpath_cost, unsigned int thread_id)
{
    if (!enable_subpath_history_table) return NULL;
    int group_index = get_bucket_index(depth);

    if (blocked_groups[group_index])
        return NULL;
    
    SubpathMap &map = subpath_maps[group_index];

    size_t val = hash<boost::dynamic_bitset<>>{}(key.bit_vector);
    int bucket_index = (val + key.first_node * 1583 + key.last_node * 1217) % num_buckets;

    spin_lock &lock = map.locks[bucket_index / COVER_AREA];
    lock.lock();

    SubpathEntry *entry = insert_subpath_entry(map, group_index, thread_id, bucket_index, key, subpath_cost);
    lock.unlock();
    return entry == NULL ? NULL : &entry->node;
}


SubpathHistoryNode *History_Table::retrieve_or_insert_subpath(SubpathKey &key, unsigned int depth, int subpath_cost, unsigned int thread_id, bool *inserted)
{
    if (!enable_subpath_history_table) return NULL;
    int group_index = get_bucket_index(depth);
    *inserted = false;
    if (!is_data_available[group_index]) return NULL;
    if (blocked_groups[group_index]) return retrieve_subpath(key, depth);
    SubpathMap &map = subpath_maps[group_index];
    
    size_t val = hash<boost::dynamic_bitset<>>{}(key.bit_vector);
    int bucket_index = (val + key.first_node * 1583 + key.last_node * 1217) % num_buckets;

    spin_lock &lock = map.locks[bucket_index / COVER_AREA];
    lock.lock();

    SubpathEntry *found_entry = search_subpath_bucket(map.buckets[bucket_index], key);
    if (found_entry != NULL) {
        lock.unlock();
        return &found_entry->node;
    }

    SubpathEntry *entry = insert_subpath_entry(map, group_index, thread_id, bucket_index, key, subpath_cost);
    lock.unlock();
    *inserted = true;
    return entry == NULL ? NULL : &entry->node;
}

SubpathHistoryNode *History_Table::retrieve_subpath(SubpathKey &key, unsigned int depth)
{
    if (!enable_subpath_history_table) return NULL;
    int group_index = get_bucket_index(depth);
    if (!is_data_available[group_index]) return NULL;
    SubpathMap &map = subpath_maps[group_index];

    size_t val = hash<boost::dynamic_bitset<>>{}(key.bit_vector);
    int bucket_index = (val + key.first_node * 1583 + key.last_node * 1217) % num_buckets;

    SubpathEntry *entry = search_subpath_bucket(map.buckets[bucket_index], key);
    if (entry == NULL) return NULL;
    return &entry->node;
}

PrefixEntry *History_Table::search_prefix_bucket(PrefixBucket &bucket, PrefixKey &key)
{
    PrefixEntry *first = bucket.load();
    if (first == NULL) return NULL;

    if (first->next == NULL) {
        if (
            key.last_node == first->key.last_node &&
            key.bit_vector == first->key.bit_vector
        ) {
            return first;
        } else {
            return NULL;
        }
    }
    for (PrefixEntry *entry = first; entry != NULL; entry = entry->next) {
        if (
            key.last_node == entry->key.last_node &&
            key.bit_vector == entry->key.bit_vector
        ) {
            return entry;
        }
    }
    return NULL;
}

SubpathEntry *History_Table::search_subpath_bucket(SubpathBucket &bucket, SubpathKey &key)
{
    SubpathEntry *first = bucket.load();
    if (first == NULL) return NULL;

    if (first->next == NULL) {
        if (
            key.first_node == first->key.first_node &&
            key.last_node == first->key.last_node &&
            key.bit_vector == first->key.bit_vector
        ) {
            return first;
        } else {
            return NULL;
        }
    }
    for (SubpathEntry *entry = first; entry != NULL; entry = entry->next) {
        if (
            key.first_node == entry->key.first_node &&
            key.last_node == entry->key.last_node &&
            key.bit_vector == entry->key.bit_vector
        ) {
            return entry;
        }
    }
    return NULL;
}

PrefixEntry *History_Table::insert_prefix_entry(PrefixMap &map, int group_index, unsigned int thread_id, size_t bucket_index, PrefixKey &key, int prefix_cost, int lower_bound, HistoryNodeState state)
{
    if (limit_insertion) return NULL;
    if (current_size >= max_size) {
        limit_insertion = true;
        double t = main_timer->get_time_seconds();
        std::cout << "Blocking insertion at time " << t << std::endl;
        return NULL;
    }

    if (thread_id % 4 == 0)
    {
        insert_count++;
        if (insert_count >= 100000)
        {
            current_size = total_ram - get_free_mem();
            insert_count = 0;
        }
    }

    if (!is_data_available[group_index] || blocked_groups[group_index])
    {
        block_count[group_index]++;
        return NULL;
    }

    PrefixEntry *entry = map.bucket_allocators[thread_id].allocate();
    if (entry == NULL)
        return NULL;
    
    entry->key = key;
    entry->node.prefix_cost = prefix_cost;
    entry->node.lower_bound = lower_bound;
    entry->node.state = state;
    entry->node.active_thread = thread_id;
    entry->next = map.buckets[bucket_index];
    map.buckets[bucket_index] = entry;
    return entry;
}

SubpathEntry *History_Table::insert_subpath_entry(SubpathMap &map, int group_index, unsigned int thread_id, size_t bucket_index, SubpathKey &key, int subpath_cost)
{
    if (limit_insertion) return NULL;
    if (current_size >= max_size) {
        limit_insertion = true;
        double t = main_timer->get_time_seconds();
        std::cout << "Blocking insertion at time " << t << std::endl;
        return NULL;
    }

    if (thread_id % 4 == 0)
    {
        insert_count++;
        if (insert_count >= 100000)
        {
            current_size = total_ram - get_free_mem();
            insert_count = 0;
        }
    }

    if (!is_data_available[group_index] || blocked_groups[group_index])
    {
        block_count[group_index]++;
        return NULL;
    }

    SubpathEntry *entry = map.bucket_allocators[thread_id].allocate();
    if (entry == NULL)
        return NULL;
    
    entry->key = key;
    entry->node.subpath_cost = subpath_cost;
    entry->next = map.buckets[bucket_index];
    map.buckets[bucket_index] = entry;
    return entry;
}


bool History_Table::check_and_manage_memory(int depth, float *updated_mem_limit, bool *is_all_table_blocked)
{
    int group_index = get_bucket_index(depth);

    for (int i = prefix_maps.size() - 1; i > 0; --i)
    {
        current_size = total_ram - get_free_mem();
        // cout << "attempted to block the group :" << group_index << "\n";
        if (current_size < *updated_mem_limit * max_size)
            return blocked_groups[group_index];

        if (blocked_groups[i]) // to prevent using unnecessary locks
            continue;

        group_locks[i].lock();
        if (current_size < *updated_mem_limit * max_size)
        {
            group_locks[i].unlock();
            return blocked_groups[group_index];
            // return group_index > i; // if the subtable number is lesser than i, we CAN insert that record by returning FALSE
        }
        if (!blocked_groups[i])
        {
            cout << "Blocking insertion in bucket: " << i + 1 << " when the mem limit is: " << *updated_mem_limit << std::endl;
            blocked_groups[i] = true;
            *is_all_table_blocked = (i == 1);

            /**
             * when the global pool is empty, we are updating the memory limit to 0.9
             * and to prevent updating the memory limit to greater than 0.9
             *
             * mem_limit was 0.7, so we are updating limit to 0.8
             * global pool is empty, so we are updating the memory limit to 0.9 (most recent limit)
             * max(limit, updated_mem_limit) = will return the most updated memory limit
             * min(0.9, max(limit, updated_mem_limit)) = will make sure that limit should not exceed 0.9
             */
            float limit = *updated_mem_limit + 0.1f;
            *updated_mem_limit = min(0.9f, max(limit, *updated_mem_limit));
            cout << "Updated memory limit: " << *updated_mem_limit << std::endl;

            // for (int j = memory_allocators.size() - 1; j >= 0; --j)
            // {
            //     cout << blocked_groups[j] << " ";
            // }
            cout << "\n";
            group_locks[i].unlock();
            return blocked_groups[i]; // if the subtable number is equal to the blocking group, we CAN"T insert that record by returning TRUE
        }
        group_locks[i].unlock();
    }
    return blocked_groups[group_index];
}

bool History_Table::free_subtable_memory(float *mem_limit)
{
    std::unique_lock<std::mutex> lock(mtx); // Lock for the entire function to synchronize threads

    for (int i = prefix_maps.size() - 1; i >= 0; --i)
    {
        // Wait until no thread is working
        cv.wait(lock, []
                { return !is_working; });

        current_size = total_ram - get_free_mem();

        /**
         * if i == 0, we don't have anything left to clear out
         * but to block the last bucket
         */
        if (current_size >= *mem_limit * max_size && i == 0)
        {
            /**
             * when our global pool is empty and our history_table size is 3
             * we are calling this function to block the last bucket
             * since our history_pool size is 3, we are blocking the last bucket only
             * we need to block all of them, that's why we are running a loop again to make sure all the groups are blocked
             */
            for (int j = blocked_groups.size() - 1; j >= 0; --j)
                blocked_groups[j] = true;

            print_curmem();
            // if i == 0 is true, we don't have anything left to clear out except
            // the last remaining bucket that's why we always return false here

            return false;
        }
        if (current_size < *mem_limit * max_size)
            return true;

        if (blocked_groups[i] && is_data_available[i]) // Ensure the group is blocked and data is there
        {

            // cout << "locking the group\n";
            group_locks[i].lock();
            is_working = true;
            lock.unlock(); // Unlock the main mutex to allow other threads to enter wait state

            current_size = total_ram - get_free_mem();
            if (current_size < *mem_limit * max_size)
            {
                is_working = false;
                group_locks[i].unlock();
                cv.notify_all(); // Notify other waiting threads
                return true;
            }

            // Free the memory if group entry is blocked and data is available
            if (blocked_groups[i] && is_data_available[i])
            {
                is_data_available[i] = false;
                cout << "starting the deallocation of memory\n";
                current_size = total_ram - get_free_mem();

                std::cout << "current used size before (in bytes): " << current_size / 1048576 << " MB" << std::endl;

                for (MemoryAllocator<PrefixEntry> &allocator : prefix_maps[i].bucket_allocators)
                {
                    allocator.free_all();
                }

                std::cout << "Freed memory for subtable: " << i + 1 << std::endl;

                current_size = total_ram - get_free_mem();
                std::cout << "current used size after freeing space (in bytes): " << current_size / 1048576 << " MB" << std::endl;
                // for (int k = 0; k < 3; k++)
                // {
                //     cout << block_count[k] << " ";
                // }
                cout << "\n";
                is_working = false;
                group_locks[i].unlock();
                cv.notify_all(); // Notify other waiting threads
                return true;
            }
            group_locks[i].unlock();
            is_working = false;
            cv.notify_all(); // Notify other waiting threads
        }
    }
    return (current_size < *mem_limit * max_size);
}
void History_Table::track_entries_and_references()
{
    long total_entries = 0;
    long total_references = 0;

    vector<long long> total_entries_single_table(3, 0);
    vector<long long> total_references_single_table(3, 0);

    for (int i = 0; i < num_of_groups; ++i)
    {
        total_entries = 0;
        total_references = 0;

        if (num_of_groups == 1 || is_data_available[i])
            for (PrefixEntry *bucket : prefix_maps[i].buckets)
            {
                if (!bucket)
                    continue; // Skip if the bucket pointer is NULL

                // Iterate over each Entry in the bucket
                for (PrefixEntry *entry = bucket; entry != NULL; entry = entry->next)
                {
                    HistoryNode *history_node = &entry->node;

                    // Check if history_node is not nullptr before accessing its members
                    if (!history_node)
                        continue;

                    // int index = history_node->level;
                    // if (index > 2)
                    //     cout << "incorrect index";

                    // total_entries_single_table[index]++;
                    // if (history_node->referred)
                    //     total_references_single_table[index]++;
                }
            }
    }
    for (int i = 0; i < 3; ++i)
    {
        cout << "Total Entries Single Table " << i << ": " << total_entries_single_table[i] << endl;
        cout << "Total References Single Table " << i << ": " << total_references_single_table[i] << endl;
    }
}

int History_Table::get_bucket_index(int depth)
{
    if (num_of_groups <= 1)
        return 0;
    if (depth <= gp_depth)
        return 0;
    else if (depth <= num_of_groups * groups_size + gp_depth)
        return std::ceil(static_cast<double>(depth - gp_depth) / groups_size) - 1;
    else
        return num_of_groups - 1;
}

void History_Table::update_gp_depth(int depth)
{
    cout << "Global pool entry depth = " << depth << endl;
    gp_depth = depth;
}