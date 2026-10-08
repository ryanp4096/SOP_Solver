#ifndef HASH_H
#define HASH_H

#include <iostream>
#include <list>

#include "memory.hpp"
#include "timer.hpp"

// #include <boost/smart_ptr/detail/spinlock.hpp>
#include <boost/dynamic_bitset.hpp>

#include "history.hpp"
#include "synchronization.hpp"

#define BUCKET_BLK_SIZE 81920
#define HIS_BLK_SIZE 81920
#define FREED_SIZE 100000000
#define COVER_AREA 10

/**
 * Allocates memory for a given type. Used when frequently allocating memory for a given type to prevent high number of memory allocations, by allocating one large block of memory at a time.
 */
template <typename T>
class MemoryAllocator {
private:
    vector<T *> blocks{};
    unsigned int items_count = 0;
    unsigned int items_per_block;
public:
    MemoryAllocator(unsigned int items_per_block = HIS_BLK_SIZE)
        : items_per_block(items_per_block) {}
    ~MemoryAllocator();
    /**
     * Allocate memory to store type T
     */
    T *allocate();

    void free_all();
};

/* A node in the linked list of entries in a bucket */
struct PrefixEntry {
    PrefixKey key;
    HistoryNode node;
    PrefixEntry *next;
};

/* A bucket in the prefix history table is an atomic pointer to the first entry in the linked list */
typedef atomic<PrefixEntry *> PrefixBucket;

/* A node in the linked list of entries in a bucket */
struct SubpathEntry {
    SubpathKey key;
    SubpathHistoryNode node;
    SubpathEntry *next;
};

/* A bucket in the subpath history table is an atomic pointer to the first entry in the linked list */
typedef atomic<SubpathEntry *> SubpathBucket;

/* Contains everything needed for retrieving and allocating to buckets of prefix entries */
struct PrefixMap {
    vector<PrefixBucket> buckets;                               // Buckets of entries
    vector<spin_lock> locks;                                    // A lock for every COVER_AREA buckets
    vector<MemoryAllocator<PrefixEntry>> bucket_allocators;     // A memory allocator for every thread
};

/* Contains everything needed for retrieving and allocating to buckets of subpath entries */
struct SubpathMap {
    vector<SubpathBucket> buckets;                              // Buckets of entries
    vector<spin_lock> locks;                                    // A lock for every COVER_AREA buckets
    vector<MemoryAllocator<SubpathEntry>> bucket_allocators;    // A memory allocator for every thread
};

/* Separate subpath table containing only LKH subpaths for faster access if lkh_subpaths_only is enabled */
class LKH_Subpath_Table {
private:
    struct LKHSubpathEntry {
        SubpathKey key;
        SubpathHistoryNode node;
    };

    bool data_available = false;
    MemoryAllocator<LKHSubpathEntry> allocator{};
    vector<vector<LKHSubpathEntry *>> data{};

public:
    void initialize(int instance_size) {
        data.resize(instance_size);
        for (int i = 0; i < instance_size; i++) {
            data[i].resize(instance_size, nullptr);
        }
    }

    SubpathHistoryNode *insert(SubpathKey &key, unsigned int length, int subpath_cost) {
        LKHSubpathEntry *entry = allocator.allocate();
        entry->key = key;
        entry->node.subpath_cost = subpath_cost;
        data[length][key.last_node] = entry;
        return &entry->node;
    }

    void complete_insertion() {
        data_available = true;
    }

    SubpathHistoryNode *retrieve(SubpathKey &key, unsigned int length, bool *can_break = NULL) {
        if (can_break != NULL) *can_break = false;
        if (!data_available) {
            if (can_break != NULL) *can_break = true;
            return NULL;
        }
        LKHSubpathEntry *entry = data[length][key.last_node];
        if (entry == nullptr) {
            if (can_break != NULL) *can_break = true;
            return NULL;
        }
        if (entry->key.first_node != key.first_node) return NULL;
        if (entry->key.bit_vector != key.bit_vector) return NULL;
        return &entry->node;
    }
};

/* A thread-safe collection of history entries for previously processed partial paths.
    For efficiency, the history table is allocated in buckets of many nodes, not individually. */
class History_Table
{
public:
    enum SubpathHistorySetting {
        SUBPATHS_OFF,           // do not store subpaths in the history table
        SUBPATHS_LKH_ONLY,      // only store subpaths of lkh's best tour in the history table. stored in a different format to speed up checking for matches
        SUBPATHS_ON             // store all subpaths in the history table
    };
private:
    SubpathHistorySetting subpath_history_setting;

    size_t num_buckets = 0;              // The number of buckets the history table should be stored in
    vector<PrefixMap> prefix_maps;       // A table of buckets, allocators, and locks for each group
    vector<SubpathMap> subpath_maps;     // A table of buckets, allocators, and locks for each group
    LKH_Subpath_Table lkh_subpath_table; // Special table only storing lkh best tour subpaths for faster access

    unsigned long total_ram = 0;        // The total amount of memory in the system, in bytes
    unsigned long max_size = 0;         // The maximum allowed size of the history table, in bytes
    atomic<unsigned long> current_size; // The current size of the history table, in bytes
    int insert_count;                   // A counter to ensure that, periodically, current_size is updated

    int num_of_groups;                  // The number of groups in history table
    int groups_size;                    // The node depth for insertion for each group

    vector<bool> blocked_groups;        // Track which groups are blocked from insertions
    vector<bool> is_data_available;     // Tracks whether data is available in each subtable
    vector<spin_lock> group_locks;      // Locks for each group for modifying an entire group

    vector<long> block_count;

    bool limit_insertion = false;       // Whether insertion is blocked as the history table is full
    timer *main_timer;

public:
    /**
     * @brief Construct a new history table object
     * @param size The number of buckets the history table should be stored in
     */
    History_Table(size_t size);

    /**
     * @brief Initialize the history table with various settings
     * @param thread_count The total number of threads
     * @param size The number of buckets the history table should be stored in
     * @param number_of_groups The number of groups 
     * @param group_size The node depth for insertion for each group
     * @param main_timer Pointer to the solver's timer
     * @param instance_size The size (number of vertices) of the instance
     * @param subpath_history_setting Whether to store no subpaths, lkh subpaths only, or all subpaths
     */
    void initialize(int thread_count, size_t size, int number_of_groups, int group_size, timer *main_timer, unsigned int instance_size, SubpathHistorySetting subpath_history_setting = SUBPATHS_OFF);
    
    /* Returns the max allowed size of the history table, in bytes. */
    size_t get_max_size();

    /* Returns the current size of the history table, in bytes. */
    size_t get_current_size();

    /* Calculates the amount of free ram, in bytes, available on system. */
    unsigned long get_free_mem();

    /* Print to console the total amount of memory exhausted so far. */
    void print_curmem();

    /**
     * @brief Insert a prefix entry into the history table
     * @param key The key specifying the bitset and last node of the prefix
     * @param depth The length of the prefix
     * @param prefix_cost The cost of the prefix
     * @param lower_bound The calculated value of the lower bound at the prefix
     * @param state Whether the node has been explored, is being explored, or is unexplored
     * @param thread_id The id of the thread inserting this entry
     * @return Pointer to the inserted history node, or NULL if table is full
     */
    HistoryNode *insert(PrefixKey &key, unsigned int depth, int prefix_cost, int lower_bound, HistoryNodeState state, unsigned int thread_id);
    
    /**
     * @brief Retrieve a prefix entry corresponding to a key
     * @param key The key specifying the bitset and last node of the prefix
     * @param depth The length of the prefix
     * @return Pointer to the found history node, or NULL if not found
     */
    HistoryNode *retrieve(PrefixKey &key, unsigned int depth);

    /**
     * @brief Try to retrieve a prefix entry and insert if not found
     * @param key The key specifying the bitset and last node of the prefix
     * @param depth The length of the prefix
     * @param prefix_cost The cost of the prefix
     * @param lower_bound The calculated value of the lower bound at the prefix
     * @param state Whether the node has been explored, is being explored, or is unexplored
     * @param thread_id The id of the thread inserting this entry
     * @param inserted A return variable containing true if an entry was not found and a new entry was created
     * @return Pointer to the inserted history node, or NULL if table is full
     */
    HistoryNode *retrieve_or_insert(PrefixKey &key, unsigned int depth, int prefix_cost, int lower_bound, HistoryNodeState state, unsigned thread_id, bool *inserted);

    /**
     * @brief Insert a subpath entry into the history table
     * @param key The key specifying the bitset and first and last nodes of the subpath
     * @param depth The length of the subpath
     * @param subpath_cost The cost of the subpath
     * @param thread_id The id of the thread inserting this entry
     * @return Pointer to the inserted history node, or NULL if table is full
     */
    SubpathHistoryNode *insert_subpath(SubpathKey &key, unsigned int depth, int subpath_cost, unsigned int thread_id);

    /**
     * @brief Retrieve a subpath entry corresponding to a key
     * @param key The key specifying the bitset and first and last nodes of the subpath
     * @param depth The length of the subpath
     * @param can_break Return variable. If lkh_subpaths_only is enabled, set to true if no other subpaths that contain this subpath can be matched
     * @return Pointer to the found history node, or NULL if not found
     */
    SubpathHistoryNode *retrieve_subpath(SubpathKey &key, unsigned int depth, bool *can_break = NULL);

    /**
     * @brief Try to retrieve a subpath entry and insert if not found
     * @param key The key specifying the bitset and first and last nodes of the subpath
     * @param depth The length of the subpath
     * @param subpath_cost The cost of the subpath
     * @param thread_id The id of the thread inserting this entry
     * @param inserted A return variable containing true if an entry was not found and a new entry was created
     * @param can_break Return variable. If lkh_subpaths_only is enabled, set to true if no other subpaths that contain this subpath can be matched
     * @return Pointer to the inserted history node, or NULL if table is full
     */
    SubpathHistoryNode *retrieve_or_insert_subpath(SubpathKey &key, unsigned int depth, int subpath_cost, unsigned int thread_id, bool *inserted, bool *can_break = NULL);

    /* Checks for full groups and blocks them, if there are multiple groups */
    bool check_and_manage_memory(int depth, float *updatedMemLimit, bool *is_all_table_blocked);

    /* Frees a group if out of memory, if there are multiple groups */
    bool free_subtable_memory(float *mem_limit);

    /* To track down the entries and its reference in history_table */
    void track_entries_and_references();

    /* Get the index of the group containing entries of a certain depth, if there are multiple groups */
    int get_bucket_index(int depth);

    /* Update the depth of entries in the global pool */
    void update_gp_depth(int gp_depth);

    /* If lkh_subpaths_only enabled, registers that all lkh subpaths have been inserted */
    void complete_lkh_subpath_insertion() { lkh_subpath_table.complete_insertion(); }

private:
    /**
     * @brief Search a prefix bucket for entries matching a key
     * @param bucket The bucket to search
     * @param key The key to search for
     * @return A pointer to the found entry, or NULL if not found
     */
    PrefixEntry *search_prefix_bucket(PrefixBucket &bucket, PrefixKey &key);

    /**
     * @brief Search a subpath bucket for entries matching a key
     * @param bucket The bucket to search
     * @param key The key to search for
     * @return A pointer to the found entry, or NULL if not found
     */
    SubpathEntry *search_subpath_bucket(SubpathBucket &bucket, SubpathKey &key);

    /**
     * @brief Insert an entry into a prefix bucket
     * @param map The map of the group to insert into
     * @param group_index The index of the group to insert into
     * @param thread_id The id of the thread inserting the entry
     * @param bucket_index The index of the bucket to insert into
     * @param key The key to be stored in the inserted entry
     * @param prefix_cost The cost to be stored in the inserted entry
     * @param lower_bound The lower bound to be stored in the inserted entry
     * @param state The state to be stored in the inserted entry
     * @return A pointer to the inserted entry, or NULL if could not be inserted
     */
    PrefixEntry *insert_prefix_entry(PrefixMap &map, int group_index, unsigned int thread_id, size_t bucket_index, PrefixKey &key, int prefix_cost, int lower_bound, HistoryNodeState state);
    
    /**
     * @brief Insert an entry into a subpath bucket
     * @param map The map of the group to insert into
     * @param group_index The index of the group to insert into
     * @param thread_id The id of the thread inserting the entry
     * @param bucket_index The index of the bucket to insert into
     * @param key The key to be stored in the inserted entry
     * @param subpath_cost The cost to be stored in the inserted entry
     * @return A pointer to the inserted entry, or NULL if could not be inserted
     */
    SubpathEntry *insert_subpath_entry(SubpathMap &map, int group_index, unsigned int thread_id, size_t bucket_index, SubpathKey &key, int subpath_cost);
};

#endif