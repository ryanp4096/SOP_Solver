#ifndef HISTORY_H
#define HISTORY_H

#include <atomic>
#include <cstdint>
#include <mutex>
#include <boost/dynamic_bitset.hpp>
#include "synchronization.hpp"

using namespace std;

struct PrefixKey
{
	boost::dynamic_bitset<> bit_vector;
	int last_node;
};

// struct HistoryContent
// {
// 	int32_t prefix_cost = -1; // the cost of the current path represented by this node
// 	int32_t lower_bound = -1; // the lower bound cost of a solution starting with this path
// };

enum HistoryNodeState : uint8_t {
	UNEXPLORED, EXPLORING, EXPLORED
};

/* A single node in the history table. */
struct HistoryNode
{
	int32_t prefix_cost = -1;
	int32_t lower_bound = -1;
	HistoryNodeState state = UNEXPLORED;
	int8_t active_thread = -1;
	spin_lock lock{};

	// atomic<bool> explored;			 // if the subspace under this node has already been fully explored
	// // atomic<bool> referred;			 // if the subspace under this node has already been referred
	// atomic<bool> is_best_suffix;
	// // int level;
	// atomic<uint8_t> active_threadID; // the thread that is exploring this subspace (there can only ever be one, because any others would be stopped)
	// atomic<HistoryContent> entry;
	// spin_lock lock;
};

struct Active_Node
{
	HistoryNode *history_link = NULL;
	int32_t total_children_cnt = 0; // total number of children of this node
	int32_t cur_children_cnt = 0;	// number of children processed so far
	atomic<int16_t> cur_threadID;
	atomic<bool> deprecated; // whether this node is still necessary to consider
	mutex nlck;				 // lock on total and current children counts
};

struct SubpathKey
{
    boost::dynamic_bitset<> bit_vector;
    int first_node;
    int last_node;
};

struct SubpathHistoryNode
{
    atomic<uint32_t> subpath_cost;
	// uint32_t active_threads;
	// bool in_best_solution;
	spin_lock lock;
};

#endif