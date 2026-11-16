/**
 * Used for expanding buckets when they are full
 */

#ifndef HISTORY_TABLE_EXPANSION_H
#define HISTORY_TABLE_EXPANSION_H

#include <immintrin.h>

#include <cstdint>
#include <cstddef>
#include <vector>

#include "hash.hpp"
#include "buckets.hpp"
#include "memory.hpp"
#include "entries.hpp"
#include <mutex>


class ExpansionPool {
private:
    HistoryExpansion *start{};
    size_t current{};
    size_t capacity{};
public:
    ExpansionPool() = default;

    ExpansionPool(HistoryExpansion *start, size_t capacity)
        : start{start}, current{0}, capacity{capacity} {}

    bool full() { return current >= capacity; }

    HistoryExpansion *allocate_expansion() {
        return &start[current++];
    }
};

class ExpansionManager {
private:
    HistoryExpansion *expansions;               // Pointer to expansion spaces
    size_t current;                             // Next free expansion space
    size_t total;                               // Total available expansion spaces
    std::mutex lock{};                          // Lock so only one thread allocates new pool at a time
    std::vector<ExpansionPool> thread_pools;    // Pool for each thread

public:
    ExpansionManager(HistoryExpansion *expansions, size_t expansion_count, size_t thread_count)
          : expansions{expansions},
            current{0},
            total{expansion_count},
            lock{},
            thread_pools(thread_count)
    {
        for (size_t i = 0; i < thread_count; i++) create_pool(i);
    }

    ExpansionManager(const ExpansionManager &other)
          : expansions{other.expansions},
            current{other.current},
            total{other.total},
            lock{},
            thread_pools{other.thread_pools} {}
    
    ExpansionManager &operator=(const ExpansionManager &other) {
        expansions = other.expansions;
        current = other.current;
        total = other.total;
        thread_pools = other.thread_pools;
        return *this;
    }

    HistoryExpansion *allocate_expansion(size_t thread_id) {
        if (thread_pools[thread_id].full()) {
            bool success = create_pool(thread_id);
            if (!success) return nullptr;
        }
        return thread_pools[thread_id].allocate_expansion();
    }

private:
    bool create_pool(size_t thread_id) {
        if (current + expansion_pool_size > total) return false;

        lock.lock();
        if (current + expansion_pool_size > total) {
            lock.unlock();
            return false;
        }

        thread_pools[thread_id] = ExpansionPool(
            &expansions[current], expansion_pool_size
        );
        current += expansion_pool_size;

        lock.unlock();
        return true;
    }
};

#endif