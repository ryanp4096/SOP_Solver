#ifndef ALLOCATOR_H
#define ALLOCATOR_H

#include <vector>
#include <boost/dynamic_bitset.hpp>

#define MEMORY_BLOCK_SIZE 5765760

/**
 * Efficiently allocate memory from large pre-allocated blocks.
 * Memory must be freed all at once when done.
 */
template <typename T = unsigned char>
class MemoryPool {
private:
    std::list<T *> blocks;
    size_t block_size;
    size_t pos;

public:
    MemoryPool(size_t block_size = MEMORY_BLOCK_SIZE) : blocks{}, block_size{block_size}, pos{block_size + 1} {}

    ~MemoryPool() { free_all(); }

    T *allocate(size_t length) {
        if (pos + length > block_size) {
            blocks.push_back(new T[block_size]);
            pos = 0;
        }
        T *p = &blocks.back()[pos];
        pos += length;
        return p;
    }

    void free_all() {
        for (T *block : blocks) {
            delete[] block;
        }
        blocks.clear();
    }
};

/**
 * Custom memory allocator.
 * If given a MemoryPool, it will efficiently allocate from large pre-allocated blocks.
 * By default, uses the default malloc/free.
 */
template <typename T>
class MemoryAllocator {
private:
    MemoryPool<> *pool;

public:
    using value_type = T;

    /**
     * Construct a default memory allocator that uses malloc/free
     */
    MemoryAllocator() : pool{} {}

    /**
     * Construct a memory allocator that allocates from a pool of pre-allocated blocks
     */
    MemoryAllocator(MemoryPool<> &pool) : pool{&pool} {}

    template <typename U>
    MemoryAllocator(const MemoryAllocator<U> &allocator) : pool{allocator.pool} {}

    T *allocate(std::size_t n) {
        if (pool == nullptr)
            return static_cast<T *>(malloc(n * sizeof(T)));
        
        return static_cast<T *>(pool->allocate(n * sizeof(T)));
    }

    void deallocate(T *p, std::size_t n) {
        if (pool == nullptr)
            free(p);
    }
};

template <typename T, typename U>
bool operator==(const MemoryAllocator<T>&, const MemoryAllocator<U>&) { return true; }

template <typename T, typename U>
bool operator!=(const MemoryAllocator<T>&, const MemoryAllocator<U>&) { return false; }

template<typename T = unsigned long>
using Bitset = boost::dynamic_bitset<T, MemoryAllocator<T>>;

#endif