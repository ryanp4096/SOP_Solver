#ifndef THREAD_STOPPING_H
#define THREAD_STOPPING_H

#include <atomic>
#include <cstddef>
#include <cstdint>

#include "history_table/history_table.hpp"

constexpr size_t MaxThreads = 32;
constexpr size_t Slots = 64;

/**
 * Packet for a single thread stop request
 */
struct ThreadStopPacket {
    int first_node;
    int last_node;
    int depth;
    Hash128 hash;
};

/**
 * Queue for all requests from a source thread to a target thread
 */
struct alignas(64) ThreadStopQueue {
    alignas(64) std::atomic<uint32_t> head{0};
    alignas(64) std::atomic<uint32_t> tail{0};
    uint32_t cached_tail{0};
    alignas(64) ThreadStopPacket slots[Slots];
};

/**
 * Allows to send stop requests to a target thread, and for target thread to process them
 */
class ThreadStopTarget {
private:
    ThreadStopQueue queues[MaxThreads];
    alignas(64) std::atomic<uint32_t> queue_mask{0};

public:
    /**
     * Send a thread stop request to source thread
     */
    bool request(int source_thread, const ThreadStopPacket &packet) {
        ThreadStopQueue &queue = queues[source_thread];
        uint32_t head = queue.head.load(std::memory_order_relaxed);
        if (head - queue.cached_tail >= Slots) { // if queue full, don't push
            queue.cached_tail = queue.tail.load(std::memory_order_acquire);
            if (head - queue.cached_tail >= Slots) return false;
        }
        queue.slots[head & (Slots - 1)] = packet;
        queue.head.store(head + 1, std::memory_order_release);
        queue_mask.fetch_or(1u << source_thread, std::memory_order_release);
        return true;
    }

    template<class Callback>
    void process(Callback &&callback) {
        uint32_t mask = queue_mask.exchange(0, std::memory_order_acquire);
        while (mask) {
            int i = __builtin_ctz(mask);
            mask &= mask - 1;
            ThreadStopQueue &queue = queues[i];
            uint32_t tail = queue.tail.load(std::memory_order_relaxed);
            uint32_t head = queue.head.load(std::memory_order_acquire);
            for (; tail != head; tail++)
                callback(queue.slots[tail & (Slots - 1)]);
            queue.tail.store(tail, std::memory_order_release);
        }
    }
};

/**
 * Global storage for all thread stopping data
 */
class ThreadStopping {
public:
    ThreadStopTarget targets[MaxThreads];
};

#endif