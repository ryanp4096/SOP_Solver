#ifndef HISTORY_TABLE_LOCK_H
#define HISTORY_TABLE_LOCK_H

#include <immintrin.h>

#include <cstdint>
#include <cstddef>

/**
 * Lock for history table buckets.
 * Locks when constructed, releases when destroyed.
 */
class BucketLock {
private:
    uint8_t* lock_byte;

public:
    /**
     * Locks bucket. Releases on deletion
     */
    explicit BucketLock(uint8_t* lock_byte) : lock_byte(lock_byte) {
        for (;;) {
            // Try to take the lock with one atomic exchange
            if (__atomic_exchange_n(lock_byte, static_cast<uint8_t>(1), __ATOMIC_ACQUIRE) == 0)
                return;
            
            // Someone else holds it. Spin with plain reads, which keep the cache
            // line shared instead of fighting over ownership, until it looks free
            while (__atomic_load_n(lock_byte, __ATOMIC_RELAXED) != 0)
                _mm_pause();
        }
    }

    /**
     * When destroyed, unlock
     */
    ~BucketLock() {
        __atomic_store_n(lock_byte, static_cast<uint8_t>(0), __ATOMIC_RELEASE);
    }
};

/**
 * Lock for history table buckets.
 * Locks two buckets at once for history table format storing entry in smallest of 2 buckets.
 */
class BucketDoubleLock {
private:
    uint8_t *lower_lock;
    uint8_t *higher_lock;
    bool higher_locked = false;

public:
    /**
     * Creates lock for two buckets and locks first bucket
     * To avoid deadlock, lower index lock must be locked first.
     */
    explicit BucketDoubleLock(uint8_t *lock_byte1, uint8_t *lock_byte2)
        : lower_lock(lock_byte1), higher_lock(lock_byte2) {
        acquire(lower_lock);
    }

    /**
     * Locks second bucket
     */
    void lock_second() {
        acquire(higher_lock);
        higher_locked = true;
    }

    /**
     * Unlocks both buckets on destruction
     */
    ~BucketDoubleLock() {
        if (higher_locked) release(higher_lock);
        release(lower_lock);
    }

private:
    /* Acquire a lock */
    static void acquire(uint8_t *lock_byte) {
        for (;;) {
            // Try to take the lock with one atomic exchange
            if (__atomic_exchange_n(lock_byte, static_cast<uint8_t>(1), __ATOMIC_ACQUIRE) == 0)
                return;
            
            // Someone else holds it. Spin with plain reads, which keep the cache
            // line shared instead of fighting over ownership, until it looks free
            while (__atomic_load_n(lock_byte, __ATOMIC_RELAXED) != 0)
                _mm_pause();
        }
    }

    /* Release a lock */
    static void release(uint8_t *lock_byte) {
        __atomic_store_n(lock_byte, static_cast<uint8_t>(0), __ATOMIC_RELEASE);
    }
};

#endif