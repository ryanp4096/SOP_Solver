#ifndef MEMORY_H
#define MEMORY_H

#include <cstdint>

/**
 * @brief Contains system memory statistics.
 */
struct MemoryInfo {
    uint64_t totalBytes;
    uint64_t freeBytes;
    uint64_t availableBytes;
};

/**
 * @brief Retrieves information about the system's physical memory.
 * 
 * @return Struct containing total, free, and available memory.
 */
MemoryInfo getSystemMemory();


/**
 * @brief Get peak resident set size (maximum memory used by process)
 * 
 * @return Peak RSS in bytes
 */
uint64_t getPeakMemoryUsage();

/**
 * @brief Get number of times memory info is looked up
 */
unsigned long long getMemoryLookupCount();
#endif