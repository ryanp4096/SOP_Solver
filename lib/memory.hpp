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

#endif