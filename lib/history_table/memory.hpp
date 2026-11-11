/**
 * Allocates one large block of zero-initialized memory for the history table
 */

#ifndef HISTORY_MEMORY_H
#define HISTORY_MEMORY_H

#include <dirent.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "../memory.hpp"
#include <iostream>


/* Size of huge pages on x86-64 in bytes */
constexpr size_t huge_page_bytes = 2 * 1024 * 1024;


/* Options for initializing history table memory */
struct HistoryTableMemoryOptions {
    size_t memory_bytes = 0;            // Size of the history table, or 0 to automatically set based on available memory
    double restrict_percent = 0.9;      // Maximum percent of available memory for history table to use
    unsigned thread_count = 1;          // Thread count is used when touching each page on initialization
};


/**
 * Represents one large block of zero-initialized memory used to store the history table
 * 
 * Implementation:
 * - Uses 2 MiB transparent huge pages if available to reduce TLB misses
 * - Interleaves pages across NUMA nodes with RAM if possible to evenly spread memory bandwidth
 * - Each page is touched on initialization so pages are obtained and zeroed at startup
 */
class HistoryTableMemory {
private:
    size_t memory_bytes{}; // Size of the history table
    void *memory_start{}; // Start of usable history table memory
    size_t reserved_bytes{}; // Size of reserved memory
    void *reserved_start{}; // Start of reserved memory

public:
    /**
     * Create history table memory region from config options
     * @param options HistoryTableMemoryOptions struct specifying settings
     */
    explicit HistoryTableMemory(const HistoryTableMemoryOptions &options = {}) {
        /* Calculate size of history table */
        size_t memory_limit = getSystemMemory().availableBytes * options.restrict_percent;
        memory_bytes = options.memory_bytes;
        if (memory_bytes == 0 || memory_bytes > memory_limit)
            memory_bytes = memory_limit;

        /* Round down to nearest complete huge page*/
        memory_bytes = (memory_bytes / huge_page_bytes) * huge_page_bytes;
        if (memory_bytes < huge_page_bytes) memory_bytes = huge_page_bytes;
        
        /* Reserve memory (extra page in case it isn't aligned) */
        reserved_bytes = memory_bytes + huge_page_bytes;
        void* mapping = mmap(nullptr, reserved_bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (mapping == MAP_FAILED) {
            std::cerr << "[HistoryTableMemory] !! Could not map " << (memory_bytes >> 20) << " MB of memory: " << std::strerror(errno) << std::endl;
            exit(1);
        }
        reserved_start = mapping;

        /* Find start of usable memory (aligned to page size) */
        const uintptr_t aligned_address = (reinterpret_cast<uintptr_t>(mapping) + huge_page_bytes - 1) & ~(huge_page_bytes - 1);
        memory_start = reinterpret_cast<void*>(aligned_address);

        /* If huge pages are available, request transparent huge pages */
        #ifdef MADV_HUGEPAGE
        if (madvise(memory_start, memory_bytes, MADV_HUGEPAGE) != 0) {
            std::cout << "[HistoryTableMemory] !! Not using huge pages: " << std::strerror(errno) << std::endl;
        }
        #endif

        /* Interleave memory accross NUMA nodes, if available */
        std::vector<int> node_ids = list_numa_nodes_with_memory();
        if (node_ids.size() > 1) {
            constexpr size_t bits_per_word = sizeof(unsigned long) * 8;
            const int highest_node_id = *std::max_element(node_ids.begin(), node_ids.end());
            const size_t word_count = (static_cast<size_t>(highest_node_id) + bits_per_word + 1) / bits_per_word;
            
            std::vector<unsigned long> node_mask(word_count, 0);
            for (const int node_id : node_ids) {
                node_mask[static_cast<size_t>(node_id) / bits_per_word] |=
                    1ul << (static_cast<size_t>(node_id) % bits_per_word);
            }

            const long result = syscall(
                SYS_mbind, reinterpret_cast<long>(memory_start), static_cast<long>(memory_bytes),
                3, reinterpret_cast<long>(node_mask.data()), static_cast<long>(word_count * bits_per_word), 0L
            );

            if (result == 0) {
                std::cout << "[HistoryTableMemory] Interleaved across NUMA nodes";
                for (const int node_id : node_ids) std::cout << " " << node_id;
                std::cout << std::endl;
            } else {
                std::cout << "[HistoryTableMemory] !! Could not interleave across NUMA nodes: " << strerror(errno) << std::endl;
            }

        } else {
            std::cout << "[HistoryTableMemory] Using one NUMA node or NUMA nodes not found" << std::endl;
        }

        /* Touch every page */
        unsigned char* const block_start = static_cast<unsigned char*>(memory_start);
        const size_t small_page_bytes = static_cast<size_t>(sysconf(_SC_PAGESIZE));
        const size_t huge_pages_per_thread = (memory_bytes / huge_page_bytes) / options.thread_count;
        const size_t bytes_per_thread = huge_pages_per_thread * huge_page_bytes;

        std::vector<std::thread> threads;
        for (unsigned i = 0; i < options.thread_count; i++) {
            const size_t first_byte = i * bytes_per_thread;
            const size_t end_byte = (i == options.thread_count - 1) ? memory_bytes : first_byte + bytes_per_thread;
            threads.emplace_back([=] {
                volatile unsigned char* bytes = block_start;
                for (size_t offset = first_byte; offset < end_byte; offset += small_page_bytes) {
                    bytes[offset] = 0; // Touch one byte in each small page
                }
            });
        }
        for (std::thread& thread : threads) thread.join();
    }

    /**
     * Get start of history table memory
     */
    inline void *get_memory_start() { return memory_start; }

    /**
     * Get size of history table memory
     */
    inline size_t get_memory_bytes() { return memory_bytes; }

    /**
     * Release allocated memory when done using
     */
    void free() {
        if (reserved_start != nullptr) munmap(reserved_start, reserved_bytes);
        reserved_start = nullptr;
        reserved_bytes = 0;
        memory_start = nullptr;
        memory_bytes = 0;
    }

    /**
     * Transfer memory region from another HistoryTableMemory to this one
     * @param other The HistoryTableMemory to transfer from
     */
    void transfer(HistoryTableMemory &other) {
        free();
        reserved_start = std::exchange(other.reserved_start, nullptr);
        reserved_bytes = std::exchange(other.reserved_bytes, 0);
        memory_start = std::exchange(other.memory_start, nullptr);
        memory_bytes = std::exchange(other.memory_bytes, 0);
    }

    /**
     * Destroy history table memory by freeing
     */
    ~HistoryTableMemory() { free(); }

    /**
     * Create history table memory by transfering memory region from another
     * @param other Another HistoryTableMemory to transfer from
     */
    HistoryTableMemory(HistoryTableMemory&& other) noexcept { transfer(other); }

    /**
     * Replace this history table memory with another
     * @param other Another HistoryTableMemory to transfer from
     */
    HistoryTableMemory& operator=(HistoryTableMemory&& other) noexcept {
        if (this != &other)
            transfer(other);
        return *this;
    }

    HistoryTableMemory(const HistoryTableMemory&) = delete;
    HistoryTableMemory& operator=(const HistoryTableMemory&) = delete;

private:
    /**
     * Helper: mask of NUMA nodes that have memory
     */
    static inline std::vector<int> list_numa_nodes_with_memory() {
        std::vector<int> node_ids;
        DIR* node_directory = opendir("/sys/devices/system/node");
        if (node_directory == nullptr) return node_ids;

        // Find all node<number> in directory
        while (const dirent* entry = readdir(node_directory)) {
            int node_id = 0;
            char unexpected_trailing_character = 0;
            if (std::sscanf(entry->d_name, "node%d%c", &node_id, &unexpected_trailing_character) != 1)
                continue;

            char meminfo_path[128];
            std::snprintf(meminfo_path, sizeof meminfo_path, "/sys/devices/system/node/node%d/meminfo", node_id);
            std::FILE* meminfo_file = std::fopen(meminfo_path, "r");
            if (meminfo_file == nullptr)
                continue;

            char line[256];
            while (std::fgets(line, sizeof line, meminfo_file) != nullptr) {
                unsigned long long total_kilobytes = 0;
                if (std::sscanf(line, "Node %*d MemTotal: %llu kB", &total_kilobytes) == 1) {
                    if (total_kilobytes > 0) node_ids.push_back(node_id);
                    break;
                }
            }
            std::fclose(meminfo_file);
        }
        closedir(node_directory);

        std::sort(node_ids.begin(), node_ids.end());

        return node_ids;
    }
};

#endif