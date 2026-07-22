#include "memory.hpp"

#if defined(__linux__)

#include <cstdint>
#include <fstream>
#include <string>

MemoryInfo getSystemMemory() {
    std::ifstream file("/proc/meminfo");

    uint64_t totalBytes = 0;
    uint64_t availableBytes = 0;

    std::string key;
    uint64_t value;
    std::string unit;

    while (file >> key >> value >> unit) {
        if (key == "MemTotal:") {
            totalBytes = value * 1024;
        } else if (key == "MemAvailable:") {
            availableBytes = value * 1024;
        }

        if (totalBytes != 0 && availableBytes != 0) {
            break;
        }
    }

    return {
        .totalBytes = totalBytes,
        .availableBytes = availableBytes
    };
}

#elif defined(__APPLE__)

#include <mach/mach.h>
#include <sys/sysctl.h>

MemoryInfo getSystemMemory() {
    // Get total physical memory
    uint64_t totalBytes = 0;
    size_t size = sizeof(totalBytes);

    sysctlbyname(
        "hw.memsize",
        &totalBytes,
        &size,
        nullptr,
        0
    );

    // Get VM page size
    vm_size_t pageSize = 0;
    host_page_size(mach_host_self(), &pageSize);

    // Get VM statistics
    vm_statistics64_data_t vmStats{};
    mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;

    host_statistics64(
        mach_host_self(),
        HOST_VM_INFO64,
        reinterpret_cast<host_info64_t>(&vmStats),
        &count
    );

    uint64_t freePages = static_cast<uint64_t>(vmStats.free_count);
    uint64_t availablePages = freePages +
        static_cast<uint64_t>(vmStats.inactive_count) +
        static_cast<uint64_t>(vmStats.speculative_count);

    uint64_t freeBytes = freePages * static_cast<uint64_t>(pageSize);
    uint64_t availableBytes = availablePages * static_cast<uint64_t>(pageSize);

    return {
        .totalBytes = totalBytes,
        .freeBytes = freeBytes,
        .availableBytes = availableBytes
    };
}

#else

#error "Unsupported operating system"

#endif