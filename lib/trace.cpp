#include "trace.hpp"

#ifdef ENABLE_TRACE

Trace::Trace() {
    opened = false;
}

void Trace::open(string path, int instance_size, TraceDetailLevel detail_level, int thread_id) {
    this->path = path;
    this->instance_size = instance_size;
    this->detail_level = detail_level;
    this->thread_id = thread_id;
    file = ofstream(path, std::ios::binary);
    if (file.is_open()) {
        opened = true;
        instance_size_shift = 0;
        int nodes_possible = 255;
        while (instance_size > nodes_possible) {
            instance_size_shift++;
            nodes_possible = (255 - instance_size_shift) + 256 * instance_size_shift;
        }
    } else {
        cout << "[Trace] Could not open file " << path << endl;
    }
}

void Trace::write_header() {
    write(TRACE_VERSION_NUMBER, 4);
    write(thread_id, 1);
    write(detail_level, 1);
    write(instance_size_shift, 1);
    write(instance_size, 4);
}

void Trace::close() {
    if (opened && file.is_open()) {
        file.close();
    }
}

bool Trace::is_open() {
    return opened;
}

void Trace::write(unsigned long long data, size_t bytes)
{
    if (!opened) return;
    
    std::uint8_t d8;
    std::uint16_t d16;
    std::uint32_t d32;
    std::uint64_t d64;
    char *ptr;

    switch (bytes) {
        case 1:
            if (data > UINT8_MAX) data = UINT8_MAX;
            d8 = static_cast<std::uint8_t>(data);
            ptr = reinterpret_cast<char *>(&d8);
            break;

        case 2:
            if (data > UINT16_MAX) data = UINT16_MAX;
            d16 = static_cast<std::uint16_t>(data);
            ptr = reinterpret_cast<char *>(&d16);
            break;

        case 4:
            if (data > UINT32_MAX) data = UINT32_MAX;
            d32 = static_cast<std::uint32_t>(data);
            ptr = reinterpret_cast<char *>(&d32);
            break;

        case 8:
            if (data > UINT64_MAX) data = UINT64_MAX;
            d64 = static_cast<std::uint64_t>(data);
            ptr = reinterpret_cast<char *>(&d64);
            break;

        default:
            return;
    }
    file.write(ptr, bytes);
}

void Trace::write_detail(unsigned long long data, size_t bytes) {
    if (detail_level == DETAIL_NORMAL) {
        write(data, bytes);
    }
}

void Trace::write_node(int node) {
    if (node < 255 - instance_size_shift) {
        write(node, 1);
    } else {
        node -= 255 - instance_size_shift;
        int first_byte = 254 - (node / 256);
        int second_byte = node % 256;
        write(first_byte, 1);
        write(second_byte, 1);
    }
}

void Trace::write_end_list() {
    write(0xff, 1);
}

#else

/*
 * If trace is disabled, trace class has no behavior.
 * This allows compiler to optimize and remove trace functions to improve runtime.
 */

Trace::Trace() {}
void Trace::open(std::string path, int instance_size, TraceDetailLevel detail_level, int thread_id) {}
void write_header() {}
void Trace::close() {}
bool Trace::is_open() { return false; }
void Trace::write(unsigned long long data, size_t bytes) {}
void Trace::write_detail(unsigned long long data, size_t bytes) {}
void Trace::write_node(int node) {}
void Trace::write_end_list() {}

#endif