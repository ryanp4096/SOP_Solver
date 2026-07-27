#include "trace.hpp"

Trace::Trace() {
    opened = false;
}

void Trace::open(string path) {
    this->path = path;
    file = ofstream(path, std::ios::binary);
    if (file.is_open()) {
        opened = true;
    } else {
        cout << "[Trace] Could not open file " << path << endl;
    }
}

void Trace::close() {
    if (opened && file.is_open()) {
        file.close();
    }
}

void Trace::write(unsigned long long data, size_t bytes)
{
    if (!opened || !file.is_open()) return;
    
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