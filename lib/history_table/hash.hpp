/**
 * Defines the hashing algorithm used to locate a key in the history table
 */

#ifndef HISTORY_HASH_H
#define HISTORY_HASH_H

#include <cstdint>
#include <cstddef>
#include <vector>
#include <random>
#include <boost/dynamic_bitset.hpp>
#include "entries.hpp"

/**
 * Represents a 128-bit hash.
 */
struct Hash128 {
    uint64_t high; // First 64 bits
    uint64_t low; // Second 64 bits

    friend constexpr Hash128 operator^(Hash128 a, Hash128 b) {
        return {a.high ^ b.high, a.low ^ b.low};
    }

    friend constexpr bool operator==(Hash128 a, Hash128 b) {
        return a.high == b.high && a.low == b.low;
    }
};


/**
 * Zobrist-style hashing for keys in the history table.
 * 
 * Features:
 * - Hash values are uniformly distributed
 * - Can modify an existing hash when a bit is changed or first/last
 *   nodes are changed without recomputing the hash
 * - The chance that two different keys will have the same hash is 2^-128,
 *   which is far less likely than regular hardware errors even under
 *   the worst conditions
 * 
 * Implementation:
 *   On initialization, a table of random 128-bit values is generated, one
 *   for each bit in the bitset and for each possible first and last node.
 *   The hash is the XOR of the random values for each present bit
 *   in the bitset and the random values for the first and last nodes.
 * 
 * Seeding:
 *   Set seed to 0 to use OS's non-deterministic random number generator
 *   Other seeds will use SplitMix64 to generate deterministic pseudo-random numbers
 */
class Hasher {
private:
    size_t instance_size; // Number of bits in bitset
    uint64_t seed; // Seed for random numbers in table

    Hash128 empty_hash_value; // Hash for empty history key (empty bitset, first = -1, last = -1)

    std::vector<Hash128> bitset_random_table; // Random value for each bit in bitset
    std::vector<Hash128> first_node_random_table; // Random value for each possible first node
    std::vector<Hash128> last_node_random_table; // Random value for each possible last node

public:
    /**
     * Uninitialized hasher
     */
    Hasher() {}

    /**
     * Create new hasher
     * @param instance_size Size of instance
     * @param seed Set to 0 (default) to use non-deterministic random numbers,
     *             otherwise use seed for deterministic pseudo-random numbers
     */
    Hasher(size_t instance_size, uint64_t seed = 0)
          : instance_size(instance_size),
            seed(seed),
            bitset_random_table(instance_size),
            first_node_random_table(instance_size + 1), // shifted by 1 to support first_node == -1 (empty has no first node)
            last_node_random_table(instance_size + 1) {
        
        std::random_device operating_system_entropy{};
        uint64_t random_state = seed;

        // Function to generate next random number
        auto next_random = [&]() -> uint64_t {
            if (seed == 0) {
                // Seed == 0: use operating system's non-determinstic random number generator
                const uint64_t high_bits = operating_system_entropy();  // 32 random bits
                const uint64_t low_bits = operating_system_entropy();   // 32 random bits
                return (high_bits << 32) | low_bits;

            } else {
                // Otherwise: use SplitMix64 to generate deterministic pseudo-random numbers based on seed
                random_state += 0x9E3779B97F4A7C15ull;
                uint64_t mixed = random_state;
                mixed = (mixed ^ (mixed >> 30)) * 0xBF58476D1CE4E5B9ull;
                mixed = (mixed ^ (mixed >> 27)) * 0x94D049BB133111EBull;
                return mixed ^ (mixed >> 31);
            }
        };

        // Fill all random tables with random numbers generated using function above
        for (Hash128& entry : bitset_random_table) {
            entry.high = next_random();
            entry.low = next_random();
        }
        for (Hash128& entry : first_node_random_table) {
            entry.high = next_random();
            entry.low = next_random();
        }
        for (Hash128& entry : last_node_random_table) {
            entry.high = next_random();
            entry.low = next_random();
        }

        empty_hash_value = first_node_random_table[0] ^ last_node_random_table[0];
    }

    /**
     * Get hash for an empty key without calculating
     * Empty key is an empty bitset, first = -1, last = -1
     * @return Hash for an empty key
     */
    Hash128 empty_hash() const {
        return empty_hash_value;
    }

    /**
     * Get hash for a key containing one bit
     * @param bit_index The index of the only present bit in the key
     * @return Hash for a key containing only bit_index
     */
    Hash128 hash_single_bit(int bit_index) const {
        return bitset_random_table[bit_index] ^ first_node_random_table[bit_index + 1] ^ last_node_random_table[bit_index + 1];
    }

    /**
     * Compute hash of a history key from scratch
     * @param bit_vector The bit vector of the key to hash
     * @param first The first node of the key to hash
     * @param last The last node of the key to hash
     * @return Hash of the key
     */
    Hash128 hash(const boost::dynamic_bitset<> bit_vector, int first, int last) const {
        Hash128 hash = first_node_random_table[first + 1] ^ last_node_random_table[last + 1];

        for (size_t i = 0; i < instance_size; i++) {
            if (bit_vector[i])
                hash = hash ^ bitset_random_table[i];
        }

        return hash;
    }

    /**
     * Compute hash of a history key from scratch
     * @param key The key to hash
     * @return Hash of the key
     */
    Hash128 hash(const HistoryKey& key) const {
        return hash(key.bit_vector, key.first_node, key.last_node);
    }

    /**
     * Modify an existing hash to toggle a bit in the bitset
     * @param hash Old hash
     * @param bit_index Index of bit to toggle in bitset
     * @return New hash after bit toggled
     */
    Hash128 toggle_bit(Hash128 hash, int bit_index) const {
        return hash ^ bitset_random_table[bit_index];
    }

    /**
     * Modify an existing hash to change the first node
     * @param hash Old hash
     * @param old_first Old first node
     * @param new_first New first node
     * @return New hash after first node changed
     */
    Hash128 change_first_node(Hash128 hash, int old_first, int new_first) const {
        return hash ^ first_node_random_table[old_first + 1] ^ first_node_random_table[new_first + 1];
    }

    /**
     * Modify an existing hash to change the last node
     * @param hash Old hash
     * @param old_last Old last node
     * @param new_last New last node
     * @return New hash after last node changed
     */
    Hash128 change_last_node(Hash128 hash, int old_last, int new_last) const {
        return hash ^ last_node_random_table[old_last + 1] ^ last_node_random_table[new_last + 1];
    }
};


/**
 * Automatically updating hash as nodes are added/removed from a history key
 * Assumes that the same node is never added to the path multiple times
 */
class HashedHistoryKey {
private:
    const Hasher &hasher;
    int first;
    int last;
    Hash128 hash;

public:
    /**
     * Creates new empty key (empty bitset, first = -1, last = -1)
     * @param hasher The hasher used to generate hashes
     */
    HashedHistoryKey(const Hasher &hasher)
          : hasher(hasher),
            first(-1),
            last(-1),
            hash(hasher.empty_hash()) {}
    
    /**
     * Creates new hashed key from history key
     * @param hasher The hasher used to generate hashes
     * @param key The key to create a hashed key from
     */
    HashedHistoryKey(const Hasher &hasher, const HistoryKey &key)
          : hasher(hasher),
            first(key.first_node),
            last(key.last_node),
            hash(hasher.hash(key)) {}
    
    /**
     * Creates new hashed key from history key parts (bit vector, first node, and last node)
     * @param hasher The hasher used to generate hashes
     * @param bit_vector The bit vector of the key
     * @param first The first node of the key
     * @param last The last node of the key
     */
    HashedHistoryKey(const Hasher &hasher, const boost::dynamic_bitset<> &bit_vector, int first, int last)
          : hasher(hasher),
            first(first),
            last(last),
            hash(hasher.hash(bit_vector, first, last)) {}
    
    /**
     * Get current hash value. Already stored, does not need to compute.
     */
    inline Hash128 get_hash() const { return hash; }
    
    /**
     * Set to a different key. Requires re-hashing the entire key
     * @param bit_vector The bit vector of the key
     * @param first_node The first node of the key
     * @param last_node The last node of the key
     */
    void set(const boost::dynamic_bitset<> &bit_vector, int first_node, int last_node) {
        first = first_node;
        last = last_node;
        hash = hasher.hash(bit_vector, first, last);
    }

    /**
     * Set to a different key. Requires re-hashing the entire key
     * @param key The key to set this hashed key to
     */
    void set(const HistoryKey &key) {
        set(key.bit_vector, key.first_node, key.last_node);
    }

    /**
     * Push a node to the end of the history key and update hash
     * @param node Node to push
     */
    void push_last(int node) {
        if (first == -1) { // if empty, the new node is both the first and last
            hash = hasher.change_first_node(hash, first, node);
            first = node;
        }
        hash = hasher.change_last_node(hash, last, node);
        last = node;
        hash = hasher.toggle_bit(hash, node);
    }

    /**
     * Push a node to the start of the history key and update hash
     * @param node Node to push
     */
    void push_first(int node) {
        if (last == -1) { // if empty, the new node is both the first and last
            hash = hasher.change_last_node(hash, last, node);
            last = node;
        }
        hash = hasher.change_first_node(hash, first, node);
        first = node;
        hash = hasher.toggle_bit(hash, node);
    }

    /**
     * Pop a node from the end of the history key and update hash
     * @param new_last The new last node after the last node is popped
     */
    void pop_last(int new_last) {
        if (last == -1) return;
        hash = hasher.toggle_bit(hash, last);
        hash = hasher.change_last_node(hash, last, new_last);
        last = new_last;
        if (new_last == -1) { // the key is now empty, so first is also -1
            hash = hasher.change_first_node(hash, first, -1);
            first = -1;
        }
    }

    /**
     * Pop a node from the start of the history key and update hash
     * @param new_first The new first node after the first node is popped
     */
    void pop_first(int new_first) {
        if (first == -1) return;
        hash = hasher.toggle_bit(hash, first);
        hash = hasher.change_first_node(hash, first, new_first);
        first = new_first;
        if (new_first == -1) { // the key is now empty, so last is also -1
            hash = hasher.change_last_node(hash, last, -1);
            last = -1;
        }
    }
};

#endif