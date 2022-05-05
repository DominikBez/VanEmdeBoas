// Copyright Dominik Bez 2022.

#pragma once

#include <cstdint>
#include <cassert>
#include <utility>

namespace veb {

/// Specialized hash table for VanEmdeBoas32 to store 8 bit keys with 64 bit values.
/// For uniform random keys, faster locate but slower insert and remove compared to bytell_hash_map.
/// For random clusters, always faster than bytell_hash_map.
class HashTable {
    using uint8_t = std::uint8_t; // using std::uint8_t inside a class is not valid
    using uint32_t = std::uint32_t;
    using uint64_t = std::uint64_t;
    [[maybe_unused]] static constexpr uint32_t lower8BitsMask = (1U << 8) - 1;
    static constexpr uint32_t minTableSize = 4;
    static constexpr uint32_t maxTableSize = 256;
    static constexpr uint8_t hash[256] = {241,  83, 109, 217, 128,  64,  76, 104,  93,  81, 163, 114, 158,
                                          191, 123,  74, 198,  44, 213, 204,  72,  73, 187,  54, 169, 227,
                                           49,  63,  79,  69, 116, 237, 175, 219,  22, 188,  31,  46, 239,
                                          185, 112, 230, 202, 247, 106, 132, 110, 208,  70, 152, 203,  21,
                                          201,  26,  37, 108,  91, 195,  10, 119, 117, 172,  59, 124, 102,
                                           94, 249,  77,  42, 224, 210, 151, 178,  78,  89, 171, 186,  57,
                                           27,  62,  84, 170, 232, 238,  86, 115,  52,   4, 101, 147,  33,
                                          122, 255,   6, 211,  85, 180,  24,  71, 144, 137,  92,  17, 157,
                                          248,  68, 161, 164,  29, 136, 145, 215, 196,  36,  45, 245, 120,
                                           19, 105,  80,  12, 244, 192, 181,  66, 235, 141, 200, 216,  30,
                                           88, 103, 107,  58,  16, 140, 126, 131,  47, 148, 179,   9,   2,
                                          153, 252, 125, 221, 165, 159,  55, 166,  14, 197, 149, 135, 173,
                                          228,   5,  90,  82, 121,  98, 190, 134,   3, 189,  38,  95, 176,
                                          146,  67, 150, 142,  48,  50,   1,  61, 253, 130, 139, 127,  18,
                                          155, 205, 223,  75, 183, 194, 206, 214, 212,   7, 111, 229,  32,
                                          240,  99,  20,  34, 118,  23, 236, 184,   8, 138, 233, 162, 193,
                                          177, 231,  87, 129,  65, 225, 222,  96,  35,  41,  39,  97, 251,
                                          167, 209,  43, 100,  13, 156, 246, 113, 133,  11, 226, 143,  60,
                                          218, 168,  53, 254,  28, 182, 242, 220,  25,  15, 234, 243,   0,
                                          207, 154,  56, 250,  40, 199, 174,  51, 160};

public:
    using Cell = std::pair<uint32_t, uint64_t>; // (key, value); value == 0 is empty cell

    HashTable() : tableSizeMinusOne(minTableSize - 1) {
        table = new Cell[tableSizeMinusOne + 1]();
    }

    ~HashTable() noexcept {
        delete[] table;
    }

    // Assumes a Cell with the given key is contained in this hash table! The key must be 8 bit!
    Cell* find(const uint32_t key) const noexcept {
        assert((key & ~lower8BitsMask) == 0 && "Must only use the lower 8 bits of the key!");
        uint32_t i = hash[key] & tableSizeMinusOne;
        while (true) {
            Cell* cell = table + i;
            if (cell->first == key)
                return cell;
            i = (i + 1) & tableSizeMinusOne;
        }
    }

    // Finds the value to the given key or inserts the key with value 0 if the key is not contained yet.
    // The key must be 8 bit!
    uint64_t& operator[](const uint32_t key) noexcept {
        assert((key & ~lower8BitsMask) == 0 && "Must only use the lower 8 bits of the key!");
        uint32_t i = hash[key] & tableSizeMinusOne;
        while (true) {
            Cell& cell = table[i];
            if (cell.second == 0) { // empty cell
                if ((elements++ << 2) >= 3 * tableSizeMinusOne && tableSizeMinusOne < maxTableSize - 1) { // increase table size
                    uint32_t oldTableSizeMinusOne = tableSizeMinusOne;
                    tableSizeMinusOne = (tableSizeMinusOne << 1) + 1U;
                    Cell* oldTable = table;
                    table = new Cell[tableSizeMinusOne + 1]();
                    for (uint32_t j = 0; j <= oldTableSizeMinusOne; ++j) {
                        const Cell& oldCell = oldTable[j];
                        if (oldCell.second != 0)
                            insert(oldCell);
                    }
                    delete[] oldTable;
                    return insert(Cell(key, 0));
                } else {
                    cell.first = key;
                    return cell.second;
                }
            } else if (cell.first == key) {
                return cell.second;
            }
            i = (i + 1) & tableSizeMinusOne;
        }
    }

    // Assumes there are at least two elements in this hash table and cell->second == 0!
    void remove(Cell* cell) noexcept {
        assert(elements >= 2 && "There must be at least two elements in the table for remove!");
        assert(cell->second == 0);
        if (--elements <= (tableSizeMinusOne >> 2)) {
            uint32_t oldTableSizeMinusOne = tableSizeMinusOne;
            tableSizeMinusOne >>= 1;
            Cell* oldTable = table;
            table = new Cell[tableSizeMinusOne + 1]();
            for (uint32_t j = 0; j <= oldTableSizeMinusOne; ++j) {
                const Cell& oldCell = oldTable[j];
                if (oldCell.second != 0)
                    insert(oldCell);
            }
            delete[] oldTable;
        } else {
            uint32_t pos = cell - table; // position of empty cell
            uint32_t i = pos + 1;
            while (i <= tableSizeMinusOne) {
                Cell* current = table + i;
                if (current->second == 0) {
                    table[pos].second = 0;
                    return; // every cell is at a suitable spot
                }
                const uint32_t currentHash = hash[current->first] & tableSizeMinusOne;
                if (currentHash > i || (currentHash < i && currentHash <= pos)) {
                    table[pos] = *current; // position of current can be improved
                    pos = i;
                }
                ++i;
            }
            i = 0; // wrap around at the end, condition changes since pos > i now
            while (true) {
                Cell* current = table + i;
                if (current->second == 0) {
                    table[pos].second = 0;
                    return; // every cell is at a suitable spot
                }
                const uint32_t currentHash = hash[current->first] & tableSizeMinusOne;
                if (currentHash > i && currentHash <= pos) {
                    table[pos] = *current; // position of current can be improved
                    pos = i;
                    break; // condition changes since pos < i in the next iteration
                }
                ++i;
            }
            ++i;
            while (true) {
                Cell* current = table + i;
                if (current->second == 0) {
                    table[pos].second = 0;
                    return; // every cell is at a suitable spot
                }
                const uint32_t currentHash = hash[current->first] & tableSizeMinusOne;
                if (currentHash > i || (currentHash < i && currentHash <= pos)) {
                    table[pos] = *current; // position of current can be improved
                    pos = i;
                }
                ++i;
            }
        }
    }

private:
    // Assumes there is an empty cell for newCell!
    uint64_t& insert(const Cell& newCell) {
        uint32_t i = hash[newCell.first] & tableSizeMinusOne;
        while (true) {
            Cell& cell = table[i];
            if (cell.second == 0) {
                cell = newCell;
                return cell.second;
            }
            i = (i + 1) & tableSizeMinusOne;
        }
    }

    uint32_t tableSizeMinusOne;
    uint32_t elements = 0;
    Cell* table; // no std::vector such that Bottom fits in a cache line
};

} // namespace veb
