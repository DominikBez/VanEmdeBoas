// Copyright Dominik Bez 2022.
//
// You may define VEB_32_USE_BYTELL to use ska::bytell_hash_map by Malte Skarupke instead of the
// specialized hash table.

#pragma once

#include <cassert>

#include "utils/utils.hpp"

#ifdef VEB_32_USE_BYTELL
#include "utils/bytell_hash_map.hpp"
#else
#include "hash_table.hpp"
#endif

namespace veb {

class alignas(64) Bottom { // sizeof(Bottom) <= 64 => fits in cache line
protected:
    using uint32_t = std::uint32_t; // using std::uint32_t inside a class is not valid
    using uint64_t = std::uint64_t;

    static constexpr int lowerK = 6;
    static constexpr uint32_t lowerMask = (1U << lowerK) - 1;
    static constexpr int upperLowerK = 6;
    static constexpr uint32_t upperLowerMask = (1U << upperLowerK) - 1;

public:
    uint32_t locate(const uint32_t value) const noexcept {
        const uint32_t upperBits = value >> lowerK;
        const uint32_t upperLowerBits = upperBits & upperLowerMask;
        const uint64_t upperLowerBitsMask = 1LLU << upperLowerBits;
        uint64_t mask = ~(upperLowerBitsMask - 1);
        uint32_t switchCase = upperBits >> upperLowerK;
        uint32_t nextTop = 0;
        uint32_t upperResult;
        int retry = 0;
        do { // alternative to goto, at most two iterations
            uint64_t maskedOut;
            switch (switchCase) {
                case 0:
                    maskedOut = top[0] & mask;
                    if (maskedOut == 0)
                        nextTop = 1;
                    else
                        upperResult = utils::countr_zero(maskedOut);
                    break;
                case 1:
                    maskedOut = top[1] & mask;
                    if (maskedOut == 0)
                        nextTop = 2;
                    else
                        upperResult = utils::countr_zero(maskedOut) + 64;
                    break;
                case 2:
                    maskedOut = top[2] & mask;
                    if (maskedOut == 0)
                        nextTop = 3;
                    else
                        upperResult = utils::countr_zero(maskedOut) + 2 * 64;
                    break;
                case 3:
                    maskedOut = top[3] & mask;
                    assert((maskedOut != 0) && "Last top in bottom is zero masked out!");
                    upperResult = utils::countr_zero(maskedOut) + 3 * 64;
                    break;
                default:
                    assert(false && "Wrong top case in Bottom::locate!");
                    return (uint32_t)-1; // suppress uninitialized variable warnings
            }
            uint64_t topValue;
            switch (nextTop) {
                case 0:
                    if (retry == 0 && (maskedOut & upperLowerBitsMask) != 0) {
                        const uint64_t bits = innerBottoms.find(upperResult)->second
                            & ~((1LLU << (value & lowerMask)) - 1);
                        if (bits != 0)
                            return (upperResult << lowerK) | utils::countr_zero(bits);
                        retry = 1; // take minimum of next non-zero innerBottom instead
                        if (upperLowerBits == upperLowerMask) {
                            ++switchCase; // we have to look at the next uint64_t in top
                            mask = ~(uint64_t)0;
                        } else {
                            mask ^= upperLowerBitsMask; // we do not want to reach this case again
                        }
                    }
                    break; // upperResult already found
                case 1:
                    topValue = top[1];
                    if (topValue != 0) {
                        upperResult = utils::countr_zero(topValue) + 64;
                        break;
                    } // else fallthrough
                    [[fallthrough]];
                case 2:
                    topValue = top[2];
                    if (topValue != 0) {
                        upperResult = utils::countr_zero(topValue) + 2 * 64;
                        break;
                    } // else fallthrough
                    [[fallthrough]];
                case 3:
                    topValue = top[3];
                    assert((topValue != 0) && "Last top in bottom is zero!");
                    upperResult = utils::countr_zero(topValue) + 3 * 64;
                    break;
                default:
                    assert(false && "Wrong nextTop in bottom!");
                    break;
            }
        } while (++retry == 2); // at most once true

        return (upperResult << lowerK) | utils::countr_zero(innerBottoms.find(upperResult)->second);
    }

    void insert(const uint32_t value) noexcept {
        min = std::min(min, value);
        max = std::max(max, value);
        const uint32_t upperBits = value >> lowerK;
        top[upperBits >> upperLowerK] |= (1LLU << (upperBits & upperLowerMask));
        innerBottoms[upperBits] |= (1LLU << (value & lowerMask));
    }

    /// Removes the given value. The value and at least one other value must be contained!
    void remove(const uint32_t value) noexcept {
        assert(min != max);
        const uint32_t upperBits = value >> lowerK;
        auto bottomCell = innerBottoms.find(upperBits);
        uint64_t& bottom = bottomCell->second;
        bottom ^= (1LLU << (value & lowerMask));
        if (bottom == 0) {
            #ifdef VEB_32_USE_BYTELL
            innerBottoms.erase(upperBits);
            #else
            innerBottoms.remove(bottomCell);
            #endif
            top[upperBits >> upperLowerK] ^= (1LLU << (upperBits & upperLowerMask));
        }
        if (value == min) {
            uint32_t upperResult;
            uint64_t topValue;
            switch (upperBits >> upperLowerK) {
                case 0:
                    topValue = top[0];
                    if (topValue != 0) {
                        upperResult = utils::countr_zero(topValue);
                        break;
                    } // else fallthrough
                    [[fallthrough]];
                case 1:
                    topValue = top[1];
                    if (topValue != 0) {
                        upperResult = utils::countr_zero(topValue) + 64;
                        break;
                    } // else fallthrough
                    [[fallthrough]];
                case 2:
                    topValue = top[2];
                    if (topValue != 0) {
                        upperResult = utils::countr_zero(topValue) + 2 * 64;
                        break;
                    } // else fallthrough
                    [[fallthrough]];
                case 3:
                    topValue = top[3];
                    upperResult = utils::countr_zero(topValue) + 3 * 64;
                    break;
                default:
                    assert(false && "Wrong top case in Bottom::remove!");
                    return; // suppress uninitialized variable warnings
            }
            min = (upperResult << lowerK) | utils::countr_zero(innerBottoms.find(upperResult)->second);
        } else if (value == max) {
            uint32_t upperResult;
            uint64_t topValue;
            switch (upperBits >> upperLowerK) { // similar as for value == min but backwards
                case 3:
                    topValue = top[3];
                    if (topValue != 0) {
                        upperResult = 63 - utils::countl_zero(topValue) + 3 * 64;
                        break;
                    } // else fallthrough
                    [[fallthrough]];
                case 2:
                    topValue = top[2];
                    if (topValue != 0) {
                        upperResult = 63 - utils::countl_zero(topValue) + 2 * 64;
                        break;
                    } // else fallthrough
                    [[fallthrough]];
                case 1:
                    topValue = top[1];
                    if (topValue != 0) {
                        upperResult = 63 - utils::countl_zero(topValue) + 64;
                        break;
                    } // else fallthrough
                    [[fallthrough]];
                case 0:
                    topValue = top[0];
                    upperResult = 63 - utils::countl_zero(topValue);
                    break;
                default:
                    assert(false && "Wrong top case in Bottom::remove!");
                    return; // suppress uninitialized variable warnings
            }
            max = (upperResult << lowerK) | (63 - utils::countl_zero(innerBottoms.find(upperResult)->second));
        }
    }

    uint32_t min = std::numeric_limits<uint32_t>::max();
    uint32_t max = std::numeric_limits<uint32_t>::min();

protected:
    uint64_t top[4] = {0};
    #ifdef VEB_32_USE_BYTELL
    ska::bytell_hash_map<uint32_t, uint64_t> innerBottoms;
    #else
    HashTable innerBottoms;
    #endif
};

} // namespace veb
