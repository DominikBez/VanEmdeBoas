/* Copyright Dominik Bez 2022.
 *
 * Implementation of a Van Emde Boas set for 32 bit keys. It is based on the paper
 * Roman Dementiev, Lutz Kettner, Jens Mehnert, and Peter Sanders. “Engineering a Sorted List
 * Data Structure for 32 Bit Key”. In: Proceedings of the Sixth Workshop on Algorithm Engineering
 * and Experiments and the First Workshop on Analytic Algorithmics and Combinatorics (ALENEX-04).
 * Ed. by Lars Arge, Giuseppe F. Italiano, and Robert Sedgewick. New Orleans, LA, USA: SIAM, 2004,
 * pp. 142–151. isbn: 0-89871-564-4
 * 
 * Notable differences to the original paper:
 * - This implements only a set and not a map as in the paper
 * - Exploitation of 64 Bit processors (it should work on 32 bit systems too, but is slower):
 *   - The top data structure works on 64 bit values. Therefore, it manages the upper log2(64^3)=18
 *     bits instead of the upper 16 bits. This increases the constant space overhead for an empty
 *     tree from about 0.5 MiB to 2 MiB.
 *   - The Bottom data structure only needs to manage the lower 14 bits. The upper 8 bits of those
 *     are used for its top data structure. Other than in the paper, this Bottom-Top data structure
 *     only consists of one level (4 uint64_t's) since the overhead of another level is not worth
 *     it on 64 bit processors. The remaining 6 bits are used as an bit index into a 64 bit unsigned
 *     integer stored in the hash table. Therefore, this implementation consists of only two levels,
 *     not three as in the original paper.
 * - If a Bottom data structure only contains a single value, this value is stored directly in the
 *   pointer (shifted one bit to the left) and no Bottom data structure must be created. The least
 *   significant bit of the pointer is set to 1 in this case. Note that by alignment of Bottom
 *   this bit is always 0 for a valid pointer. In the original paper, it was not possible to store
 *   the value itself in the pointer since there is associated data which may be too large.
 * - The hash table uses the lower bits (and therefore an "and" operation) instead of the upper bits
 *   (and therefore a shift operation) of the lookup hash value.
 * 
 * The total amount of used memory of this data structure in the worst case (evenly distributed
 * values) is at most 2 MiB + (number of values) * 96 Byte
 * The 2MiB are for the top data structure. Each Bottom data structure needs 64 Bytes (usually one
 * cache line) and contains at least two values (as a single value is stored in the pointer)
 * => at most 32 Bytes per element. The hash table of the Bottom data structure always has a
 * fill rate of at least 1/4. Each Cell in the table consumes 16 Bytes => at most 64 Bytes per
 * element. Overall, the data structure never consumes more than 1 GB since there are at most
 * 2^18 Bottom data structures with at most (64 + 256*16) Bytes each. However, this memory
 * might be scattered all over the heap leading to fragmentation. For equidistant values, the
 * 1 GB limit is reached at about 2^32 / 64 * 3/8 = 25 million values where the last factor
 * comes from the load factor 3/4 of the hash table that is doubled on resize. 25 million
 * values is approximately also the point where std::set consumes more memory for uniform
 * random values than VanEmdeBoas32. Note that the memory consumption of std::set grows way
 * beyond 1 GB for many keys. For non-uniform-random values, the space consumption of VanEmdeBoas32
 * is much smaller than in the worst case. In the best case (many succeeding values) it is as low
 * as 0.25 Bytes per element.
 * 
 * You may define VEB_32_USE_BYTELL to use ska::bytell_hash_map by Malte Skarupke instead of the
 * specialized hash table. It might have a lower space overhead and sometimes have faster
 * insert and delete for uniform random keys. The specialized hash table seems faster for locate.
 */

#pragma once

#include <cstdint>
#include <cstddef>
#include <cassert>
#include <type_traits>

#include "van_emde_boas32_bottom.hpp"
#include "utils/utils.hpp"

namespace veb {

/// An optimized van Emde Boas tree for 32 bit integers or floats in IEEE format.
/// Implements a set and is optimized for 64 bit processors.
/// If T is integral, the value std::numeric_limits<T>::max() is reserved and cannot be inserted or located.
/// If T is float, NaN is reserved and cannot be inserted or located.
/// @tparam T The type to use. Must be uint32_t, int32_t or float in IEEE format.
template<class T = std::uint32_t>
class VanEmdeBoas32 {
    using uint32_t = std::uint32_t; // using std::uint32_t inside a class is not valid
    using uint64_t = std::uint64_t;
    using size_t = std::size_t;
    static_assert(std::is_same_v<T, uint32_t> || std::is_same_v<T, int32_t> || std::is_same_v<T, float>,
        "T must be uint32_t, int32_t or float!");
    static_assert(sizeof(T) == 4 && (!std::is_same_v<T, float> || std::numeric_limits<T>::is_iec559),
        "Floating point type must be in 32 Bit IEE format!");
    
    // If true, all 32 bits fit in tagged pointer, otherwise only bottomK bits. True on 64 bit systems.
    static constexpr bool fullValueInPtr = sizeof(Bottom*) >= 5;
    static constexpr bool signedType = std::is_signed<T>::value;
    static constexpr bool floatType = std::is_floating_point<T>::value;
    static constexpr auto MAX_UT = std::numeric_limits<uint32_t>::max();
    static constexpr auto MIN_UT = std::numeric_limits<uint32_t>::min();
    static constexpr auto MAX_T = std::numeric_limits<T>::max();
    static constexpr auto INVALID = floatType ? std::numeric_limits<T>::quiet_NaN() : MAX_T;
    static constexpr uint32_t signBit = 1U << 31;
    static constexpr int topK = 18;
    static constexpr uint32_t maxUpperBitsValue = (1U << topK) - 1;
    static constexpr int bottomK = 14;
    static constexpr uint32_t bottomMask = (1U << bottomK) - 1;
    static constexpr size_t topSizeLevel1 = 1LLU << topK; // in bits
    static constexpr size_t topSizeLevel2 = topSizeLevel1 / 64;

public:
    static constexpr int BitsK = 32;
    using TypeT = T;

    VanEmdeBoas32() : topLevel1(topSizeLevel1 / sizeof(uint64_t), 0), bottoms(topSizeLevel1, nullptr) { }

    ~VanEmdeBoas32() noexcept {
        uint32_t value = 0;
        while (!empty()) {
            value = locateTop(value);
            Bottom* bottom = bottoms[value];
            if (((uintptr_t)bottom & 1U) == 0)
                delete bottom;
            removeTop(value);
        }
    }

    /// If T is integral: returns std::numeric_limits<T>::max() if value > max or the tree is empty.
    /// If T is floating point: returns NaN if value > max or the tree is empty.
    /// @param sValue The value to search for. If T is integral, sValue must not be std::numeric_limits<T>::max()
    ///              If T is floating point, sValue must not be NaN.
    /// @return The smallest contained number not smaller (i.e., greater or equal) than sValue.
    T locate(const T sValue) const noexcept {
        assert(((floatType && !std::isnan(sValue)) || (!floatType && sValue != MAX_T)) && "sValue is invalid!");
        const uint32_t value = convertToUnsigned(sValue);
        if (value > max || empty()) return INVALID;

        const uint32_t upperBits = value >> bottomK;
        const Bottom* bottom = bottoms[upperBits];
        if (bottom == nullptr) {
            // do nothing here
        } else if (((uintptr_t)bottom & 1U) != 0) {
            uint32_t bottomValue = (uint32_t)((uintptr_t)bottom >> 1);
            if constexpr (!fullValueInPtr) {
                bottomValue = (value & ~bottomMask) | bottomValue;
            }
            if (value <= bottomValue)
                return convertToSigned(bottomValue);
            // do nothing here
        } else if ((value & bottomMask) > bottom->max) {
            // do nothing here
        } else {
            // locate value in bottom
            return convertToSigned((value & ~bottomMask) | bottom->locate(value & bottomMask));
        }

        // locate top, return min of found bottom
        const uint32_t upperResult = locateTop(upperBits + 1);
        bottom = bottoms[upperResult];
        uint32_t lowerResult;
        if (((uintptr_t)bottom & 1U) != 0) {
            lowerResult = (uint32_t)((uintptr_t)bottom >> 1);
            if constexpr (fullValueInPtr) {
                return convertToSigned(lowerResult);
            }
        } else {
            lowerResult = bottom->min;
        }
        return convertToSigned((upperResult << bottomK) | lowerResult);
    }
    /// Alias to locate
    T lower_bound(const T value) const noexcept {
        return locate(value);
    }

    /// Inserts the given value. Does nothing if the value is already contained.
    void insert(const T sValue) noexcept {
        assert(((floatType && !std::isnan(sValue)) || (!floatType && sValue != MAX_T)) && "sValue is invalid!");
        const uint32_t value = convertToUnsigned(sValue);
        max = std::max(max, value);
        const uint32_t upperBits = value >> bottomK;
        Bottom*& bottom = bottoms[upperBits];
        if (bottom == nullptr) {
            if constexpr (fullValueInPtr)
                bottom = (Bottom*)(((uintptr_t)value << 1) | 1U);
            else
                bottom = (Bottom*)((((uintptr_t)value & bottomMask) << 1) | 1U);
            insertTop(upperBits);
        } else if (((uintptr_t)bottom & 1U) != 0) {
            uint32_t oldValue = (uint32_t)((uintptr_t)bottom >> 1);
            if constexpr (fullValueInPtr) {
                if (oldValue == value)
                    return; // value already contained
                oldValue = oldValue & bottomMask; // only insert bottomK bits
            } else if (oldValue == (value & bottomMask)) {
                return; // value already contained
            }
            bottom = new Bottom;
            bottom->insert(oldValue);
            bottom->insert(value & bottomMask);
        } else {
            bottom->insert(value & bottomMask);
        }
    }

    /// Removes the given value. The value must be contained in this tree!
    void remove(const T sValue) noexcept {
        assert(((floatType && !std::isnan(sValue)) || (!floatType && sValue != MAX_T)) && "sValue is invalid!");
        const uint32_t value = convertToUnsigned(sValue);
        const uint32_t upperBits = value >> bottomK;
        Bottom*& bottom = bottoms[upperBits];
        if (((uintptr_t)bottom & 1U) != 0) {
            bottom = nullptr;
            removeTop(upperBits);
        } else {
            bottom->remove(value & bottomMask);
            if (bottom->min == bottom->max) {
                uint32_t oldValue = bottom->min;
                if constexpr (fullValueInPtr) {
                    oldValue |= value & ~bottomMask;
                }
                delete bottom;
                bottom = (Bottom*)(((uintptr_t)oldValue << 1) | 1U);
            }
        }

        if (value == max) {
            if (empty()) {
                max = MIN_UT;
            } else {
                const uint32_t upperBits = maxTop();
                Bottom* bottom = bottoms[upperBits];
                if (((uintptr_t)bottom & 1U) != 0) {
                    const uint32_t lowerBits = (uint32_t)((uintptr_t)bottom >> 1);
                    if constexpr (fullValueInPtr) {
                        max = lowerBits;
                    } else {
                        max = (upperBits << bottomK) | lowerBits;
                    }
                } else {
                    max = (upperBits << bottomK) | bottom->max;
                }
            }
        }
    }
    /// Alias to remove
    void erase(const T value) noexcept {
        remove(value);
    }

    bool empty() const noexcept {
        return topLevel3 == 0;
    }

private:
    static constexpr uint32_t convertToUnsigned(const T value) noexcept {
        if constexpr (floatType) {
            return reinterpret_cast<const uint32_t&&>(std::move(value)) ^ signBit; // just flip the sign
        } else if constexpr (signedType) {
            return reinterpret_cast<const uint32_t&&>(value) ^ signBit; // just flip the sign
        } else {
            return value;
        }
    }
    static constexpr T convertToSigned(const uint32_t value) noexcept {
        if constexpr (signedType) {
            return reinterpret_cast<const T&&>(value ^ signBit); // just flip the sign
        } else {
            return value;
        }
    }

    uint32_t locateTop(const uint32_t value) const noexcept {
        constexpr int bits = 6;
        constexpr uint32_t mask = (1U << bits) - 1;
        uint32_t valueShifted = value >> bits;
        uint64_t current = topLevel1[valueShifted] & ~((1ULL << (value & mask)) - 1);
        if (current != 0)
            return (value & ~mask) | utils::countr_zero(current);
        ++valueShifted;
        const uint32_t valueShiftedTwice = valueShifted >> bits;
        current = topLevel2[valueShiftedTwice] & ~((1ULL << (valueShifted & mask)) - 1);
        if (current != 0) {
            const uint32_t upperBits = (valueShifted & ~mask) | utils::countr_zero(current);
            return (upperBits << bits) | utils::countr_zero(topLevel1[upperBits]);
        }
        current = topLevel3 & ~((1ULL << (valueShiftedTwice + 1)) - 1);
        assert(current != 0);
        uint32_t upperBits = utils::countr_zero(current);
        upperBits = (upperBits << bits) | utils::countr_zero(topLevel2[upperBits]);
        return (upperBits << bits) | utils::countr_zero(topLevel1[upperBits]);
    }

    void insertTop(const uint32_t value) noexcept {
        constexpr int bits = 6;
        constexpr uint32_t mask = (1U << bits) - 1;
        const uint32_t valueShifted = value >> bits;
        uint64_t& level1 = topLevel1[valueShifted];
        if (level1 == 0) {
            const uint32_t valueShiftedTwice = value >> (2 * bits);
            uint64_t& level2 = topLevel2[valueShiftedTwice];
            if (level2 == 0) {
                topLevel3 |= (1LLU << valueShiftedTwice);
            }
            level2 |= (1LLU << (valueShifted & mask));
        }
        level1 |= (1LLU << (value & mask));
    }

    // Removes the value from the top data structure. Only the lower 18 bits of value are
    // considered and the value is assumed to be contained in the top data structure.
    void removeTop(const uint32_t value) noexcept {
        constexpr int bits = 6;
        constexpr uint32_t mask = (1U << bits) - 1;
        const uint32_t valueShifted = value >> bits;
        uint64_t& level1 = topLevel1[valueShifted];
        level1 ^= (1LLU << (value & mask));
        if (level1 != 0) return;
        const uint32_t valueShiftedTwice = value >> (2 * bits);
        uint64_t& level2 = topLevel2[valueShiftedTwice];
        level2 ^= (1LLU << (valueShifted & mask));
        if (level2 != 0) return;
        topLevel3 ^= (1LLU << valueShiftedTwice);
    }

    /// Returns the maximum value in the top data structure. Assumes it is not empty!
    uint32_t maxTop() const noexcept {
        constexpr int bits = 6;
        uint32_t result = 63 - utils::countl_zero(topLevel3);
        result = (result << bits) | (63 - utils::countl_zero(topLevel2[result]));
        return (result << bits) | (63 - utils::countl_zero(topLevel1[result]));
    }

    uint32_t max = MIN_UT;
    std::vector<uint64_t> topLevel1;
    uint64_t topLevel2[topSizeLevel2 / sizeof(uint64_t)] = {0}; // use a vector if sizeof(VanEmdeBoas32) is too large
    uint64_t topLevel3 = 0;
    std::vector<Bottom*> bottoms;
};

} // namespace veb
