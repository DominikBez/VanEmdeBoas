/* Copyright Dominik Bez 2022.
 *
 * Parallel implementation of a Van Emde Boas set for 32 bit keys. It is based on the paper
 * Roman Dementiev, Lutz Kettner, Jens Mehnert, and Peter Sanders. “Engineering a Sorted List
 * Data Structure for 32 Bit Key”. In: Proceedings of the Sixth Workshop on Algorithm Engineering
 * and Experiments and the First Workshop on Analytic Algorithmics and Combinatorics (ALENEX-04).
 * Ed. by Lars Arge, Giuseppe F. Italiano, and Robert Sedgewick. New Orleans, LA, USA: SIAM, 2004,
 * pp. 142–151. isbn: 0-89871-564-4
 * 
 * For an explanation, see the sequential version in van_emde_boas32.hpp.
 * This implementation here is a relatively naive version where the whole top data structure is
 * locked on every access.
 * 
 * You may define VEB_32_USE_BYTELL to use ska::bytell_hash_map by Malte Skarupke instead of the
 * specialized hash table. Furthermore, you may define VEB_32_LOCKED_TOP_USE_STD_MUTEX to use
 * std::mutex instead of sf::contention_free_shared_mutex<>
 * (https://github.com/AlexeyAB/object_threadsafe). With std::mutex, locate acquires a unique_lock
 * and therefore only a single thread can locate at a time. sf::contention_free_shared_mutex<>
 * allows to use a shared_lock in locate. std::shared_mutex would allow the same but is much
 * slower than std::mutex (probably the writers are starving). sf::contention_free_shared_mutex<>
 * seems faster even with a single reader.
 */

#pragma once

#include <cstdint>
#include <cstddef>
#include <cassert>
#include <type_traits>
#include <mutex>

#include "van_emde_boas32_bottom.hpp"
#include "utils/utils.hpp"

#ifndef VEB_32_LOCKED_TOP_USE_STD_MUTEX
#include <shared_mutex> // for shared_lock; shared_mutex is much slower than mutex
#include "utils/safe_ptr.h"
#endif

namespace veb {

/// A concurrent (sequentially consistent) van Emde Boas tree for 32 bit integers or floats in IEEE format.
/// Implements a set and is optimized for 64 bit processors.
/// If T is integral, the value std::numeric_limits<T>::max() is reserved and cannot be inserted or located.
/// If T is float, NaN is reserved and cannot be inserted or located.
/// @tparam T The type to use. Must be uint32_t, int32_t or float in IEEE format.
template<class T = std::uint32_t>
class VanEmdeBoas32LockedTop {
    using uint8_t = std::uint8_t; // using std::uint8_t inside a class is not valid
    using uint32_t = std::uint32_t;
    using uint64_t = std::uint64_t;
    using size_t = std::size_t;
    static_assert(std::is_same_v<T, uint32_t> || std::is_same_v<T, int32_t> || std::is_same_v<T, float>,
        "T must be uint32_t, int32_t or float!");
    static_assert(sizeof(T) == 4 && (!std::is_same_v<T, float> || std::numeric_limits<T>::is_iec559),
        "Floating point type must be in 32 Bit IEE format!");
    
    // Specialization of Bottom since min/max for insert are updated outside of Bottom.
    // sizeof(BottomMinMaxOutside) <= 64 => fits in cache line
    class alignas(64) BottomMinMaxOutside : public Bottom {
    public:
        void insert(const uint32_t value) noexcept {
            const uint32_t upperBits = value >> lowerK;
            top[upperBits >> upperLowerK] |= (1LLU << (upperBits & upperLowerMask));
            innerBottoms[upperBits] |= (1LLU << (value & lowerMask));
        }
    };

    // If true, all 32 bits fit in tagged pointer, otherwise only bottomK bits. True on 64 bit systems.
    static constexpr bool fullValueInPtr = sizeof(BottomMinMaxOutside*) >= 5;
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

    VanEmdeBoas32LockedTop()
        : topLevel1(topSizeLevel1 / sizeof(uint64_t), 0),
          bottoms(topSizeLevel1, nullptr),
          bottomMutexes(topSizeLevel1)
        { }

    ~VanEmdeBoas32LockedTop() noexcept {
        uint32_t value = 0;
        while (topLevel3 != 0) {
            value = locateTop(value);
            BottomMinMaxOutside* bottom = bottoms[value];
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
        #ifdef VEB_32_LOCKED_TOP_USE_STD_MUTEX
        std::unique_lock topLock(topMutex);
        #else
        std::shared_lock topLock(topMutex);
        #endif
        if (value > max || topLevel3 == 0) return INVALID;

        const uint32_t upperBits = value >> bottomK;
        const BottomMinMaxOutside* bottom = bottoms[upperBits];
        if (bottom == nullptr) {
            // do nothing here
        } else if (((uintptr_t)bottom & 1U) != 0) {
            uint32_t bottomValue = (uint32_t)((uintptr_t)bottom >> 1);
            if constexpr (!fullValueInPtr) {
                bottomValue = (value & ~bottomMask) | bottomValue;
            }
            if (value <= bottomValue) {
                topLock.unlock();
                return convertToSigned(bottomValue);
            }
            // do nothing here
        } else if ((value & bottomMask) > bottom->max) {
            // do nothing here
        } else {
            // locate value in bottom
            std::lock_guard bottomLock(bottomMutexes[upperBits]);
            topLock.unlock(); // we don't need it anymore since bottom is found and locked
            return convertToSigned((value & ~bottomMask) | bottom->locate(value & bottomMask));
        }

        // locate top, return min of found bottom
        const uint32_t upperResult = locateTop(upperBits + 1);
        bottom = bottoms[upperResult];
        uint32_t lowerResult;
        if (((uintptr_t)bottom & 1U) != 0) {
            lowerResult = (uint32_t)((uintptr_t)bottom >> 1);
            if constexpr (fullValueInPtr) {
                topLock.unlock();
                return convertToSigned(lowerResult);
            }
        } else {
            lowerResult = bottom->min;
        }
        topLock.unlock();
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
        const uint32_t upperBits = value >> bottomK;
        BottomMinMaxOutside*& bottom = bottoms[upperBits];
        std::unique_lock topLock(topMutex);
        max = std::max(max, value);
        if (bottom == nullptr) {
            if constexpr (fullValueInPtr)
                bottom = (BottomMinMaxOutside*)(((uintptr_t)value << 1) | 1U);
            else
                bottom = (BottomMinMaxOutside*)((((uintptr_t)value & bottomMask) << 1) | 1U);
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
            const uint32_t lowerBits = value & bottomMask;
            std::lock_guard bottomLock(bottomMutexes[upperBits]);
            bottom = new BottomMinMaxOutside;
            bottom->min = std::min(oldValue, lowerBits);
            bottom->max = std::max(oldValue, lowerBits);
            topLock.unlock();
            bottom->insert(oldValue);
            bottom->insert(lowerBits);
        } else {
            const uint32_t lowerBits = value & bottomMask;
            bottom->min = std::min(bottom->min, lowerBits); // must be updated with held topMutex
            bottom->max = std::max(bottom->max, lowerBits);
            std::lock_guard bottomLock(bottomMutexes[upperBits]);
            topLock.unlock();
            bottom->insert(lowerBits);
        }
    }

    /// Removes the given value. The value must be contained in this tree!
    void remove(const T sValue) noexcept {
        assert(((floatType && !std::isnan(sValue)) || (!floatType && sValue != MAX_T)) && "sValue is invalid!");
        const uint32_t value = convertToUnsigned(sValue);
        const uint32_t upperBits = value >> bottomK;
        BottomMinMaxOutside*& bottom = bottoms[upperBits];
        std::lock_guard topLock(topMutex);
        if (((uintptr_t)bottom & 1U) != 0) {
            bottom = nullptr;
            removeTop(upperBits);
        } else {
            // Actually, we do not really need to lock the mutex, we just need to wait until it is available
            std::lock_guard bottomLock(bottomMutexes[upperBits]);
            bottom->remove(value & bottomMask);
            if (bottom->min == bottom->max) {
                uint32_t oldValue = bottom->min;
                if constexpr (fullValueInPtr) {
                    oldValue |= value & ~bottomMask;
                }
                delete bottom;
                bottom = (BottomMinMaxOutside*)(((uintptr_t)oldValue << 1) | 1U);
            }
        }

        if (value == max) {
            if (topLevel3 == 0) {
                max = MIN_UT;
            } else {
                const uint32_t upperBits = maxTop();
                BottomMinMaxOutside* bottom = bottoms[upperBits];
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
        #ifdef VEB_32_LOCKED_TOP_USE_STD_MUTEX
        std::lock_guard topLock(topMutex);
        #else
        std::shared_lock topLock(topMutex);
        #endif
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
    // use a vector for topLevel2 if sizeof(VanEmdeBoas32LockedTop) is too large
    uint64_t topLevel2[topSizeLevel2 / sizeof(uint64_t)] = {0};
    uint64_t topLevel3 = 0;
    std::vector<BottomMinMaxOutside*> bottoms;

    #ifdef VEB_32_LOCKED_TOP_USE_STD_MUTEX
    mutable std::mutex topMutex; // also responsible for max, the pointers in bottoms and BottomMinMaxOutside::min/max
    #else
    mutable sf::contention_free_shared_mutex<> topMutex; // see above
    //mutable std::vector<sf::contention_free_shared_mutex<>> bottomMutexes; // does not seem to be worth it
    #endif
    mutable std::vector<std::mutex> bottomMutexes; // one for every of the possible 2^18 bottom data structures
    // std::shared_mutex seems to starve the writer. At least the "BackwardsInsertTwoProducers" test takes 10x more time.
};

} // namespace veb
