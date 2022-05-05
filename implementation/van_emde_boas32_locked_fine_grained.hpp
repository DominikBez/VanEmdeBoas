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
 * This implementation here is a version where the top data structure is locked fine-grained.
 * 
 * You may define VEB_32_USE_BYTELL to use ska::bytell_hash_map by Malte Skarupke instead of the
 * specialized hash table. Furthermore, you may define VEB_32_LOCKED_FINE_GRAINED_USE_STD_SHARED_MUTEX
 * to use std::shared_mutex instead of sf::contention_free_shared_mutex<>
 * (https://github.com/AlexeyAB/object_threadsafe).
 */

#pragma once

#include <cstdint>
#include <cstddef>
#include <cassert>
#include <type_traits>
#include <atomic>
#include <mutex>
#include <shared_mutex>

#include "van_emde_boas32_bottom.hpp"
#include "utils/utils.hpp"

#ifndef VEB_32_LOCKED_FINE_GRAINED_USE_STD_SHARED_MUTEX
#include "utils/safe_ptr.h"
#endif

namespace veb {

/// An optimized concurrent (not sequentially consistent) van Emde Boas tree for 32 bit integers or floats
/// in IEEE format. Implements a set and is optimized for 64 bit processors.
/// If T is integral, the value std::numeric_limits<T>::max() is reserved and cannot be inserted or located.
/// If T is float, NaN is reserved and cannot be inserted or located.
/// @tparam T The type to use. Must be uint32_t, int32_t or float in IEEE format.
template<class T = std::uint32_t>
class VanEmdeBoas32LockedFineGrained {
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

    VanEmdeBoas32LockedFineGrained()
        : topLevel1(topSizeLevel1 / sizeof(uint64_t), 0),
          bottoms(topSizeLevel1, nullptr),
          topMutexes(topSizeLevel1 / sizeof(uint64_t))
        { }

    ~VanEmdeBoas32LockedFineGrained() noexcept {
        uint32_t value = 0;
        std::unique_lock<std::mutex> lockDummy;
        while (topLevel3.load(std::memory_order_relaxed) != 0) {
            value = locateTop<false>(value, lockDummy);
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
        std::shared_lock removeLock(removeMutex);
        if (value > max.load(std::memory_order_acquire) || topLevel3.load(std::memory_order_acquire) == 0)
            return INVALID;

        const uint32_t upperBits = value >> bottomK;
        const Bottom* bottom = bottoms[upperBits];
        const uint32_t topMutexIndex = upperBits >> 6;
        std::unique_lock topLock(topMutexes[topMutexIndex]);
        if (bottom == nullptr) {
            // do nothing here
        } else if (((uintptr_t)bottom & 1U) != 0) {
            uint32_t bottomValue = (uint32_t)((uintptr_t)bottom >> 1);
            if constexpr (!fullValueInPtr) {
                bottomValue = (value & ~bottomMask) | bottomValue;
            }
            if (value <= bottomValue) {
                topLock.unlock();
                removeLock.unlock();
                return convertToSigned(bottomValue);
            }
            // do nothing here
        } else if ((value & bottomMask) > bottom->max) {
            // do nothing here
        } else {
            // locate value in bottom
            return convertToSigned((value & ~bottomMask) | bottom->locate(value & bottomMask));
        }

        // locate top, return min of found bottom
        const uint32_t nextUpperBits = upperBits + 1;
        const uint32_t nextTopMutexIndex = nextUpperBits >> 6;
        if (nextTopMutexIndex != topMutexIndex) {
            std::unique_lock(topMutexes[nextTopMutexIndex]).swap(topLock);
        }
        const uint32_t upperResult = locateTop(nextUpperBits, topLock);
        bottom = bottoms[upperResult];
        uint32_t lowerResult;
        if (((uintptr_t)bottom & 1U) != 0) {
            lowerResult = (uint32_t)((uintptr_t)bottom >> 1);
            if constexpr (fullValueInPtr) {
                topLock.unlock();
                removeLock.unlock();
                return convertToSigned(lowerResult);
            }
        } else {
            lowerResult = bottom->min;
        }
        topLock.unlock();
        removeLock.unlock();
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
        Bottom*& bottom = bottoms[upperBits];
        std::shared_lock removeLock(removeMutex);
        std::unique_lock topLock(topMutexes[upperBits >> 6]);
        if (bottom == nullptr) {
            if constexpr (fullValueInPtr)
                bottom = (Bottom*)(((uintptr_t)value << 1) | 1U);
            else
                bottom = (Bottom*)((((uintptr_t)value & bottomMask) << 1) | 1U);
            insertTop(upperBits, topLock);
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
        uint32_t oldMax = max.load(std::memory_order_acquire);
        while (value > oldMax && !max.compare_exchange_weak(oldMax, value, std::memory_order_release))
            ;
    }

    /// Removes the given value. The value must be contained in this tree!
    /// remove locks the tree completely.
    void remove(const T sValue) noexcept {
        assert(((floatType && !std::isnan(sValue)) || (!floatType && sValue != MAX_T)) && "sValue is invalid!");
        const uint32_t value = convertToUnsigned(sValue);
        const uint32_t upperBits = value >> bottomK;
        Bottom*& bottom = bottoms[upperBits];
        std::lock_guard removeLock(removeMutex);
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

        if (value == max.load(std::memory_order_relaxed)) {
            if (topLevel3.load(std::memory_order_relaxed) == 0) {
                max.store(MIN_UT, std::memory_order_relaxed);
            } else {
                const uint32_t upperBits = maxTop();
                Bottom* bottom = bottoms[upperBits];
                if (((uintptr_t)bottom & 1U) != 0) {
                    const uint32_t lowerBits = (uint32_t)((uintptr_t)bottom >> 1);
                    if constexpr (fullValueInPtr) {
                        max.store(lowerBits, std::memory_order_relaxed);
                    } else {
                        max.store((upperBits << bottomK) | lowerBits, std::memory_order_relaxed);
                    }
                } else {
                    max.store((upperBits << bottomK) | bottom->max, std::memory_order_relaxed);
                }
            }
        }
    }
    /// Alias to remove
    void erase(const T value) noexcept {
        remove(value);
    }

    bool empty() const noexcept {
        std::shared_lock removeLock(removeMutex);
        return topLevel3.load(std::memory_order_acquire) == 0;
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

    // The corresponding top mutex must be locked and may be swapped to lock the correct mutex!
    template<bool useLock = true>
    uint32_t locateTop(const uint32_t value, std::unique_lock<std::mutex>& topLock) const noexcept {
        constexpr int bits = 6;
        constexpr uint32_t mask = (1U << bits) - 1;
        uint32_t valueShifted = value >> bits;
        uint64_t current = topLevel1[valueShifted] & ~((1ULL << (value & mask)) - 1);
        if (current != 0)
            return (value & ~mask) | utils::countr_zero(current);
        if constexpr (useLock)
            topLock.unlock();
        ++valueShifted;
        const uint32_t valueShiftedTwice = valueShifted >> bits;
        current = topLevel2[valueShiftedTwice].load(std::memory_order_acquire)
            & ~((1ULL << (valueShifted & mask)) - 1);
        if (current != 0) {
            const uint32_t upperBits = (valueShifted & ~mask) | utils::countr_zero(current);
            if constexpr (useLock)
                std::unique_lock(topMutexes[upperBits]).swap(topLock);
            return (upperBits << bits) | utils::countr_zero(topLevel1[upperBits]);
        }
        current = topLevel3.load(std::memory_order_acquire) & ~((1ULL << (valueShiftedTwice + 1)) - 1);
        assert(current != 0);
        uint32_t upperBits = utils::countr_zero(current);
        upperBits = (upperBits << bits) | utils::countr_zero(topLevel2[upperBits].load(std::memory_order_acquire));
        if constexpr (useLock)
            std::unique_lock(topMutexes[upperBits]).swap(topLock);
        return (upperBits << bits) | utils::countr_zero(topLevel1[upperBits]);
    }

    // The corresponding top mutex must be locked and will be unlocked!
    void insertTop(const uint32_t value, std::unique_lock<std::mutex>& topLock) noexcept {
        constexpr int bits = 6;
        constexpr uint32_t mask = (1U << bits) - 1;
        const uint32_t valueShifted = value >> bits;
        uint64_t& level1 = topLevel1[valueShifted];
        uint64_t oldLevel1 = level1;
        level1 |= (1LLU << (value & mask));
        topLock.unlock();
        if (oldLevel1 == 0) {
            const uint32_t valueShiftedTwice = value >> (2 * bits);
            std::atomic<uint64_t>& level2 = topLevel2[valueShiftedTwice];
            // difference to sequential: first set level 2
            if (level2.fetch_or(1LLU << (valueShifted & mask), std::memory_order_release) == 0) {
                topLevel3.fetch_or(1LLU << valueShiftedTwice, std::memory_order_release);
            }
        }
    }

    /// Removes the value from the top data structure. Only the lower 18 bits of value are
    /// considered and the value is assumed to be contained in the top data structure.
    /// The top tree must be locked for this!
    void removeTop(const uint32_t value) noexcept {
        // only used in remove and destructor => std::memory_order_relaxed
        constexpr int bits = 6;
        constexpr uint32_t mask = (1U << bits) - 1;
        const uint32_t valueShifted = value >> bits;
        uint64_t& level1 = topLevel1[valueShifted];
        level1 ^= (1LLU << (value & mask));
        if (level1 != 0) return;
        const uint32_t valueShiftedTwice = value >> (2 * bits);
        std::atomic<uint64_t>& level2 = topLevel2[valueShiftedTwice];
        level2.store(level2.load(std::memory_order_relaxed) ^ (1LLU << (valueShifted & mask)),
                     std::memory_order_relaxed);
        if (level2 != 0) return;
        topLevel3.store(topLevel3.load(std::memory_order_relaxed) ^ (1LLU << valueShiftedTwice),
                        std::memory_order_relaxed);
    }

    /// Returns the maximum value in the top data structure. Assumes it is not empty!
    /// The top tree must be locked for this!
    uint32_t maxTop() const noexcept { // only used in remove => std::memory_order_relaxed
        constexpr int bits = 6;
        uint32_t result = 63 - utils::countl_zero(topLevel3.load(std::memory_order_relaxed));
        result = (result << bits) | (63 - utils::countl_zero(topLevel2[result].load(std::memory_order_relaxed)));
        return (result << bits) | (63 - utils::countl_zero(topLevel1[result]));
    }

    std::atomic<uint32_t> max = MIN_UT;
    std::vector<uint64_t> topLevel1;
    // use a vector for topLevel2 if sizeof(VanEmdeBoas32LockedFineGrained) is too large
    std::atomic<uint64_t> topLevel2[topSizeLevel2 / sizeof(uint64_t)] = {0};
    std::atomic<uint64_t> topLevel3 = 0;
    std::vector<Bottom*> bottoms;

    mutable std::vector<std::mutex> topMutexes; // one for every uint64_t in topLevel1
    #ifdef VEB_32_LOCKED_FINE_GRAINED_USE_STD_SHARED_MUTEX
    mutable std::shared_mutex removeMutex;
    #else
    mutable sf::contention_free_shared_mutex<> removeMutex;
    #endif
};

} // namespace veb
