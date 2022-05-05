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
 * This implementation here is a version where insert and locate are lockless. Remove obtains a
 * unique lock.
 * 
 * You may define VEB_32_USE_BYTELL to use ska::bytell_hash_map by Malte Skarupke instead of the
 * specialized hash table. Furthermore, you may define VEB_32_LOCKLESS_USE_STD_SHARED_MUTEX
 * to use std::shared_mutex instead of sf::contention_free_shared_mutex<>
 * (https://github.com/AlexeyAB/object_threadsafe).
 */

#pragma once

#include <cstdint>
#include <cstddef>
#include <cassert>
#include <type_traits>
#include <thread>
#include <atomic>
#include <mutex>
#include <shared_mutex>

#include "van_emde_boas32_bottom.hpp"
#include "utils/utils.hpp"

#ifndef VEB_32_LOCKLESS_USE_STD_SHARED_MUTEX
#include "utils/safe_ptr.h"
#endif

#ifdef VEB_32_USE_BYTELL
#include "utils/bytell_hash_map.hpp"
#else
#include "hash_table.hpp"
#endif

namespace veb {

/// An optimized concurrent (not sequentially consistent) van Emde Boas tree for 32 bit integers or floats
/// in IEEE format. Implements a set and is optimized for 64 bit processors.
/// If T is integral, the value std::numeric_limits<T>::max() is reserved and cannot be inserted or located.
/// If T is float, NaN is reserved and cannot be inserted or located.
/// @tparam T The type to use. Must be uint32_t, int32_t or float in IEEE format.
template<class T = std::uint32_t>
class VanEmdeBoas32Lockless {
    using uint32_t = std::uint32_t; // using std::uint32_t inside a class is not valid
    using uint64_t = std::uint64_t;
    using size_t = std::size_t;
    static_assert(std::is_same_v<T, uint32_t> || std::is_same_v<T, int32_t> || std::is_same_v<T, float>,
        "T must be uint32_t, int32_t or float!");
    static_assert(sizeof(T) == 4 && (!std::is_same_v<T, float> || std::numeric_limits<T>::is_iec559),
        "Floating point type must be in 32 Bit IEE format!");

    class alignas(64) BottomAtomicMinMax { // sizeof(Bottom) <= 64 => fits in cache line
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
            uint32_t oldMin = min.load(std::memory_order_relaxed);
            while (value < oldMin && !min.compare_exchange_weak(oldMin, value, std::memory_order_relaxed))
                ;
            const uint32_t upperBits = value >> lowerK;
            top[upperBits >> upperLowerK] |= (1LLU << (upperBits & upperLowerMask));
            innerBottoms[upperBits] |= (1LLU << (value & lowerMask));
            uint32_t oldMax = max.load(std::memory_order_relaxed);
            while (value > oldMax && !max.compare_exchange_weak(oldMax, value, std::memory_order_release))
                ;
        }

        /// Removes the given value. The value and at least one other value must be contained!
        void remove(const uint32_t value) noexcept { // remove locks the tree => std::memory_order_relaxed
            assert(min.load(std::memory_order_relaxed) != max.load(std::memory_order_relaxed));
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
            if (value == min.load(std::memory_order_relaxed)) {
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
                min.store((upperResult << lowerK) | utils::countr_zero(innerBottoms.find(upperResult)->second),
                    std::memory_order_relaxed);
            } else if (value == max.load(std::memory_order_relaxed)) {
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
                max.store((upperResult << lowerK) | (63 - utils::countl_zero(innerBottoms.find(upperResult)->second)),
                    std::memory_order_relaxed);
            }
        }

        std::atomic<uint32_t> min = std::numeric_limits<uint32_t>::max(); // std::memory_order_relaxed is sufficient
        std::atomic<uint32_t> max = std::numeric_limits<uint32_t>::min();

    protected:
        uint64_t top[4] = {0};
        #ifdef VEB_32_USE_BYTELL
        ska::bytell_hash_map<uint32_t, uint64_t> innerBottoms;
        #else
        HashTable innerBottoms;
        #endif
    };
    
    // If true, all 32 bits fit in tagged pointer, otherwise only bottomK bits. True on 64 bit systems.
    static constexpr bool fullValueInPtr = sizeof(BottomAtomicMinMax*) >= 5;
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

    // 2 is never a valid BottomAtomicMinMax pointer
    // We cannot store a BottomAtomicMinMax pointer directly as constexpr since reinterpret_cast
    // to pointer is not possible here.
    static constexpr uintptr_t lockedBottom = 2U;

public:
    static constexpr int BitsK = 32;
    using TypeT = T;

    VanEmdeBoas32Lockless() {
        topLevel1 = new std::atomic<uint64_t>[topSizeLevel1 / sizeof(uint64_t)]();
        bottoms = new std::atomic<BottomAtomicMinMax*>[topSizeLevel1]();
    }

    ~VanEmdeBoas32Lockless() noexcept {
        uint32_t value = 0;
        while (topLevel3.load(std::memory_order_relaxed) != 0) {
            value = locateTop(value);
            BottomAtomicMinMax* bottom = bottoms[value].load(std::memory_order_relaxed);
            if (((uintptr_t)bottom & 1U) == 0)
                delete bottom;
            removeTop(value);
        }
        delete[] topLevel1;
        delete[] bottoms;
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
        std::atomic<BottomAtomicMinMax*>& bottomRef = bottoms[upperBits];
        BottomAtomicMinMax* bottom = loadBottom(bottomRef);
        if (bottom == nullptr) {
            // do nothing here
        } else if (((uintptr_t)bottom & 1U) != 0) {
            uint32_t bottomValue = (uint32_t)((uintptr_t)bottom >> 1);
            if constexpr (!fullValueInPtr) {
                bottomValue = (value & ~bottomMask) | bottomValue;
            }
            if (value <= bottomValue) {
                removeLock.unlock();
                return convertToSigned(bottomValue);
            }
            // do nothing here
        } else if ((value & bottomMask) > bottom->max.load(std::memory_order_acquire)) {
            // do nothing here
        } else {
            // locate value in bottom
            const uint32_t locateSearchValue = value & bottomMask;
            bottom = lockBottom(bottomRef);
            // We still know that bottom is a valid pointer since insert cannot change that fact
            // and remove locks the whole data structure.
            const uint32_t lowerBits = bottom->locate(locateSearchValue);
            unlockBottom(bottomRef, bottom);
            return convertToSigned((value & ~bottomMask) | lowerBits);
        }

        // locate top, return min of found bottom
        const uint32_t upperResult = locateTop(upperBits + 1);
        bottom = loadBottom(bottoms[upperResult]);
        uint32_t lowerResult;
        if (((uintptr_t)bottom & 1U) != 0) {
            lowerResult = (uint32_t)((uintptr_t)bottom >> 1);
            if constexpr (fullValueInPtr) {
                removeLock.unlock();
                return convertToSigned(lowerResult);
            }
        } else {
            lowerResult = bottom->min.load(std::memory_order_relaxed);
        }
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
        std::atomic<BottomAtomicMinMax*>& bottomRef = bottoms[upperBits];
        std::shared_lock removeLock(removeMutex);
        
        BottomAtomicMinMax* bottom = lockBottom(bottomRef);
        if (bottom == nullptr) {
            if constexpr (fullValueInPtr)
                unlockBottom(bottomRef, (BottomAtomicMinMax*)(((uintptr_t)value << 1) | 1U));
            else
                unlockBottom(bottomRef, (BottomAtomicMinMax*)((((uintptr_t)value & bottomMask) << 1) | 1U));
            insertTop(upperBits);
        } else if (((uintptr_t)bottom & 1U) != 0) {
            uint32_t oldValue = (uint32_t)((uintptr_t)bottom >> 1);
            if constexpr (fullValueInPtr) {
                if (oldValue == value) {
                    unlockBottom(bottomRef, bottom);
                    return; // value already contained
                }
                oldValue = oldValue & bottomMask; // only insert bottomK bits
            } else if (oldValue == (value & bottomMask)) {
                unlockBottom(bottomRef, bottom);
                return; // value already contained
            }
            bottom = new BottomAtomicMinMax;
            bottom->insert(oldValue);
            bottom->insert(value & bottomMask);
            unlockBottom(bottomRef, bottom);
        } else {
            bottom->insert(value & bottomMask);
            unlockBottom(bottomRef, bottom);
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
        std::atomic<BottomAtomicMinMax*>& bottom = bottoms[upperBits];
        std::lock_guard removeLock(removeMutex);
        BottomAtomicMinMax* bottomValue = bottom.load(std::memory_order_relaxed);
        if (((uintptr_t)bottomValue & 1U) != 0) {
            bottom.store(nullptr, std::memory_order_relaxed);
            removeTop(upperBits);
        } else {
            bottomValue->remove(value & bottomMask);
            if (bottomValue->min == bottomValue->max) {
                uint32_t oldValue = bottomValue->min;
                if constexpr (fullValueInPtr) {
                    oldValue |= value & ~bottomMask;
                }
                delete bottomValue;
                bottom.store((BottomAtomicMinMax*)(((uintptr_t)oldValue << 1) | 1U), std::memory_order_relaxed);
            }
        }

        if (value == max.load(std::memory_order_relaxed)) {
            if (topLevel3.load(std::memory_order_relaxed) == 0) {
                max.store(MIN_UT, std::memory_order_relaxed);
            } else {
                const uint32_t upperBits = maxTop();
                BottomAtomicMinMax* bottom = bottoms[upperBits].load(std::memory_order_relaxed);
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

    uint32_t locateTop(const uint32_t value) const noexcept {
        constexpr int bits = 6;
        constexpr uint32_t mask = (1U << bits) - 1;
        uint32_t valueShifted = value >> bits;
        uint64_t current = topLevel1[valueShifted].load(std::memory_order_acquire) & ~((1ULL << (value & mask)) - 1);
        if (current != 0)
            return (value & ~mask) | utils::countr_zero(current);
        ++valueShifted;
        const uint32_t valueShiftedTwice = valueShifted >> bits;
        current = topLevel2[valueShiftedTwice].load(std::memory_order_acquire)
            & ~((1ULL << (valueShifted & mask)) - 1);
        if (current != 0) {
            const uint32_t upperBits = (valueShifted & ~mask) | utils::countr_zero(current);
            return (upperBits << bits) | utils::countr_zero(topLevel1[upperBits].load(std::memory_order_acquire));
        }
        current = topLevel3.load(std::memory_order_acquire) & ~((1ULL << (valueShiftedTwice + 1)) - 1);
        assert(current != 0);
        uint32_t upperBits = utils::countr_zero(current);
        upperBits = (upperBits << bits) | utils::countr_zero(topLevel2[upperBits].load(std::memory_order_acquire));
        return (upperBits << bits) | utils::countr_zero(topLevel1[upperBits].load(std::memory_order_acquire));
    }

    void insertTop(const uint32_t value) noexcept {
        constexpr int bits = 6;
        constexpr uint32_t mask = (1U << bits) - 1;
        const uint32_t valueShifted = value >> bits;
        // difference to sequential: first set lower levels
        if (topLevel1[valueShifted].fetch_or((1LLU << (value & mask)), std::memory_order_release) == 0) {
            const uint32_t valueShiftedTwice = value >> (2 * bits);
            std::atomic<uint64_t>& level2 = topLevel2[valueShiftedTwice];
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
        std::atomic<uint64_t>& level1 = topLevel1[valueShifted];
        level1.store(level1.load(std::memory_order_relaxed) ^ (1LLU << (value & mask)), std::memory_order_relaxed);
        if (level1.load(std::memory_order_relaxed) != 0) return;
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
        return (result << bits) | (63 - utils::countl_zero(topLevel1[result].load(std::memory_order_relaxed)));
    }

    // Does not lock it.
    BottomAtomicMinMax* loadBottom(std::atomic<BottomAtomicMinMax*>& bottom) const noexcept {
        BottomAtomicMinMax* result;
        while ((uintptr_t)(result = bottom.load(std::memory_order_acquire)) == lockedBottom)
            std::this_thread::yield();
        return result;
    }

    BottomAtomicMinMax* lockBottom(std::atomic<BottomAtomicMinMax*>& bottom) const noexcept {
        BottomAtomicMinMax* locked = (BottomAtomicMinMax*)lockedBottom;
        BottomAtomicMinMax* result = bottom.load(std::memory_order_acquire);
        do {
            while (result == locked) {
                std::this_thread::yield();
                result = bottom.load(std::memory_order_acquire);
            }
        } while (!bottom.compare_exchange_weak(result, locked, std::memory_order_acquire));
        return result;
    }

    void unlockBottom(std::atomic<BottomAtomicMinMax*>& bottom, BottomAtomicMinMax* value) const noexcept {
        bottom.store(value, std::memory_order_release);
    }

    std::atomic<uint32_t> max = MIN_UT;
    std::atomic<uint64_t>* topLevel1; // not a vector because atomics are not trivially copyable
    // use a vector for topLevel2 if sizeof(VanEmdeBoas32LockedFineGrained) is too large
    std::atomic<uint64_t> topLevel2[topSizeLevel2 / sizeof(uint64_t)] = {0};
    std::atomic<uint64_t> topLevel3 = 0;
    std::atomic<BottomAtomicMinMax*>* bottoms; // not a vector because atomics are not trivially copyable

    #ifdef VEB_32_LOCKLESS_USE_STD_SHARED_MUTEX
    mutable std::shared_mutex removeMutex;
    #else
    mutable sf::contention_free_shared_mutex<> removeMutex;
    #endif
};

} // namespace veb
