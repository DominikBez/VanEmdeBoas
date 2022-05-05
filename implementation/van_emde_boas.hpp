/* Copyright Dominik Bez 2022.
 *
 * Implementation of a generic Van Emde Boas set for integers and floating point values. The
 * recursion ends if at most 6 bits are left. In this case, a single 64 bit unsigned integer
 * is sufficient to store all 2^6 = 64 values by using the value as an bit index and using
 * least/most significant bit operations to find minimum or maximum values.
 * 
 * The recursion for the bottom data structure ends as well if it only contains a single value.
 * In this case, the value is stored in the pointer itself, if possible, or on the heap directly
 * without the bottom data structure otherwise. The least significant bit of the pointer is set
 * to 1 to indicate a single value. Note that by alignment of the Bottom data structure as well
 * as the values which may be stored directly on the heap, the least significant bit of the
 * pointer is always 0 for a valid pointer.
 * 
 * You may define VEB_USED_HASHMAP to 0 to use std::unordered_map, to 1 to use
 * ska::bytell_hash_map by Malte Skarupke (this is the standard behavior) or to 2
 * to use ska::flat_hash_map by Malte Skarupke.
 */

#pragma once

#include <cstdint>
#include <type_traits>
#include <limits>
#include <cassert>
#include <algorithm>

#ifndef NDEBUG
#include <cmath>
#endif

#include "utils/utils.hpp"

#define VEB_USE_STD_UNORDERED_MAP 0
#define VEB_USE_SKA_BYTELL_HASH_MAP 1
#define VEB_USE_SKA_FLAT_HASH_MAP 2

#ifndef VEB_USED_HASHMAP
#define VEB_USED_HASHMAP VEB_USE_SKA_BYTELL_HASH_MAP
#endif

#if VEB_USED_HASHMAP == VEB_USE_SKA_BYTELL_HASH_MAP
#include "utils/bytell_hash_map.hpp"
#elif VEB_USED_HASHMAP == VEB_USE_SKA_FLAT_HASH_MAP
#include "utils/flat_hash_map.hpp"
#else
#include <unordered_map>
#endif

namespace veb {

inline constexpr int BASE_THRESHOLD = 6;

template<int K>
struct selectType {
    static_assert(K > 0, "K must be greater than zero!");
    static_assert(K <= 128, "Use custom type T for K > 128!");
    #ifndef __SIZEOF_INT128__
        static_assert(K <= 64, "Use custom type T for K > 64 if your compiler doesn't support __uint128_t!"),
    #endif

    typedef typename std::conditional<(K > 64), __uint128_t,
            typename std::conditional<(K > 32), uint64_t,
            typename std::conditional<(K > 16), uint32_t,
            typename std::conditional<(K > 8), uint16_t, uint8_t>::type>::type>::type>::type type;
    // It is possible to add wider types here.
};

template<bool> struct Range;

/**
 * A generic van Emde Boas tree implementing a set. The type T must be an integral or float type with at least K bits.
 * If T is not explicitly set, it is an unsigned type with at least K bits.
 * If T is signed, it must use the two's complement. This should always be true for the standard types,
 * especially for C++20. Furthermore, if T is signed or a floating point type, K must be equal to sizeof(T) * 8.
 * If T is integral, the value std::numeric_limits<T>::max() is reserved and cannot be inserted or located.
 * If T is floating point, NaN is reserved and cannot be inserted or located.
 * Only float and double in IEEE format are supported. Other formats and long double may not work correctly.
 * 
 * @tparam K The number of bits to use in the values. Must be equal to 8*sizeof(K) if T is signed or floating point.
 *           If K > 128 (or > 256 if the compiler supports __uint128_t), appropriate unsigned integer types of all
 *           powers of two from 128 (or 256 with __uint128_t) to K/2 must be added to veb::selectType::type;
 * @tparam T The type to use. Must be integral or floating point and have at least K bits.
 */
template<int K, class T = typename selectType<K>::type, class = Range<true>>
class alignas(16) VanEmdeBoas {
    template<int FriendK, class FriendT, class>
    friend class VanEmdeBoas;
    // The following assert may fail for T = __uint128_t in GCC with std=C++17
    static_assert(std::is_unsigned<T>::value || K == sizeof(T) * 8,
        "The type T must be an unsigned integral type or K must be the number of bits in T!");
    static_assert(8 * sizeof(T) >= K, "T must have at least K bits!");

    static constexpr bool recursionBase = false;
    static constexpr bool signedType = std::is_signed<T>::value;
    static constexpr bool floatType = std::is_floating_point<T>::value;
    static_assert(!floatType || std::numeric_limits<T>::is_iec559, "Floating point type must be in IEE format!");
    using UT = typename std::conditional<signedType, typename selectType<8 * sizeof(T)>::type, T>::type;
    static_assert(sizeof(T) == sizeof(UT), "sizeof(T) must be equal to sizeof(UT)!");
    static constexpr auto MAX_T = std::numeric_limits<T>::max();
    static constexpr auto INVALID = floatType ? std::numeric_limits<T>::quiet_NaN() : MAX_T;
    static constexpr auto MAX_UT = std::numeric_limits<UT>::max();
    static constexpr auto MIN_UT = std::numeric_limits<UT>::min();
    static constexpr UT signBit = UT(1) << sizeof(T) * 8 - 1;

    static constexpr int bottomK = K / 2;
    static constexpr UT bottomMask = (UT(1) << bottomK) - 1;
    using Bottom = VanEmdeBoas<bottomK>;
    using BottomUT = typename Bottom::UT;
    static constexpr bool storeBottomInline = Bottom::recursionBase;
    using BottomPtr = typename std::conditional<storeBottomInline, Bottom, Bottom*>::type;
    static constexpr int topK = (K + 1) / 2;
    using Top = VanEmdeBoas<topK>;
    using TopUT = typename Top::UT;
    static constexpr bool storeSingleValueInPtr = bottomK < sizeof(Bottom*) * 8;
    static_assert(storeSingleValueInPtr || alignof(BottomUT) >= 2,
        "Alignment of BottomUT must be at least 2 for pointer tagging!");

public:
    static constexpr int BitsK = K;
    using TypeT = T;

    ~VanEmdeBoas() noexcept {
        assert((((uintptr_t)this & 1ULL) == 0) && "Least significant bit of address must be zero for pointer tagging!");
        if constexpr (!storeBottomInline) {
            for(auto& [key, value] : bottoms) {
                if constexpr (storeSingleValueInPtr) {
                    if (((uintptr_t)value & 1U) == 0)
                        delete value;
                } else {
                    if (((uintptr_t)value & 1U) == 0) {
                        delete value;
                    } else {
                        delete (BottomUT*)((uintptr_t)value & ~(uintptr_t)1U);
                    }
                }
            }
        }
    }

    /// If T is integral: returns std::numeric_limits<T>::max() if value > max or the tree is empty.
    /// If T is floating point: returns NaN if value > max or the tree is empty.
    /// @param value The value to search for. If T is integral, value must not be std::numeric_limits<T>::max()
    ///              If T is floating point, value must not be NaN.
    /// @return The smallest contained number not smaller (i.e., greater or equal) than value.
    T locate(const T value) const noexcept {
        assert(((floatType && !std::isnan(value)) || (!floatType && value != MAX_T)) && "value is invalid!");
        const UT uValue = convertToUnsigned(value);
        if (uValue > max || empty()) return INVALID;
        return convertToSigned(locateInner(uValue));
    }
    /// Alias to locate
    T lower_bound(const T value) const noexcept {
        return locate(value);
    }

    /// Inserts the given value. Does nothing if the value is already contained.
    void insert(const T value) noexcept {
        assert(((floatType && !std::isnan(value)) || (!floatType && value != MAX_T)) && "value is invalid!");
        const UT uValue = convertToUnsigned(value);
        insertInner(uValue);
    }

    /// Removes the given value. The value must be contained in this tree!
    void remove(const T value) noexcept {
        assert(((floatType && !std::isnan(value)) || (!floatType && value != MAX_T)) && "value is invalid!");
        const UT uValue = convertToUnsigned(value);
        removeInner(uValue);
    }
    /// Alias to remove
    void erase(const T value) noexcept {
        remove(value);
    }

    bool empty() const noexcept {
        return max < min;
    }
    
private:
    static constexpr UT convertToUnsigned(const T value) noexcept {
        if constexpr (floatType) {
            return reinterpret_cast<const UT&&>(std::move(value)) ^ signBit; // just flip the sign
        } else if constexpr (signedType) {
            return reinterpret_cast<const UT&&>(value) ^ signBit; // just flip the sign
        } else {
            return value;
        }
    }
    static constexpr T convertToSigned(const UT value) noexcept {
        if constexpr (signedType) {
            return reinterpret_cast<const T&&>(value ^ signBit); // just flip the sign
        } else {
            return value;
        }
    }

    static constexpr BottomUT getMin(const BottomPtr& bottom) noexcept {
        if constexpr (storeBottomInline) {
            return bottom.getMin();
        } else {
            if (((uintptr_t)bottom & 1U) != 0) {
                if constexpr (storeSingleValueInPtr) {
                    return (BottomUT)((uintptr_t)bottom >> 1);
                } else {
                    return *(BottomUT*)((uintptr_t)bottom & ~(uintptr_t)1U);
                }
            } else {
                return bottom->getMin();
            }
        }
    }

    UT getMin() const noexcept {
        return min;
    }
    UT getMax() const noexcept {
        return max;
    }

    template<bool outerCall = false>
    UT locateInner(const UT value) const noexcept {
        #pragma GCC diagnostic ignored "-Wshift-count-overflow" // because of the static_assert above and the || below
        assert(((K == 8*sizeof(UT)) || (value & ~((UT(1) << K) - 1)) == 0) && "Must only use the lower K bits!");
        #pragma GCC diagnostic pop

        if constexpr (outerCall) {
            // The second check is necessary since value may be 0.
            if (value > max || empty()) return MAX_UT;
        } else {
            assert(!empty());
            assert(value <= max);
        }

        TopUT upperBits = value >> bottomK;
        BottomUT lowerBits;
        auto bottomIt = bottoms.find(upperBits);

        if (bottomIt == bottoms.end()) {
            upperBits = top.locateInner(upperBits + 1);
            lowerBits = getMin(bottoms.find(upperBits)->second);
        } else {
            lowerBits = BottomUT(value & bottomMask);
            const BottomPtr& bottom = bottomIt->second;
            if constexpr (storeBottomInline) {
                if (lowerBits > bottom.getMax()) {
                    upperBits = top.locateInner(upperBits + 1);
                    lowerBits = getMin(bottoms.find(upperBits)->second);
                } else {
                    lowerBits = bottom.locateInner(lowerBits);
                }
            } else {
                if (((uintptr_t)bottom & 1U) != 0) {
                    BottomUT bottomValue;
                    if constexpr (storeSingleValueInPtr) {
                        bottomValue = (BottomUT)((uintptr_t)bottom >> 1);
                    } else {
                        bottomValue = *(BottomUT*)((uintptr_t)bottom & ~(uintptr_t)1U);
                    }
                    if (lowerBits > bottomValue) {
                        upperBits = top.locateInner(upperBits + 1);
                        lowerBits = getMin(bottoms.find(upperBits)->second);
                    } else {
                        lowerBits = bottomValue;
                    }
                } else {
                    if (lowerBits > bottom->getMax()) {
                        upperBits = top.locateInner(upperBits + 1);
                        lowerBits = getMin(bottoms.find(upperBits)->second);
                    } else {
                        lowerBits = bottom->locateInner(lowerBits);
                    }
                }
            }
        }

        return ((UT)upperBits << bottomK) | lowerBits;
    }

    void insertInner(const UT value) noexcept {
        #pragma GCC diagnostic ignored "-Wshift-count-overflow" // because of the static_assert above and the || below
        assert(((K == 8*sizeof(UT)) || (value & ~((UT(1) << K) - 1)) == 0) && "Must only use the lower K bits!");
        #pragma GCC diagnostic pop

        max = std::max(max, value);
        min = std::min(min, value);
        const TopUT upperBits = value >> bottomK;
        const BottomUT lowerBits = value & bottomMask;
        auto& bottom = bottoms[upperBits];
        if constexpr (storeBottomInline) {
            constexpr Bottom emptyBottom = Bottom();
            if (bottom == emptyBottom)
                top.insertInner(upperBits);
            bottom.insertInner(lowerBits);
        } else {
            if (bottom == nullptr) {
                top.insertInner(upperBits);
                if constexpr (storeSingleValueInPtr) {
                    bottom = (Bottom*)(((uintptr_t)lowerBits << 1) | 1U);
                } else {
                    BottomUT *newBottom = new BottomUT(lowerBits);
                    bottom = (Bottom*)((uintptr_t)newBottom | 1U);
                }
            } else if (((uintptr_t)bottom & 1U) != 0) {
                BottomUT oldValue;
                if constexpr (storeSingleValueInPtr) {
                    oldValue = (BottomUT)((uintptr_t)bottom >> 1);
                } else {
                    BottomUT *oldBottom = (BottomUT*)((uintptr_t)bottom & ~(uintptr_t)1U);
                    oldValue = *oldBottom;
                    delete oldBottom;
                }
                bottom = new Bottom;
                bottom->insertInner(oldValue);
                bottom->insertInner(lowerBits);
                // Inserting them one by one may cause reallocation on several recursion levels
                // if !storeSingleValueInPtr, especially if oldValue and lowerBits are close to
                // each other. The problem could be mitigated by having a bulk insert to insert
                // both values at once. This way, saving oldValue directly on the heap which is
                // deleted immediately after to create a new VanEmdeBoas containing both values
                // can be avoided. However, since this is only a problem for K > 128 on a 64 bit
                // system, I spare me the time to implement a bulk insert.
            } else {
                bottom->insertInner(lowerBits);
            }
        }
    }

    void removeInner(const UT value) noexcept {
        #pragma GCC diagnostic ignored "-Wshift-count-overflow" // because of the static_assert above and the || below
        assert(((K == 8*sizeof(UT)) || (value & ~((UT(1) << K) - 1)) == 0) && "Must only use the lower K bits!");
        #pragma GCC diagnostic pop
        assert(locateInner<true>(value) == value && "A value must be contained to be removed!");

        const TopUT upperBits = value >> bottomK;
        auto& bottom = bottoms[upperBits];
        if constexpr (storeBottomInline) {
            bottom.removeInner(value & bottomMask);
            if (bottom.empty()) {
                bottoms.erase(upperBits);
                top.removeInner(upperBits);
            }
        } else {
            if (((uintptr_t)bottom & 1U) != 0) {
                if constexpr (!storeSingleValueInPtr) {
                    delete (BottomUT*)((uintptr_t)bottom & ~(uintptr_t)1U);
                }
                bottoms.erase(upperBits);
                top.removeInner(upperBits);
            } else {
                bottom->removeInner(value & bottomMask);
                if (bottom->empty()) {
                    // This is necessary since we do not delete a bottom tree after a remove if it only contains one element
                    delete bottom;
                    bottoms.erase(upperBits);
                    top.removeInner(upperBits);
                }
            }
        }

        if (value == min) {
            if (min == max) {
                min = MAX_UT;
                max = MIN_UT;
            } else {
                const TopUT upperBits = top.getMin();
                auto& bottom = bottoms[upperBits];
                BottomUT lowerBits;
                if constexpr (storeBottomInline) {
                    lowerBits = bottom.getMin();
                } else {
                    if (((uintptr_t)bottom & 1U) != 0) {
                        if constexpr (storeSingleValueInPtr) {
                            lowerBits = (BottomUT)((uintptr_t)bottom >> 1);
                        } else {
                            lowerBits = *(BottomUT*)((uintptr_t)bottom & ~(uintptr_t)1U);
                        }
                    } else {
                        lowerBits = bottom->getMin();
                    }
                }
                min = ((UT)upperBits << bottomK) | lowerBits;
            }
        } else if (value == max) {
            const TopUT upperBits = top.getMax();
            auto& bottom = bottoms[upperBits];
            BottomUT lowerBits;
            if constexpr (storeBottomInline) {
                lowerBits = bottom.getMax();
            } else {
                if (((uintptr_t)bottom & 1U) != 0) {
                    if constexpr (storeSingleValueInPtr) {
                        lowerBits = (BottomUT)((uintptr_t)bottom >> 1);
                    } else {
                        lowerBits = *(BottomUT*)((uintptr_t)bottom & ~(uintptr_t)1U);
                    }
                } else {
                    lowerBits = bottom->getMax();
                }
            }
            max = ((UT)upperBits << bottomK) | lowerBits;
        }

        assert(locateInner<true>(0) == min);
        assert(empty() || locateInner<true>(min) == min);
        assert(empty() || locateInner<true>(max) == max);
        #pragma GCC diagnostic ignored "-Wshift-count-overflow" // because of the static_assert above and the || below
        assert(max == MAX_UT || !((K == 8*sizeof(UT)) || ((max + 1) & ~((UT(1) << K) - 1)) == 0)
                || locateInner<true>(max + 1) == MAX_UT);
        #pragma GCC diagnostic pop
    }

    UT min = MAX_UT;
    UT max = MIN_UT;
    Top top;

#if VEB_USED_HASHMAP == VEB_USE_SKA_BYTELL_HASH_MAP
    ska::bytell_hash_map<TopUT, BottomPtr> bottoms;
#elif VEB_USED_HASHMAP == VEB_USE_SKA_FLAT_HASH_MAP
    ska::flat_hash_map<TopUT, BottomPtr> bottoms;
#else
    std::unordered_map<TopUT, BottomPtr> bottoms;
#endif
};

// Recursion base template specialization. Consists of a single integer.
template<int K, class T>
class VanEmdeBoas<K, T, Range<(K <= BASE_THRESHOLD)>> {
    template<int FriendK, class FriendT, class>
    friend class VanEmdeBoas;
    // The following assert may fail for T = __uint128_t in GCC with std=C++17
    static_assert(std::is_unsigned<T>::value || K == sizeof(T) * 8,
        "The type T must be an unsigned integral type or K must be the number of bits in T!");
    static_assert(8 * sizeof(T) >= K && "T must have at least K bits!");

    static constexpr bool recursionBase = true;
    static constexpr bool signedType = std::is_signed<T>::value;
    using UT = typename std::conditional<signedType, typename std::make_unsigned<T>::type, T>::type;
    static constexpr int storageBits = 1 << K;
    using Storage = typename selectType<storageBits>::type;
    using StorageAtLeast32 = typename std::conditional<(storageBits < 32), std::uint32_t, Storage>::type;
    static constexpr auto MAX_T = std::numeric_limits<T>::max();
    static constexpr auto MAX_UT = std::numeric_limits<UT>::max();
    static constexpr auto MIN_UT = std::numeric_limits<UT>::min();
    static constexpr UT signBit = UT(1) << sizeof(T) * 8 - 1;

public:
    static constexpr int BitsK = K;
    using TypeT = T;

    T locate(const T value) const noexcept {
        assert((value != MAX_T) && "value is invalid!");
        const UT uValue = convertToUnsigned(value);
        return convertToSigned(locateInner<true>(uValue));
    }
    T lower_bound(const T value) const noexcept { // alias to locate
        return locate(value);
    }

    void insert(const T value) noexcept {
        assert((value != MAX_T) && "value is invalid!");
        const UT uValue = convertToUnsigned(value);
        insertInner(uValue);
    }

    void remove(const T value) noexcept {
        assert((value != MAX_T) && "value is invalid!");
        const UT uValue = convertToUnsigned(value);
        removeInner(uValue);
    }
    void erase(const T value) noexcept { // alias to remove
        remove(value);
    }

    bool empty() const noexcept {
        return storage == 0;
    }

    bool operator==(const VanEmdeBoas<K, T, Range<(K <= BASE_THRESHOLD)>>& other) {
        return storage == other.storage;
    }
    
private:
    static constexpr UT convertToUnsigned(const T value) noexcept {
        if constexpr (signedType) {
            return reinterpret_cast<const UT&&>(value) ^ signBit; // just flip the sign
        } else {
            return value;
        }
    }
    static constexpr T convertToSigned(const UT value) noexcept {
        if constexpr (signedType) {
            return reinterpret_cast<const T&&>(value ^ signBit); // just flip the sign
        } else {
            return value;
        }
    }

    static constexpr UT getMin(const Storage value) noexcept {
        return utils::countr_zero((StorageAtLeast32)value);
    }

    UT getMin() const noexcept {
        return getMin(storage);
    }
    UT getMax() const noexcept {
        constexpr UT bits = 8 * sizeof(StorageAtLeast32);
        return bits - 1 - utils::countl_zero((StorageAtLeast32)storage);
    }

    template<bool outerCall = false>
    T locateInner(const T value) const noexcept {
        #pragma GCC diagnostic ignored "-Wshift-count-overflow" // because of the static_assert above and the || below
        assert(((K == 8*sizeof(UT)) || (value & ~((UT(1) << K) - 1)) == 0) && "Must only use the lower K bits!");
        #pragma GCC diagnostic pop

        Storage maskedStorage = storage & ~((Storage(1) << (value)) - 1);
        if constexpr (outerCall) {
            if (maskedStorage == 0) return MAX_UT;
        } else {
            assert(maskedStorage != 0);
        }
        return getMin(maskedStorage);
    }

    void insertInner(const UT value) noexcept {
        assert(value < storageBits && "value must be at most 2^K - 1");

        storage |= Storage(1) << value;
    }

    void removeInner(const T value) noexcept {
        #pragma GCC diagnostic ignored "-Wshift-count-overflow" // because of the static_assert above and the || below
        assert(((K == 8*sizeof(UT)) || (value & ~((UT(1) << K) - 1)) == 0) && "Must only use the lower K bits!");
        #pragma GCC diagnostic pop

        storage &= ~(Storage(1) << value);
    }

    Storage storage = 0;
};

} // namespace veb
