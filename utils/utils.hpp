//  Copyright Dominik Bez 2022.

#pragma once

#include <random>
#include <vector>
#include <cassert>

namespace utils {

/** C++20 <bit>:
template<class T>
constexpr int countr_zero(const T value) noexcept {
    return std::countr_zero(value);
}
*/

// Undefined for value == 0
constexpr int countr_zero(const unsigned int value) noexcept {
    assert(value != 0 && "__builtin_ctz is undefined for 0!");
    return __builtin_ctz(value);
}
// Undefined for value == 0
constexpr int countr_zero(const unsigned long value) noexcept {
    assert(value != 0 && "__builtin_ctzl is undefined for 0!");
    return __builtin_ctzl(value);
}
// Undefined for value == 0
constexpr int countr_zero(const unsigned long long value) noexcept {
    assert(value != 0 && "__builtin_ctzll is undefined for 0!");
    return __builtin_ctzll(value);
}

/** C++20 <bit>:
template<class T>
constexpr int countl_zero(const T value) noexcept {
    return std::countl_zero(value);
}
*/

// Undefined for value == 0
constexpr int countl_zero(const unsigned int value) noexcept {
    assert(value != 0 && "__builtin_clz is undefined for 0!");
    return __builtin_clz(value);
}
// Undefined for value == 0
constexpr int countl_zero(const unsigned long value) noexcept {
    assert(value != 0 && "__builtin_clzl is undefined for 0!");
    return __builtin_clzl(value);
}
// Undefined for value == 0
constexpr int countl_zero(const unsigned long long value) noexcept {
    assert(value != 0 && "__builtin_clzll is undefined for 0!");
    return __builtin_clzll(value);
}


template<class T>
void shuffle(std::vector<T>& v, const int seed = 42) {
    std::mt19937_64 re(seed);
    std::shuffle(v.begin(), v.end(), re);
}

} // namespace utils