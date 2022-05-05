// Copyright Dominik Bez 2022.

#pragma once

#include <vector>
#include <cstddef>
#include <random>
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <cmath>
#include <array>

template<class T, class Gen>
std::vector<T> generateRandomVector(const std::size_t size, Gen& gen) {
    std::vector<T> numbers(size);
    for (std::size_t i = 0; i < size; ++i) {
        numbers[i] = gen.rand();
    }
    return numbers;
}

// Might run indefinitely!
template<class T, class Gen>
std::vector<T> generateUniqueRandomVector(const std::size_t size, Gen& gen) {
    std::unordered_set<T> set;
    std::vector<T> numbers(size);
    for (std::size_t i = 0; i < size; ++i) {
        do {
            numbers[i] = gen.rand();
        } while (set.find(numbers[i]) != set.end());
        set.insert(numbers[i]);
    }
    return numbers;
}

template<class T>
class UniformIntGenerator {
public:
    UniformIntGenerator(const T seed = 42, const T min = std::numeric_limits<T>::min(),
                        const T max = std::numeric_limits<T>::max())
        : re(seed), dis(min, max) { }
    
    T rand() noexcept {
        return dis(re);
    }

private:
    std::mt19937_64 re;
    std::uniform_int_distribution<T> dis;
};

template<class T>
std::vector<T> generateUniformIntVector(const std::size_t size, const T seed = 42, 
        const T min = std::numeric_limits<T>::min(),
        const T max = std::numeric_limits<T>::max()) {
    UniformIntGenerator<T> gen(seed, min, max);
    return generateRandomVector<T>(size, gen);
}

// Generates a vector with size random unique elements. May run indefinitely!
template<class T>
std::vector<T> generateUniqueUniformIntVector(const std::size_t size, const T seed = 42, 
        const T min = std::numeric_limits<T>::min(),
        const T max = std::numeric_limits<T>::max()) {
    UniformIntGenerator<T> gen(seed, min, max);
    return generateUniqueRandomVector<T>(size, gen);
}

// Generates clusters of succeeding integers
template<class T>
class ClusterIntGenerator {
public:
    ClusterIntGenerator(const T seed = 42, const T min = std::numeric_limits<T>::min(),
                        const T max = std::numeric_limits<T>::max(), const std::size_t clusterSize = 1000)
        : re(seed), dis(min, max), min(min), max(max), clusterSize(clusterSize), current(min) { }
    
    T rand() noexcept {
        if (counter >= clusterSize || current < min || current > max) {
            current = dis(re);
            counter = 0;
        }
        ++counter;
        return current++;
    }

private:
    std::mt19937_64 re;
    std::uniform_int_distribution<T> dis;
    const T min;
    const T max;
    const std::size_t clusterSize;
    T current;
    std::size_t counter = clusterSize;
};

template<class T>
std::vector<T> generateClusterIntVector(const std::size_t size, const T seed = 42, 
        const T min = std::numeric_limits<T>::min(),
        const T max = std::numeric_limits<T>::max(), const std::size_t clusterSize = 1000) {
    ClusterIntGenerator<T> gen(seed, min, max, clusterSize);
    return generateRandomVector<T>(size, gen);;
}

// Generates a vector with size random unique elements. May run indefinitely!
template<class T>
std::vector<T> generateUniqueClusterIntVector(const std::size_t size, const T seed = 42, 
        const T min = std::numeric_limits<T>::min(),
        const T max = std::numeric_limits<T>::max(), const std::size_t clusterSize = 1000) {
    ClusterIntGenerator<T> gen(seed, min, max, clusterSize);
    return generateUniqueRandomVector<T>(size, gen);
}

template<class T>
class NormalIntGenerator {
    static_assert(sizeof(T) <= 4, "T must have at most 32 bits");
public:
    NormalIntGenerator(const T seed, const T mean, const T std)
        : re(seed), dis((double)mean, (double)std) { }
    
    T rand() noexcept {
        return static_cast<T>(std::lround(dis(re)));
    }

private:
    std::mt19937_64 re;
    std::normal_distribution<double> dis;
};

template<class T>
std::vector<T> generateNormalIntVector(const std::size_t size, const T seed, 
        const T mean, const T std) {
    NormalIntGenerator<T> gen(seed, mean, std);
    return generateRandomVector<T>(size, gen);
}

// Generates a vector with size random normal distributed unique elements. May run indefinitely!
template<class T>
std::vector<T> generateUniqueNormalIntVector(const std::size_t size, const T seed, 
        const T mean, const T std) {
    NormalIntGenerator<T> gen(seed, mean, std);
    return generateUniqueRandomVector<T>(size, gen);
}

template<class T, bool increasing = true>
class LinearProbabilityIntGenerator {
    static_assert(sizeof(T) <= 4, "T must have at most 32 bits");

    static constexpr std::array<double, 2> weightsIncreasing{0., 1.};
    static constexpr std::array<double, 2> weightsDecreasing{1., 0.};
public:
    LinearProbabilityIntGenerator(const T seed = 42, const T min = std::numeric_limits<T>::min(),
            const T max = std::numeric_limits<T>::max()) : re(seed) {
        const std::array<double, 2> intervals{static_cast<double>(min), static_cast<double>(max)};
        auto weights = increasing ? weightsIncreasing.begin() : weightsDecreasing.begin();
        dis = std::piecewise_linear_distribution<double>(intervals.begin(), intervals.end(), weights);
    }
    
    T rand() noexcept {
        return dis(re);
    }

private:
    std::mt19937_64 re;
    std::piecewise_linear_distribution<double> dis;
};

template<class T, bool increasing = true>
std::vector<T> generateLinearProbabilityIntVector(const std::size_t size, const T seed = 42, 
        const T min = std::numeric_limits<T>::min(),
        const T max = std::numeric_limits<T>::max()) {
    LinearProbabilityIntGenerator<T, increasing> gen(seed, min, max);
    return generateRandomVector<T>(size, gen);
}

// Generates a vector with size random unique elements. May run indefinitely!
template<class T, bool increasing = true>
std::vector<T> generateUniqueLinearProbabilityIntVector(const std::size_t size, const T seed = 42, 
        const T min = std::numeric_limits<T>::min(),
        const T max = std::numeric_limits<T>::max()) {
    LinearProbabilityIntGenerator<T, increasing> gen(seed, min, max);
    return generateUniqueRandomVector<T>(size, gen);
}

template<class T>
class UniformRealGenerator {
public:
    UniformRealGenerator(const T seed = 42, const T min = std::numeric_limits<T>::min(),
                        const T max = std::numeric_limits<T>::max())
        : re(seed), dis(min, max) { }
    
    T rand() noexcept {
        return dis(re);
    }

private:
    std::mt19937_64 re;
    std::uniform_real_distribution<T> dis;
};

template<class T>
std::vector<T> generateUniformRealVector(const std::size_t size, const T seed = 42, 
        const T min = std::numeric_limits<T>::min(),
        const T max = std::numeric_limits<T>::max()) {
    UniformRealGenerator<T> gen(seed, min, max);
    return generateRandomVector<T>(size, gen);;
}
