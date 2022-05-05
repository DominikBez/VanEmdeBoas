// Copyright Dominik Bez 2022.

#include <fstream>
#include <chrono>
#include <vector>
#include <queue>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <utility>
#include <string>
#include <type_traits>
#include <omp.h>

#include "utils/csv_printer.hpp"
#include "utils/generators.hpp"
#include "utils/utils.hpp"
#include "implementation/van_emde_boas_locked.hpp"
#include "implementation/van_emde_boas32_locked.hpp"
#include "implementation/van_emde_boas32_locked_top.hpp"
#include "implementation/van_emde_boas32_locked_fine_grained.hpp"
#include "implementation/van_emde_boas32_lockless.hpp"

using std::size_t, std::uint8_t, std::int8_t, std::uint32_t, std::int32_t,
    std::uint64_t, std::int64_t;

template<int K, class T = typename veb::VanEmdeBoasLocked<K>::TypeT>
using VEBL = veb::VanEmdeBoasLocked<K, T>;
template<class T = typename veb::VanEmdeBoas32Locked<>::TypeT>
using VEB32L = veb::VanEmdeBoas32Locked<T>;
template<class T = typename veb::VanEmdeBoas32LockedTop<>::TypeT>
using VEB32LT = veb::VanEmdeBoas32LockedTop<T>;
template<class T = typename veb::VanEmdeBoas32LockedFineGrained<>::TypeT>
using VEB32LFG = veb::VanEmdeBoas32LockedFineGrained<T>;
template<class T = typename veb::VanEmdeBoas32Lockless<>::TypeT>
using VEB32LL = veb::VanEmdeBoas32Lockless<T>;

constexpr int ITERATIONS = 10;
constexpr size_t MIN_SIZE = 1LLU << 5;
constexpr size_t MAX_SIZE = 1LLU << 21;
constexpr int LOOKUPS = 10000;
constexpr int MAX_THREADS = 8;

// Returns per insert/lookup/remove timings
template<class Set, class T>
std::tuple<double, double, double> timing(const size_t size,
                                          const std::vector<T>& values,
                                          const std::vector<T>& lookups) {
    Set tree;

    auto insert0 = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for (size_t i = 0; i < size; ++i) {
        tree.insert(values[i]);
    }

    auto insert1 = std::chrono::high_resolution_clock::now();
    auto lookup0 = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for (const T e : lookups) {
        auto result = tree.lower_bound(e);
        (void)result;
    }

    auto lookup1 = std::chrono::high_resolution_clock::now();
    auto remove0 = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for (size_t i = 0; i < size; ++i) {
        tree.erase(values[i]);
    }

    auto remove1 = std::chrono::high_resolution_clock::now();
    return {std::chrono::duration_cast<std::chrono::nanoseconds>(insert1 - insert0).count() / (double)size,
        std::chrono::duration_cast<std::chrono::nanoseconds>(lookup1 - lookup0).count() / (double)lookups.size(),
        std::chrono::duration_cast<std::chrono::nanoseconds>(remove1 - remove0).count() / (double)size};
}

template<int K, class Set, class T, class Gen>
void runTiming(std::ofstream& out, const std::string& name, const std::string& distribution, Gen& gen) {
    constexpr T SEED = 42;
    constexpr T MIN = std::numeric_limits<T>::min();
    constexpr T MAX = std::numeric_limits<T>::max();
    constexpr T MAX_VALUE = (K == 8 * sizeof(T))
#ifdef __SIZEOF_INT128__
        ? MAX - 1 : T((__uint128_t(1) << K) - 1);
    constexpr size_t NUM_VALUES = std::min((__uint128_t)MAX_SIZE,
        (__uint128_t)(std::is_signed_v<T> ? MAX_VALUE : MAX_VALUE / 2));
#else
        ? MAX - 1 : T((1LLU << K) - 1);
    constexpr size_t NUM_VALUES = std::min(MAX_SIZE,
        (size_t)(std::is_signed_v<T> ? MAX_VALUE : MAX_VALUE / 2));
#endif
    for (int threads = 1; threads <= MAX_THREADS; ++threads) {
        omp_set_num_threads(threads);
        for (int i = 0; i < ITERATIONS; ++i) {
            auto values = generateUniqueRandomVector<T>(NUM_VALUES, gen);
            utils::shuffle(values, SEED - 1);
            const auto lookups = generateUniformIntVector<T>(LOOKUPS, SEED + i, MIN, MAX_VALUE);
            for (size_t size = MIN_SIZE; size <= NUM_VALUES; size *= 2) {
                const auto [insert, lookup, remove] = timing<Set, T>(size, values, lookups);
                csvPrint(out, name, K, distribution, threads, i, size, insert, lookup, remove);
            }
        }
    }
}

template<int K, class Set, class T>
void run(std::ofstream& out, const std::string& name) {
    constexpr T SEED = 11;
    constexpr T MIN = std::numeric_limits<T>::min();
    constexpr T MAX = std::numeric_limits<T>::max();
    constexpr T MAX_VALUE = (K == 8 * sizeof(T))
#ifdef __SIZEOF_INT128__
        ? MAX - 1 : T((__uint128_t(1) << K) - 1);
#else
        ? MAX - 1 : T((1LLU << K) - 1);
#endif

    UniformIntGenerator<T> uniformGen(SEED, MIN, MAX_VALUE);
    runTiming<K, Set, T>(out, name, "uniform", uniformGen);
    constexpr size_t CLUSTER_SIZE = 1000;
    ClusterIntGenerator<T> clusterGen(SEED, MIN, MAX_VALUE, CLUSTER_SIZE);
    runTiming<K, Set, T>(out, name, "cluster", clusterGen);
    if constexpr (sizeof(T) == 4 && K == 32) {
        constexpr T MEAN = (MIN + MAX_VALUE) / 2;
        NormalIntGenerator<T> normalGen(SEED, MEAN, (MAX_VALUE - MEAN) / 10);
        runTiming<K, Set, T>(out, name, "normal", normalGen);
        LinearProbabilityIntGenerator<T, true> incProbGen(SEED, MIN, MAX_VALUE);
        runTiming<K, Set, T>(out, name, "incProb", incProbGen);
        LinearProbabilityIntGenerator<T, false> decProbGen(SEED, MIN, MAX_VALUE);
        runTiming<K, Set, T>(out, name, "decProb", decProbGen);
    }
}

int main() {
    std::ofstream out("parallelVeb.csv");
    csvPrint(out, "tree", "K", "distribution", "threads", "iteration", "size", "insert", "lookup", "remove");

    run<32, VEBL<32>, uint32_t>(out, "32uVEBL");
    run<32, VEB32L<>, uint32_t>(out, "uVEB32L");
    run<32, VEB32LT<>, uint32_t>(out, "uVEB32LT");
    run<32, VEB32LFG<>, uint32_t>(out, "uVEB32LFG");
    run<32, VEB32LL<>, uint32_t>(out, "uVEB32LL");
}
