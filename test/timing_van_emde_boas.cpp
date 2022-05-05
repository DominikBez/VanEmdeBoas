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
#include <set>

#include "utils/csv_printer.hpp"
#include "utils/generators.hpp"
#include "utils/utils.hpp"
#include "implementation/van_emde_boas.hpp"
#include "implementation/van_emde_boas_locked.hpp"
#include "implementation/van_emde_boas32.hpp"
#include "implementation/van_emde_boas32_locked.hpp"
#include "implementation/van_emde_boas32_locked_top.hpp"
#include "implementation/van_emde_boas32_locked_fine_grained.hpp"
#include "implementation/van_emde_boas32_lockless.hpp"

using std::size_t, std::uint8_t, std::int8_t, std::uint32_t, std::int32_t,
    std::uint64_t, std::int64_t;

template<int K, class T = typename veb::VanEmdeBoas<K>::TypeT>
using VEB = veb::VanEmdeBoas<K, T>;
template<int K, class T = typename veb::VanEmdeBoasLocked<K>::TypeT>
using VEBL = veb::VanEmdeBoasLocked<K, T>;
template<class T = typename veb::VanEmdeBoas32<>::TypeT>
using VEB32 = veb::VanEmdeBoas32<T>;
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
constexpr size_t MAX_SIZE = 1LLU << 23;
constexpr int LOOKUPS = 10000;

// Returns per insert/lookup/remove timings
template<class Set, class T>
std::tuple<double, double, double> timing(const size_t size,
                                          const std::vector<T>& values,
                                          const std::vector<T>& lookups) {
    Set tree;

    auto insert0 = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < size; ++i) {
        tree.insert(values[i]);
    }

    auto insert1 = std::chrono::high_resolution_clock::now();
    auto lookup0 = std::chrono::high_resolution_clock::now();

    for (const T e : lookups) {
        auto result = tree.lower_bound(e);
        (void)result;
    }

    auto lookup1 = std::chrono::high_resolution_clock::now();
    auto remove0 = std::chrono::high_resolution_clock::now();

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
    for (int i = 0; i < ITERATIONS; ++i) {
        const auto values = generateUniqueRandomVector<T>(NUM_VALUES, gen);
        const auto lookups = generateUniformIntVector<T>(LOOKUPS, SEED + i, MIN, MAX_VALUE);
        for (size_t size = MIN_SIZE; size <= NUM_VALUES; size *= 2) {
            const auto [insert, lookup, remove] = timing<Set, T>(size, values, lookups);
            csvPrint(out, name, K, distribution, i, size, insert, lookup, remove);
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

template<int K>
void runUnsignedVEB(std::ofstream& out) {
    using UVEB = VEB<K>;
    run<K, UVEB, typename UVEB::TypeT>(out, std::to_string(K) + "uVEB");
}

template<int K>
void runSignedVEB(std::ofstream& out) {
    using UVEB = VEB<K>;
    using T = typename std::make_signed_t<typename UVEB::TypeT>;
    run<K, VEB<K, T>, T>(out, std::to_string(K) + "VEB");
}

int main() {
    std::ofstream out("veb.csv");
    csvPrint(out, "tree", "K", "distribution", "iteration", "size", "insert", "lookup", "remove");
    run<8, std::set<uint8_t>, uint8_t>(out, "8std::set");
    run<16, std::set<uint16_t>, uint16_t>(out, "16std::set");
    run<32, std::set<uint32_t>, uint32_t>(out, "32std::set");
    run<64, std::set<uint64_t>, uint64_t>(out, "64std::set");
    runUnsignedVEB<6>(out);
    runUnsignedVEB<8>(out);
    runSignedVEB<8>(out);
    runUnsignedVEB<12>(out);
    runUnsignedVEB<13>(out);
    runUnsignedVEB<14>(out);
    runUnsignedVEB<16>(out);
    runSignedVEB<16>(out);
    runUnsignedVEB<32>(out);
    runSignedVEB<32>(out);
    runUnsignedVEB<64>(out);
    runSignedVEB<64>(out);
#ifdef __SIZEOF_INT128__
    runUnsignedVEB<127>(out);
    runUnsignedVEB<128>(out);
    runSignedVEB<128>(out);
#endif

    run<32, VEBL<32>, uint32_t>(out, "32uVEBL");

    run<32, VEB32<>, uint32_t>(out, "uVEB32");
    run<32, VEB32<int32_t>, int32_t>(out, "VEB32");

    run<32, VEB32L<>, uint32_t>(out, "uVEB32L");
    run<32, VEB32LT<>, uint32_t>(out, "uVEB32LT");
    run<32, VEB32LFG<>, uint32_t>(out, "uVEB32LFG");
    run<32, VEB32LL<>, uint32_t>(out, "uVEB32LL");
}
