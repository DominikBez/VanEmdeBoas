// Copyright Dominik Bez 2022.

#include <limits>
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <type_traits>
#include <unordered_set>

#include <gtest/gtest.h>

#include "implementation/van_emde_boas.hpp"
#include "implementation/van_emde_boas_locked.hpp"
#include "implementation/van_emde_boas32.hpp"
#include "implementation/van_emde_boas32_locked.hpp"
#include "implementation/van_emde_boas32_locked_top.hpp"
#include "implementation/van_emde_boas32_locked_fine_grained.hpp"
#include "implementation/van_emde_boas32_lockless.hpp"
#include "utils/generators.hpp"
#include "utils/utils.hpp"

using std::uint32_t, std::uint64_t, std::int8_t, std::int16_t, std::int32_t, std::int64_t;

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

#define EXPECT_INVALID(x) \
    if constexpr (floatType) { \
        const auto result = (x); \
        EXPECT_TRUE(std::isnan(result)) << "Result is " << result; \
    } else { \
        EXPECT_EQ((x), MAX); \
    }

template<class Tree>
class SmallTest : public ::testing::Test {
protected:
    static constexpr int K = Tree::BitsK;
    using T = typename Tree::TypeT;
    static_assert(K >= 5, "SmallTest assumes K >= 5");
};
using SmallTestTypes = ::testing::Types<VEB<5>, VEB<5, uint32_t>,
    VEB<6>, VEB<6, uint64_t>, VEB<8>, VEB<8, int8_t>, VEB<9>, VEB<15>,
    VEB<16>, VEB<16, int16_t>, VEB<31>, VEB<32>, VEB<32, int32_t>, VEB<63>,
    VEB<64>, VEB<64, int64_t>, VEBL<64>,
#ifdef __SIZEOF_INT128__
    VEB<100>, VEB<127>, VEB<128>, VEB<128, __int128_t>,
    //typename std::conditional<(sizeof(long double) == 16), VEB<128, long double>, void>::type,
    // long double probably doesn't work
#endif
    VEB<32, float>, VEB<64, double>, VEB32<>, VEB32<int32_t>, VEB32<float>,
    VEB32LT<>, VEB32LT<int32_t>, VEB32LT<float>,
    VEB32L<>, VEB32L<int32_t>, VEB32L<float>,
    VEB32LFG<>, VEB32LFG<int32_t>, VEB32LFG<float>,
    VEB32LL<>, VEB32LL<int32_t>, VEB32LL<float>
    >;
TYPED_TEST_SUITE(SmallTest, SmallTestTypes);

TYPED_TEST(SmallTest, ArbitraryInsert) {
    if constexpr (std::is_same_v<TypeParam, void>) return;
    using T = typename TestFixture::T;
    constexpr bool floatType = std::is_floating_point_v<T>;
    [[maybe_unused]] constexpr auto MAX = std::numeric_limits<T>::max();
    TypeParam tree;
    for (T i = 0; i < 32; ++i)
        EXPECT_INVALID(tree.locate(i));
    
    tree.insert(17);
    for (T i = 0; i <= 17; ++i)
        EXPECT_EQ(tree.locate(i), 17);
    for (T i = 18; i < 32; ++i)
        EXPECT_INVALID(tree.locate(i));

    tree.insert(30);
    for (T i = 0; i <= 17; ++i)
        EXPECT_EQ(tree.locate(i), 17);
    for (T i = 18; i <= 30; ++i)
        EXPECT_EQ(tree.locate(i), 30);
    for (T i = 31; i < 32; ++i)
        EXPECT_INVALID(tree.locate(i));
    
    tree.insert(3);
    for (T i = 0; i <= 3; ++i)
        EXPECT_EQ(tree.locate(i), 3);
    for (T i = 4; i <= 17; ++i)
        EXPECT_EQ(tree.locate(i), 17);
    for (T i = 18; i <= 30; ++i)
        EXPECT_EQ(tree.locate(i), 30);
    for (T i = 31; i < 32; ++i)
        EXPECT_INVALID(tree.locate(i));
    
    tree.insert(31);
    for (T i = 0; i <= 3; ++i)
        EXPECT_EQ(tree.locate(i), 3);
    for (T i = 4; i <= 17; ++i)
        EXPECT_EQ(tree.locate(i), 17);
    for (T i = 18; i <= 30; ++i)
        EXPECT_EQ(tree.locate(i), 30);
    for (T i = 31; i < 32; ++i)
        EXPECT_EQ(tree.locate(i), 31);
    
    tree.insert(0);
    EXPECT_EQ(tree.locate(0), 0);
    for (T i = 1; i <= 3; ++i)
        EXPECT_EQ(tree.locate(i), 3);
    for (T i = 4; i <= 17; ++i)
        EXPECT_EQ(tree.locate(i), 17);
    for (T i = 18; i <= 30; ++i)
        EXPECT_EQ(tree.locate(i), 30);
    for (T i = 31; i < 32; ++i)
        EXPECT_EQ(tree.locate(i), 31);
    
    tree.remove(0);
    for (T i = 0; i <= 3; ++i)
        EXPECT_EQ(tree.locate(i), 3);
    for (T i = 4; i <= 17; ++i)
        EXPECT_EQ(tree.locate(i), 17);
    for (T i = 18; i <= 30; ++i)
        EXPECT_EQ(tree.locate(i), 30);
    for (T i = 31; i < 32; ++i)
        EXPECT_EQ(tree.locate(i), 31);
    
    tree.remove(17);
    for (T i = 0; i <= 3; ++i)
        EXPECT_EQ(tree.locate(i), 3);
    for (T i = 4; i <= 30; ++i)
        EXPECT_EQ(tree.locate(i), 30);
    for (T i = 31; i < 32; ++i)
        EXPECT_EQ(tree.locate(i), 31);
    
    tree.remove(3);
    for (T i = 0; i <= 30; ++i)
        EXPECT_EQ(tree.locate(i), 30);
    for (T i = 31; i < 32; ++i)
        EXPECT_EQ(tree.locate(i), 31);
    
    tree.remove(31);
    for (T i = 0; i <= 30; ++i)
        EXPECT_EQ(tree.locate(i), 30);
    for (T i = 31; i < 32; ++i)
        EXPECT_INVALID(tree.locate(i));
    
    tree.remove(30);
    for (T i = 0; i < 32; ++i)
        EXPECT_INVALID(tree.locate(i));
}

TYPED_TEST(SmallTest, ForwardInsert) {
    if constexpr (std::is_same_v<TypeParam, void>) return;
    using T = typename TestFixture::T;
    constexpr bool floatType = std::is_floating_point_v<T>;
    [[maybe_unused]] constexpr auto MAX = std::numeric_limits<T>::max();
    TypeParam tree;
    for (T i = 0; i < 32; ++i) {
        tree.insert(i);
        for (T k = 0; k <= i; ++k)
            EXPECT_EQ(tree.locate(k), k);
        for (T k = i + 1; k < 32; ++k)
            EXPECT_INVALID(tree.locate(k));
    }

    for (T i = 0; i < 31; ++i) {
        tree.remove(i);
        for (T k = 0; k <= i; ++k)
            EXPECT_EQ(tree.locate(k), i + 1);
        for (T k = i + 1; k < 32; ++k)
            EXPECT_EQ(tree.locate(k), k);
    }
    tree.remove(31);
    for (T i = 0; i < 32; ++i)
        EXPECT_INVALID(tree.locate(i));
}

TYPED_TEST(SmallTest, BackwardsInsert) {
    if constexpr (std::is_same_v<TypeParam, void>) return;
    using T = typename TestFixture::T;
    constexpr bool floatType = std::is_floating_point_v<T>;
    [[maybe_unused]] constexpr auto MAX = std::numeric_limits<T>::max();
    TypeParam tree;
    for (T j = 32; j > 0; --j) {
        T i = j - 1;
        tree.insert(i);
        for (T k = 0; k <= i; ++k)
            EXPECT_EQ(tree.locate(k), i);
        for (T k = j; k < 32; ++k)
            EXPECT_EQ(tree.locate(k), k);
    }

    for (T j = 32; j > 0; --j) {
        T i = j - 1;
        tree.remove(i);
        for (T k = 0; k < i; ++k)
            EXPECT_EQ(tree.locate(k), k);
        for (T k = i; k < 32; ++k)
            EXPECT_INVALID(tree.locate(k));
    }
}

template<class Tree>
class RandomTest : public ::testing::Test {
protected:
    static constexpr int K = Tree::BitsK;
    using T = typename Tree::TypeT;
};
using RandomTestTypes = ::testing::Types<VEB<1>, VEB<2>, VEB<3>,
    VEB<4>, VEB<5>, VEB<5, uint32_t>,
    VEB<6>, VEB<6, uint64_t>, VEB<8>, VEB<8, int8_t>, VEB<9>, VEB<15>,
    VEB<16>, VEB<16, int16_t>, VEB<31>, VEB<32>, VEB<32, int32_t>, VEB<63>,
    VEB<64>, VEB<64, int64_t>, VEBL<64>,
#ifdef __SIZEOF_INT128__
    VEB<100>, VEB<127>, VEB<128>, VEB<128, __int128_t>,
#endif
    VEB32<>, VEB32<int32_t>,
    VEB32LT<>, VEB32LT<int32_t>,
    VEB32L<>, VEB32L<int32_t>,
    VEB32LFG<>, VEB32LFG<int32_t>,
    VEB32LL<>, VEB32LL<int32_t>
    >;
TYPED_TEST_SUITE(RandomTest, RandomTestTypes);

TYPED_TEST(RandomTest, UniformCompareWithStdSet) {
    constexpr std::size_t N = 100'000;
    constexpr int TESTS = 100;
    using T = typename TestFixture::T;
    constexpr T SEED = 2;
    constexpr T MIN = std::numeric_limits<T>::min();
    constexpr T MAX = std::numeric_limits<T>::max();
    constexpr T MAX_VALUE = (TestFixture::K == 8 * sizeof(T))
#ifdef __SIZEOF_INT128__
        ? MAX - 1 : T((__uint128_t(1) << TestFixture::K) - 1);
#else
        ? MAX - 1 : T((1LLU << TestFixture::K) - 1);
#endif
    TypeParam tree;
    std::set<T> set;
    UniformIntGenerator<T> gen(SEED, MIN, MAX_VALUE);

    for (std::size_t i = 0; i < N; ++i) {
        const T e = gen.rand();
        tree.insert(e);
        set.insert(e);

        for (int j = 0; j < TESTS; ++j) {
            const T test = gen.rand();
            const auto setIt = set.lower_bound(test);
            const T setResult = (setIt == set.end()) ? MAX : *setIt;
            EXPECT_EQ(tree.locate(test), setResult);
        }
    }

    std::vector<T> elements(set.begin(), set.end());
    utils::shuffle(elements, SEED + 1);

    for (const T e : elements) {
        EXPECT_EQ(tree.locate(e), e);
    }

    for (const T e : elements) {
        tree.remove(e);
        set.erase(e);

        for (int i = 0; i < TESTS; ++i) {
            const T test = gen.rand();
            const auto setIt = set.lower_bound(test);
            const T setResult = (setIt == set.end()) ? MAX : *setIt;
            EXPECT_EQ(tree.locate(test), setResult);
        }
    }
}

TYPED_TEST(RandomTest, ClusterCompareWithStdSet) {
    constexpr std::size_t N = 100'000;
    constexpr int TESTS = 100;
    using T = typename TestFixture::T;
    constexpr T SEED = 5;
    constexpr T MIN = std::numeric_limits<T>::min();
    constexpr T MAX = std::numeric_limits<T>::max();
    constexpr T MAX_VALUE = (TestFixture::K == 8 * sizeof(T))
#ifdef __SIZEOF_INT128__
        ? MAX - 1 : T((__uint128_t(1) << TestFixture::K) - 1);
#else
        ? MAX - 1 : T((1LLU << TestFixture::K) - 1);
#endif
    constexpr std::size_t CLUSTER_SIZE = 100;
    TypeParam tree;
    std::set<T> set;
    ClusterIntGenerator<T> gen(SEED, MIN, MAX_VALUE, CLUSTER_SIZE);
    UniformIntGenerator<T> testGen(SEED + 1, MIN, MAX_VALUE);

    for (std::size_t i = 0; i < N; ++i) {
        const T e = gen.rand();
        tree.insert(e);
        set.insert(e);

        for (int j = 0; j < TESTS; ++j) {
            const T test = testGen.rand();
            const auto setIt = set.lower_bound(test);
            const T setResult = (setIt == set.end()) ? MAX : *setIt;
            EXPECT_EQ(tree.locate(test), setResult);
        }
    }

    std::vector<T> elements(set.begin(), set.end());
    utils::shuffle(elements, SEED + 2);

    for (const T e : elements) {
        EXPECT_EQ(tree.locate(e), e);
    }

    for (const T e : elements) {
        tree.remove(e);
        set.erase(e);

        for (int i = 0; i < TESTS; ++i) {
            const T test = testGen.rand();
            const auto setIt = set.lower_bound(test);
            const T setResult = (setIt == set.end()) ? MAX : *setIt;
            EXPECT_EQ(tree.locate(test), setResult);
        }
    }
}

TYPED_TEST(RandomTest, BigClusterCompareWithStdSet) {
    constexpr std::size_t N = 200'000;
    constexpr int TESTS = 50;
    using T = typename TestFixture::T;
    constexpr T SEED = 8;
    constexpr T MIN = std::numeric_limits<T>::min();
    constexpr T MAX = std::numeric_limits<T>::max();
    constexpr T MAX_VALUE = (TestFixture::K == 8 * sizeof(T))
#ifdef __SIZEOF_INT128__
        ? MAX - 1 : T((__uint128_t(1) << TestFixture::K) - 1);
#else
        ? MAX - 1 : T((1LLU << TestFixture::K) - 1);
#endif
    constexpr std::size_t CLUSTER_SIZE = N / 2;
    TypeParam tree;
    std::set<T> set;
    std::vector<T> elements = generateClusterIntVector(N, SEED, MIN, MAX_VALUE, CLUSTER_SIZE);
    utils::shuffle(elements, SEED + 1);
    UniformIntGenerator<T> testGen(SEED + 2, MIN, MAX_VALUE);

    for (const T e : elements) {
        tree.insert(e);
        set.insert(e);

        for (int j = 0; j < TESTS; ++j) {
            const T test = testGen.rand();
            const auto setIt = set.lower_bound(test);
            const T setResult = (setIt == set.end()) ? MAX : *setIt;
            EXPECT_EQ(tree.locate(test), setResult);
        }
    }

    elements.resize(set.size());
    elements.assign(set.begin(), set.end());
    utils::shuffle(elements, SEED + 3);

    for (const T e : elements) {
        EXPECT_EQ(tree.locate(e), e);
    }

    for (const T e : elements) {
        tree.remove(e);
        set.erase(e);

        for (int i = 0; i < TESTS; ++i) {
            const T test = testGen.rand();
            const auto setIt = set.lower_bound(test);
            const T setResult = (setIt == set.end()) ? MAX : *setIt;
            EXPECT_EQ(tree.locate(test), setResult);
        }
    }
}

template<class Tree>
class FloatRandomTest : public ::testing::Test {
protected:
    static constexpr int K = Tree::BitsK;
    using T = typename Tree::TypeT;
};
using FloatRandomTestTypes = ::testing::Types<VEB<32, float>, VEB<64, double>, VEBL<64, double>,
#ifdef __SIZEOF_INT128__
    // typename std::conditional<(sizeof(long double) == 16), VEB<128, long double>, void>::type,
    // long double probably doesn't work
#endif
    VEB32<float>,
    VEB32L<float>,
    VEB32LT<float>,
    VEB32LFG<float>,
    VEB32LL<float>
    >;
TYPED_TEST_SUITE(FloatRandomTest, FloatRandomTestTypes);

TYPED_TEST(FloatRandomTest, UniformCompareWithStdSet) {
    if constexpr (std::is_same_v<TypeParam, void>) return;
    constexpr std::size_t N = 100'000;
    constexpr int TESTS = 100;
    using T = typename TestFixture::T;
    constexpr bool floatType = true;
    constexpr T SEED = 11;
    constexpr T MIN = std::numeric_limits<T>::min();
    constexpr T MAX = std::numeric_limits<T>::max();
    TypeParam tree;
    std::set<T> set;
    UniformRealGenerator<T> gen(SEED, MIN, MAX);

    for (std::size_t i = 0; i < N; ++i) {
        const T e = gen.rand();
        tree.insert(e);
        set.insert(e);

        for (int j = 0; j < TESTS; ++j) {
            const T test = gen.rand();
            const auto setIt = set.lower_bound(test);
            if (setIt == set.end()) {
                EXPECT_INVALID(tree.locate(test));
            } else {
                EXPECT_EQ(tree.locate(test), *setIt);
            }
        }
    }

    std::vector<T> elements(set.begin(), set.end());
    utils::shuffle(elements, SEED + 1);

    for (const T e : elements) {
        EXPECT_EQ(tree.locate(e), e);
    }

    for (const T e : elements) {
        tree.remove(e);
        set.erase(e);

        for (int i = 0; i < TESTS; ++i) {
            const T test = gen.rand();
            const auto setIt = set.lower_bound(test);
            if (setIt == set.end()) {
                EXPECT_INVALID(tree.locate(test));
            } else {
                EXPECT_EQ(tree.locate(test), *setIt);
            }
        }
    }
}
