// Copyright Dominik Bez 2022.

#include <limits>
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <type_traits>
#include <thread>
#include <unordered_set>
#include <mutex>

#include <gtest/gtest.h>

#include "implementation/van_emde_boas_locked.hpp"
#include "implementation/van_emde_boas32_locked.hpp"
#include "implementation/van_emde_boas32_locked_top.hpp"
#include "implementation/van_emde_boas32_locked_fine_grained.hpp"
#include "implementation/van_emde_boas32_lockless.hpp"
#include "utils/generators.hpp"
#include "utils/utils.hpp"

using std::uint32_t, std::uint64_t, std::int8_t, std::int16_t, std::int32_t, std::int64_t;

template<class T = typename veb::VanEmdeBoasLocked<32>::TypeT>
using VEBL = veb::VanEmdeBoasLocked<32, T>;
template<class T = typename veb::VanEmdeBoas32Locked<>::TypeT>
using VEB32L = veb::VanEmdeBoas32Locked<T>;
template<class T = typename veb::VanEmdeBoas32LockedTop<>::TypeT>
using VEB32LT = veb::VanEmdeBoas32LockedTop<T>;
template<class T = typename veb::VanEmdeBoas32LockedFineGrained<>::TypeT>
using VEB32LFG = veb::VanEmdeBoas32LockedFineGrained<T>;
template<class T = typename veb::VanEmdeBoas32Lockless<>::TypeT>
using VEB32LL = veb::VanEmdeBoas32Lockless<T>;

template<class Tree>
class ParallelTest : public ::testing::Test {
protected:
    static constexpr int K = Tree::BitsK;
    using T = typename Tree::TypeT;
};
using ParallelTestTypes = ::testing::Types<
    VEBL<>, VEBL<int32_t>,
    VEB32L<>, VEB32L<int32_t>,
    VEB32LT<>, VEB32LT<int32_t>,
    VEB32LFG<>, VEB32LFG<int32_t>,
    VEB32LL<>, VEB32LL<int32_t>
    >;
TYPED_TEST_SUITE(ParallelTest, ParallelTestTypes);

TYPED_TEST(ParallelTest, SingleProducerSingleConsumer) {
    constexpr std::size_t N = 10'000'000;
    constexpr int TESTS = 8;
    using T = typename TestFixture::T;
    using UT = typename std::make_unsigned_t<T>;
    constexpr T SEED = 20;
    constexpr T MIN = std::numeric_limits<T>::min();
    constexpr T MAX = std::numeric_limits<T>::max();
    constexpr T POISON_PILL = MAX - 1;
    constexpr T MAX_VALUE = MAX - 2;
    TypeParam tree;
    std::vector<T> values = generateUniformIntVector(N, SEED, MIN, MAX_VALUE);
    std::unordered_multiset<T> set(values.begin(), values.end());
    std::unordered_multiset<T> removed;

    std::thread producer([&]{
        for (const T e : values) {
            tree.insert(e);
        }
        tree.insert(POISON_PILL);
    });

    std::thread consumer([&]{
        UniformIntGenerator<T> gen(SEED + 1, MIN, MAX_VALUE);
        T results[TESTS];
        bool poisoned = false;
        while (true) {
            for (int i = 0; i < TESTS; ++i) {
                const T test = gen.rand();
                const T result = tree.locate(test);
                if (result == POISON_PILL) {
                    poisoned = true;
                    break;
                } else if (result != MAX) {
                    EXPECT_NE(set.find(result), set.end());
                    EXPECT_GE(result, test);
                }
                results[i] = result;
            }
            if (poisoned)
                break;
            const T removing = results[(UT)gen.rand() % TESTS];
            if (removing != MAX) {
                removed.insert(removing);
                EXPECT_LE(removed.count(removing), set.count(removing));
                tree.remove(removing);
                EXPECT_GE(tree.locate(removing), removing); // producer might have added it again => not _GT
            }
        }
        // Producer is done
        EXPECT_EQ(tree.locate(POISON_PILL), POISON_PILL);
        tree.remove(POISON_PILL);
        EXPECT_EQ(tree.locate(POISON_PILL), MAX);
        T removing = MIN;
        while ((removing = tree.locate(removing)) != MAX) {
            removed.insert(removing);
            EXPECT_LE(removed.count(removing), set.count(removing));
            EXPECT_EQ(tree.locate(removing), removing);
            tree.remove(removing);
            EXPECT_GT(tree.locate(removing), removing);
        }

        for (const T e : values) {
            EXPECT_GT(removed.count(e), 0) << "Lost element " << e;
        }
    });

    producer.join();
    consumer.join();
}

TYPED_TEST(ParallelTest, UniqueSingleProducerSingleConsumer) {
    constexpr std::size_t N = 10'000'000;
    constexpr int TESTS = 8;
    using T = typename TestFixture::T;
    using UT = typename std::make_unsigned_t<T>;
    constexpr T SEED = 25;
    constexpr T MIN = std::numeric_limits<T>::min();
    constexpr T MAX = std::numeric_limits<T>::max();
    constexpr T POISON_PILL = MAX - 1;
    constexpr T MAX_VALUE = MAX - 2;
    TypeParam tree;
    std::vector<T> values = generateUniqueUniformIntVector(N, SEED, MIN, MAX_VALUE);
    std::unordered_set<T> set(values.begin(), values.end());
    std::unordered_set<T> removed;

    std::thread producer([&]{
        for (const T e : values) {
            tree.insert(e);
        }
        tree.insert(POISON_PILL);
    });

    std::thread consumer([&]{
        UniformIntGenerator<T> gen(SEED + 1, MIN, MAX_VALUE);
        T results[TESTS];
        bool poisoned = false;
        while (true) {
            for (int i = 0; i < TESTS; ++i) {
                const T test = gen.rand();
                const T result = tree.locate(test);
                if (result == POISON_PILL) {
                    poisoned = true;
                    break;
                } else if (result != MAX) {
                    EXPECT_NE(set.find(result), set.end());
                    EXPECT_GE(result, test);
                }
                results[i] = result;
            }
            if (poisoned)
                break;
            const T removing = results[(UT)gen.rand() % TESTS];
            if (removing != MAX) {
                EXPECT_EQ(removed.count(removing), 0);
                removed.insert(removing);
                tree.remove(removing);
                EXPECT_GT(tree.locate(removing), removing);
            }
        }
        // Producer is done
        EXPECT_EQ(tree.locate(POISON_PILL), POISON_PILL);
        tree.remove(POISON_PILL);
        EXPECT_EQ(tree.locate(POISON_PILL), MAX);
        T removing = MIN;
        while ((removing = tree.locate(removing)) != MAX) {
            EXPECT_EQ(tree.locate(removing), removing);
            removed.insert(removing);
            tree.remove(removing);
            EXPECT_GT(tree.locate(removing), removing);
        }

        for (const T e : values) {
            EXPECT_EQ(removed.count(e), 1) << "Lost element " << e;
        }
    });

    producer.join();
    consumer.join();
}

TYPED_TEST(ParallelTest, ClusterSingleProducerSingleConsumer) {
    constexpr std::size_t N = 10'000'000;
    constexpr int TESTS = 8;
    using T = typename TestFixture::T;
    using UT = typename std::make_unsigned_t<T>;
    constexpr T SEED = 30;
    constexpr T MIN = std::numeric_limits<T>::min();
    constexpr T MAX = std::numeric_limits<T>::max();
    constexpr T POISON_PILL = MAX - 1;
    constexpr T MAX_VALUE = MAX - 2;
    constexpr size_t CLUSTER_SIZE = 1000;
    TypeParam tree;
    std::vector<T> values = generateClusterIntVector(N, SEED, MIN, MAX_VALUE, CLUSTER_SIZE);
    std::unordered_multiset<T> set(values.begin(), values.end());
    std::unordered_multiset<T> removed;

    std::thread producer([&]{
        for (const T e : values) {
            tree.insert(e);
        }
        tree.insert(POISON_PILL);
    });

    std::thread consumer([&]{
        UniformIntGenerator<T> gen(SEED + 1, MIN, MAX_VALUE);
        T results[TESTS];
        bool poisoned = false;
        while (true) {
            for (int i = 0; i < TESTS; ++i) {
                const T test = gen.rand();
                const T result = tree.locate(test);
                if (result == POISON_PILL) {
                    poisoned = true;
                    break;
                } else if (result != MAX) {
                    EXPECT_NE(set.find(result), set.end());
                    EXPECT_GE(result, test);
                }
                results[i] = result;
            }
            if (poisoned)
                break;
            const T removing = results[(UT)gen.rand() % TESTS];
            if (removing != MAX) {
                removed.insert(removing);
                EXPECT_LE(removed.count(removing), set.count(removing));
                tree.remove(removing);
                EXPECT_GE(tree.locate(removing), removing); // producer might have added it again => not _GT
            }
        }
        // Producer is done
        EXPECT_EQ(tree.locate(POISON_PILL), POISON_PILL);
        tree.remove(POISON_PILL);
        EXPECT_EQ(tree.locate(POISON_PILL), MAX);
        T removing = MIN;
        while ((removing = tree.locate(removing)) != MAX) {
            removed.insert(removing);
            EXPECT_LE(removed.count(removing), set.count(removing));
            EXPECT_EQ(tree.locate(removing), removing);
            tree.remove(removing);
            EXPECT_GT(tree.locate(removing), removing);
        }

        for (const T e : values) {
            EXPECT_GT(removed.count(e), 0) << "Lost element " << e;
        }
    });

    producer.join();
    consumer.join();
}

TYPED_TEST(ParallelTest, ClusterUniqueSingleProducerSingleConsumer) {
    constexpr std::size_t N = 10'000'000;
    constexpr int TESTS = 8;
    using T = typename TestFixture::T;
    using UT = typename std::make_unsigned_t<T>;
    constexpr T SEED = 35;
    constexpr T MIN = std::numeric_limits<T>::min();
    constexpr T MAX = std::numeric_limits<T>::max();
    constexpr T POISON_PILL = MAX - 1;
    constexpr T MAX_VALUE = MAX - 2;
    constexpr size_t CLUSTER_SIZE = 1000;
    TypeParam tree;
    std::vector<T> values = generateUniqueClusterIntVector(N, SEED, MIN, MAX_VALUE, CLUSTER_SIZE);
    std::unordered_set<T> set(values.begin(), values.end());
    std::unordered_set<T> removed;

    std::thread producer([&]{
        for (const T e : values) {
            tree.insert(e);
        }
        tree.insert(POISON_PILL);
    });

    std::thread consumer([&]{
        UniformIntGenerator<T> gen(SEED + 1, MIN, MAX_VALUE);
        T results[TESTS];
        bool poisoned = false;
        while (true) {
            for (int i = 0; i < TESTS; ++i) {
                const T test = gen.rand();
                const T result = tree.locate(test);
                if (result == POISON_PILL) {
                    poisoned = true;
                    break;
                } else if (result != MAX) {
                    EXPECT_NE(set.find(result), set.end());
                    EXPECT_GE(result, test);
                }
                results[i] = result;
            }
            if (poisoned)
                break;
            const T removing = results[(UT)gen.rand() % TESTS];
            if (removing != MAX) {
                EXPECT_EQ(removed.count(removing), 0);
                removed.insert(removing);
                tree.remove(removing);
                EXPECT_GT(tree.locate(removing), removing);
            }
        }
        // Producer is done
        EXPECT_EQ(tree.locate(POISON_PILL), POISON_PILL);
        tree.remove(POISON_PILL);
        EXPECT_EQ(tree.locate(POISON_PILL), MAX);
        T removing = MIN;
        while ((removing = tree.locate(removing)) != MAX) {
            EXPECT_EQ(tree.locate(removing), removing);
            removed.insert(removing);
            tree.remove(removing);
            EXPECT_GT(tree.locate(removing), removing);
        }

        for (const T e : values) {
            EXPECT_EQ(removed.count(e), 1) << "Lost element " << e;
        }
    });

    producer.join();
    consumer.join();
}

// The following two tests test a part of sequential consistency which
// is provided by all implementations, even the not sequentially consistent ones.
TYPED_TEST(ParallelTest, BackwardsInsertSingleProducer) {
    constexpr std::size_t N = 10'000'000;
    using T = typename TestFixture::T;
    constexpr T MIN_T = std::numeric_limits<T>::min();
    constexpr T MIN = 35;
    constexpr T MAX = MIN + N - 1;
    TypeParam tree;

    std::thread producer([&]{
        for (T i = MAX; i >= MIN; --i) {
            tree.insert(i);
        }
    });

    std::thread reader([&]{
        T e;
        T last = MAX;
        while (tree.empty())
            ;

        do {
            e = tree.locate(MIN_T);
            EXPECT_LE(e, last);
            EXPECT_LE(e, MAX);
            EXPECT_GE(e, MIN);
            last = e;
        } while (e != MIN);
    });

    std::thread remover([&]{
        T e;
        T last = MAX;
        while (tree.empty())
            ;

        do {
            e = tree.locate(MIN_T);
            EXPECT_LE(e, last);
            EXPECT_LE(e, MAX);
            EXPECT_GE(e, MIN);
            if (e != last) {
                tree.remove(last);
                last = e;
            }
        } while (e != MIN);
    });

    producer.join();
    reader.join();
    remover.join();
}

TYPED_TEST(ParallelTest, BackwardsInsertTwoProducers) {
    using T = typename TestFixture::T;
    constexpr T N = 10'000'000; // per producer
    constexpr T MIN = 35;
    constexpr T MAX = MIN + N - 1;
    TypeParam tree;

    std::thread producers[2];
    std::thread readers[2];
    std::thread removers[2];

    for (T i = 0; i < 2; ++i) {
        producers[i] = std::thread([&, i]{
            for (T e = MAX + i*N; e >= MIN + i*N; --e) {
                tree.insert(e);
            }
        });
    }

    for (T i = 0; i < 2; ++i) {
        readers[i] = std::thread([&, i]{
            T e;
            T last = MAX + i*N;
            while (tree.locate(MIN + i*N) > MAX + i*N)
                ;

            do {
                e = tree.locate(MIN + i*N);
                EXPECT_LE(e, last);
                EXPECT_LE(e, MAX + i*N);
                EXPECT_GE(e, MIN + i*N);
                last = e;
            } while (e != MIN + i*N);
        });
    }

    for (T i = 0; i < 2; ++i) {
        removers[i] = std::thread([&, i]{
            T e;
            T last = MAX + i*N;
            while (tree.locate(MIN + i*N) > MAX + i*N)
                ;

            do {
                e = tree.locate(MIN + i*N);
                EXPECT_LE(e, last);
                EXPECT_LE(e, MAX + i*N);
                EXPECT_GE(e, MIN + i*N);
                if (e != last) {
                    tree.remove(last);
                    last = e;
                }
            } while (e != MIN + i*N);
        });
    }

    for (int i = 0; i < 2; ++i) {
        producers[i].join();
        readers[i].join();
        removers[i].join();
    }
}

TYPED_TEST(ParallelTest, UniqueMultipleProducersMultipleReaders) {
    constexpr std::size_t N = 3'000'000; // per producer
    constexpr std::size_t NUM_PRODUCERS = 4;
    constexpr std::size_t NUM_READERS = 4;
    constexpr int TESTS = 8; // for remover
    using T = typename TestFixture::T;
    using UT = typename std::make_unsigned_t<T>;
    constexpr T SEED = 40;
    constexpr T MIN = std::numeric_limits<T>::min();
    constexpr T MAX = std::numeric_limits<T>::max();
    constexpr T POISON_PILL = MAX - 1;
    constexpr T MAX_VALUE = MAX - 2;
    TypeParam tree;
    std::vector<T> values = generateUniqueUniformIntVector(N * NUM_PRODUCERS, SEED, MIN, MAX_VALUE);
    const std::unordered_set<T> set(values.begin(), values.end());

    std::thread producers[NUM_PRODUCERS];
    std::thread readers[NUM_READERS];

    for (std::size_t i = 0; i < NUM_PRODUCERS; ++i) {
        producers[i] = std::thread([&, i]{
            for (std::size_t j = i * N; j < (i + 1) * N; ++j) {
                tree.insert(values[j]);
            }
        });
    }

    for (std::size_t i = 0; i < NUM_READERS; ++i) {
        readers[i] = std::thread([&, i]{
            UniformIntGenerator<T> gen(SEED + 2 + i, MIN, MAX_VALUE);
            while (true) {
                const T test = gen.rand();
                const T result = tree.locate(test);
                if (result == POISON_PILL) {
                    break;
                } else if (result != MAX) {
                    EXPECT_NE(set.find(result), set.end());
                    EXPECT_GE(result, test);
                }
            }
        });
    }

    std::thread remover([&]{
        std::unordered_set<T> removed;
        UniformIntGenerator<T> gen(SEED + 1, MIN, MAX_VALUE);
        T results[TESTS];
        bool poisoned = false;
        while (true) {
            for (int i = 0; i < TESTS; ++i) {
                const T test = gen.rand();
                const T result = tree.locate(test);
                if (result == POISON_PILL) {
                    poisoned = true;
                    break;
                } else if (result != MAX) {
                    EXPECT_NE(set.find(result), set.end());
                    EXPECT_GE(result, test);
                }
                results[i] = result;
            }
            if (poisoned)
                break;
            const T removing = results[(UT)gen.rand() % TESTS];
            if (removing < POISON_PILL) {
                EXPECT_EQ(removed.count(removing), 0);
                removed.insert(removing);
                tree.remove(removing);
                EXPECT_GT(tree.locate(removing), removing);
            }
        }
        // Producers are done
        EXPECT_EQ(tree.locate(POISON_PILL), POISON_PILL);
        T removing = MIN;
        while ((removing = tree.locate(removing)) < POISON_PILL) {
            EXPECT_NE(set.find(removing), set.end());
            EXPECT_EQ(tree.locate(removing), removing);
            removed.insert(removing);
            tree.remove(removing);
            EXPECT_GT(tree.locate(removing), removing);
        }

        for (const T e : values) {
            EXPECT_EQ(removed.count(e), 1) << "Lost element " << e;
        }
    });

    for (std::size_t i = 0; i < NUM_PRODUCERS; ++i) {
        producers[i].join();
    }
    tree.insert(POISON_PILL);
    for (std::size_t i = 0; i < NUM_READERS; ++i) {
        readers[i].join();
    }
    remover.join();
}

TYPED_TEST(ParallelTest, UniqueMultipleProducersMultipleRemovers) {
    constexpr std::size_t N = 3'000'000; // per producer
    constexpr std::size_t NUM_PRODUCERS = 4;
    constexpr std::size_t NUM_REMOVERS = 4;
    using T = typename TestFixture::T;
    constexpr T SEED = 50;
    constexpr T MIN = std::numeric_limits<T>::min();
    constexpr T MAX = std::numeric_limits<T>::max();
    constexpr T POISON_PILL = MAX - 1;
    constexpr T MAX_VALUE = MAX - 2;
    TypeParam tree;
    std::vector<T> values = generateUniqueUniformIntVector(N * NUM_PRODUCERS, SEED, MIN, MAX_VALUE);
    const std::unordered_set<T> set(values.begin(), values.end());
    std::unordered_set<T> removed;
    std::mutex removedMutex;

    std::thread producers[NUM_PRODUCERS];
    std::thread removers[NUM_REMOVERS];

    for (std::size_t i = 0; i < NUM_PRODUCERS; ++i) {
        producers[i] = std::thread([&, i]{
            for (std::size_t j = i * N; j < (i + 1) * N; ++j) {
                tree.insert(values[j]);
            }
        });
    }

    for (std::size_t i = 0; i < NUM_REMOVERS; ++i) {
        removers[i] = std::thread([&, i]{
            UniformIntGenerator<T> gen(SEED + 1 + i, MIN, MAX_VALUE);
            std::vector<T> localRemoved;
            while (true) {
                const T test = gen.rand();
                const T result = tree.locate(test);
                if (result == POISON_PILL) {
                    break;
                } else if (result != MAX) {
                    EXPECT_NE(set.find(result), set.end());
                    EXPECT_GE(result, test);
                    if (result % NUM_REMOVERS == i) {
                        localRemoved.push_back(result);
                        tree.remove(result);
                        EXPECT_GT(tree.locate(result), result);
                    }
                }
            }
            // Producers are done
            EXPECT_EQ(tree.locate(POISON_PILL), POISON_PILL);
            T removing = MIN;
            while ((removing = tree.locate(removing)) < POISON_PILL) {
                EXPECT_NE(set.find(removing), set.end());
                if (removing % NUM_REMOVERS == i) {
                    EXPECT_EQ(tree.locate(removing), removing);
                    localRemoved.push_back(removing);
                    tree.remove(removing);
                    EXPECT_GT(tree.locate(removing), removing);
                }
                ++removing;
            }

            std::lock_guard lock(removedMutex);
            for (const T e : localRemoved) {
                EXPECT_EQ(removed.find(e), removed.end());
                removed.insert(e);
            }
        });
    }

    for (std::size_t i = 0; i < NUM_PRODUCERS; ++i) {
        producers[i].join();
    }
    tree.insert(POISON_PILL);
    for (std::size_t i = 0; i < NUM_REMOVERS; ++i) {
        removers[i].join();
    }
    for (const T e : values) {
        EXPECT_EQ(removed.count(e), 1) << "Lost element " << e;
    }
}

template<class Tree>
class ParallelSequentialConsistencyTest : public ::testing::Test {
protected:
    static constexpr int K = Tree::BitsK;
    using T = typename Tree::TypeT;
};
using ParallelSequentialConsistencyTestTypes = ::testing::Types<
    VEBL<>, VEBL<int32_t>,
    VEB32L<>, VEB32L<int32_t>,
    VEB32LT<>, VEB32LT<int32_t> // The others are not sequentially consistent
    >;
TYPED_TEST_SUITE(ParallelSequentialConsistencyTest, ParallelSequentialConsistencyTestTypes);

TYPED_TEST(ParallelSequentialConsistencyTest, SingleProducerMultipleReaders) {
    constexpr std::size_t NUM_READERS = 6;
    using T = typename TestFixture::T;
    constexpr T SPACING = 1U << 14; // bottom size (one element for each bottom)
    constexpr T MIN = std::numeric_limits<T>::min();
    constexpr T MAX = std::numeric_limits<T>::max();
    constexpr T MAX_VALUE = MAX - 1;
    constexpr T END = MAX - SPACING;
    TypeParam tree;
    tree.insert(MAX_VALUE); // Such that locate does not return immediately

    std::thread producer([&]{
        for (T e = MIN; e < END; e += SPACING) {
            tree.insert(e);
        }
    });

    std::thread readers[NUM_READERS];
    for (std::size_t i = 0; i < NUM_READERS; ++i) {
        readers[i] = std::thread([&]{
            for (T e = MIN; e < END; e += SPACING) {
                T result;
                while ((result = tree.locate(e)) != e) {
                    ASSERT_EQ(result, MAX_VALUE) << "For element " << e;
                }
            }
        });
    }

    producer.join();
    for (std::size_t i = 0; i < NUM_READERS; ++i) {
        readers[i].join();
    }
}

TYPED_TEST(ParallelSequentialConsistencyTest, SingleProducerMultipleReadersHalfSpacing) {
    constexpr std::size_t NUM_READERS = 6;
    using T = typename TestFixture::T;
    constexpr T SPACING = 1U << 13; // half bottom size (two elements for each bottom)
    constexpr T MIN = std::numeric_limits<T>::min();
    constexpr T MAX = std::numeric_limits<T>::max();
    constexpr T MAX_VALUE = MAX - 1;
    constexpr T END = MAX - SPACING;
    TypeParam tree;
    tree.insert(MAX_VALUE); // Such that locate does not return immediately

    std::thread producer([&]{
        for (T e = MIN; e < END; e += SPACING) {
            tree.insert(e);
        }
    });

    std::thread readers[NUM_READERS];
    for (std::size_t i = 0; i < NUM_READERS; ++i) {
        readers[i] = std::thread([&]{
            for (T e = MIN; e < END; e += SPACING) {
                T result;
                while ((result = tree.locate(e)) != e) {
                    ASSERT_EQ(result, MAX_VALUE) << "For element " << e;
                }
            }
        });
    }

    producer.join();
    for (std::size_t i = 0; i < NUM_READERS; ++i) {
        readers[i].join();
    }
}
