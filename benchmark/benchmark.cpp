//
// Created by Kenneth Balslev on 09/03/2023.
//

#include <benchmark/benchmark.h>
#include <numerixx.hpp>

using namespace nxx::deriv;
static auto func = [](double x) { return std::log(x) + 2 * x; };

static void BM_Order1CentralRichardson(benchmark::State& state) {
    for (auto _ : state) {
        auto result = *derivative< Order1CentralRichardson >(func, std::numbers::e);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order1CentralRichardson);

static void BM_Order1Central3Point(benchmark::State& state) {
    for (auto _ : state) {
        auto result = *derivative< Order1Central3Point >(func, std::numbers::e);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order1Central3Point);

static void BM_Order1Central5Point(benchmark::State& state) {
    for (auto _ : state) {
        auto result = *derivative< Order1Central5Point >(func, std::numbers::e);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order1Central5Point);


static void BM_Order1ForwardRichardson(benchmark::State& state) {
    for (auto _ : state) {
        auto result = *derivative< Order1ForwardRichardson >(func, std::numbers::e);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order1ForwardRichardson);

static void BM_Order1Forward2Point(benchmark::State& state) {
    for (auto _ : state) {
        auto result = *derivative< Order1Forward2Point >(func, std::numbers::e);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order1Forward2Point);

static void BM_Order1Forward3Point(benchmark::State& state) {
    for (auto _ : state) {
        auto result = *derivative< Order1Forward3Point >(func, std::numbers::e);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order1Forward3Point);


static void BM_Order1BackwardRichardson(benchmark::State& state) {
    for (auto _ : state) {
        auto result = *derivative< Order1BackwardRichardson >(func, std::numbers::e);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order1BackwardRichardson);

static void BM_Order1Backward2Point(benchmark::State& state) {
    for (auto _ : state) {
        auto result = *derivative< Order1Backward2Point >(func, std::numbers::e);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order1Backward2Point);

static void BM_Order1Backward3Point(benchmark::State& state) {
    for (auto _ : state) {
        auto result = *derivative< Order1Backward3Point >(func, std::numbers::e);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order1Backward3Point);



static void BM_Order2Central3Point(benchmark::State& state) {
    for (auto _ : state) {
        auto result = *derivative< Order2Central3Point >(func, std::numbers::e);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order2Central3Point);

static void BM_Order2Central5Point(benchmark::State& state) {
    for (auto _ : state) {
        auto result = *derivative< Order2Central5Point >(func, std::numbers::e);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order2Central5Point);

static void BM_Order2Forward3Point(benchmark::State& state) {
    for (auto _ : state) {
        auto result = *derivative< Order2Forward3Point >(func, std::numbers::e);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order2Forward3Point);

static void BM_Order2Forward4Point(benchmark::State& state) {
    for (auto _ : state) {
        auto result = *derivative< Order2Forward4Point >(func, std::numbers::e);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order2Forward4Point);

static void BM_Order2Backward3Point(benchmark::State& state) {
    for (auto _ : state) {
        auto result = *derivative< Order2Backward3Point >(func, std::numbers::e);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order2Backward3Point);

static void BM_Order2Backward4Point(benchmark::State& state) {
    for (auto _ : state) {
        auto result = *derivative< Order2Backward4Point >(func, std::numbers::e);
        benchmark::DoNotOptimize(result);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order2Backward4Point);