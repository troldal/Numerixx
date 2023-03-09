//
// Created by Kenneth Balslev on 09/03/2023.
//

#include <benchmark/benchmark.h>
#include <numerixx.hpp>

using namespace nxx::deriv;
static auto func = [](double x) { return std::log(x) + 2 * x; };

static void BM_Order1CentralRichardson(benchmark::State& state) {
    double result;
    for (auto _ : state) {
        result = *derivative< Order1CentralRichardson >(func, std::numbers::e);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order1CentralRichardson);

static void BM_Order1Central3Point(benchmark::State& state) {
    double result;
    for (auto _ : state) {
        result = *derivative< Order1Central3Point >(func, std::numbers::e);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order1Central3Point);

static void BM_Order1Central5Point(benchmark::State& state) {
    double result;
    for (auto _ : state) {
        result = *derivative< Order1Central5Point >(func, std::numbers::e);
        benchmark::ClobberMemory();
    }
}
// Register the function as a benchmark
BENCHMARK(BM_Order1Central5Point);

BENCHMARK_MAIN();