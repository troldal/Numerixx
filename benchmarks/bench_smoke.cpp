// Placeholder benchmark that keeps the google/benchmark wiring built. The root-finding benchmarks, which report
// function evaluations as a counter, arrive in phase 3.
#include <numerixx/numerixx.hpp>

#include <benchmark/benchmark.h>

namespace
{
    void bm_version(benchmark::State& state)
    {
        int major = nxx::version.major;
        for (auto _ : state) benchmark::DoNotOptimize(major);
        state.counters["fevals"] = 0;
    }
}    // namespace

BENCHMARK(bm_version);
