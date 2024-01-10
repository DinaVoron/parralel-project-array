#pragma once
#include <thread>
#include <omp.h>

const uint32_t A = 22695477;
const uint32_t B = 1;


static unsigned g_thread_num = std::thread::hardware_concurrency();


unsigned get_num_threads() {
	return g_thread_num;
}

void set_num_threads(unsigned T) {
	g_thread_num = T;
	omp_set_num_threads(T);
}

struct profiling_results_t {
	double result;
	double time;
	double speedup, efficiency;
	unsigned T;
};

class lc_t {
public:
    uint32_t A = 1, B = 0;
    lc_t(uint32_t a = 1, uint32_t b = 0) : A(a), B(b) {}

    lc_t& operator *= (const lc_t& x) {
        if (A == 1 && B == 0) {
            A = x.A;
            B = x.B;
        }
        else {
            A *= x.A;
            B += A * x.B;
        }

        return *this;
    }


    auto operator() (uint32_t seed, uint32_t min_value, uint32_t max_value) const {
        if (max_value - min_value + 1 != 0) {
            return (A * seed + B) % (max_value - min_value) + min_value;
        }
        else {
            return A * seed + B;
        }
    }

};

template <typename F>
auto run_experiment(F f, const uint32_t* v, std::size_t n, uint32_t seed)
    requires std::is_invocable_r_v<F, const uint32_t*, std::size_t, uint32_t> {
    std::vector<profiling_results_t> r;
    std::size_t T_max = get_num_threads();
    for (std::size_t T = 1; T <= T_max; ++T) {
        set_num_threads(T);
        using namespace std::chrono;
        auto t0 = std::chrono::steady_clock::now();
        auto rr_result = f(v, n, seed);
        auto t1 = std::chrono::steady_clock::now();
        profiling_results_t rr;
        r.push_back(rr);
        unsigned int times = duration_cast<nanoseconds> (t1 - t0).count();
        r[T - 1].time = times;
        r[T - 1].result = rr_result;
        r[T - 1].T = T;
        r[T - 1].speedup = r[0].time / r[T - 1].time;
        r[T - 1].efficiency = r[T - 1].speedup / T;
    }
    return r;
}