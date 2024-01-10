#include<vector>
#include<iostream>
#include<concepts>
#include<stdint.h>
#include<thread>
#include<ranges>
#include<mutex>
#include"Header.h"




template <class T, std::unsigned_integral u>
auto my_pow(T x, u n) requires requires (T x) { T(1); x *= x; } {
    T r = T(1);
    while (n > 0) {
        if (n & 1)
            r *= x;
        x *= x;
        n = n >> 1;
    }
    return r;
}


double randomize_vector(std::vector<uint32_t>& V, size_t n, uint32_t seed, uint32_t min_val = 0, uint32_t max_val = UINT32_MAX) {
    double res = 0;
    if (max_val < min_val) {
        exit(__LINE__);
    }
    lc_t g = lc_t(A, B);
    for (int i = 0; i < n; ++i) {
        g *= g;
        V[i] = g(seed, min_val, max_val);
        res += V[i];
    }
    return res / n;
};


double randomize_vector_par(uint32_t* V, size_t n, uint32_t seed, uint32_t min_val = 0, uint32_t max_val = UINT32_MAX) {
    double res = 0;

    unsigned T = get_num_threads();
    std::vector<std::thread> workers;
    std::mutex mtx;

    auto worker_proc = [V, n, seed, T, min_val, max_val, &res, &mtx](unsigned t) {
        double partial = 0;

        size_t b = n % T, e = n / T;
        if (t < b)
            b = t * ++e;
        else
            b += t * e;
        e += b;


        auto generator = my_pow(lc_t(A, B), my_pow(2u, b + 1));
        for (int i = b; i < e; ++i) {
            generator *= generator;
            V[i] = generator(seed, min_val, max_val);
            partial += V[i];
        }

        {
            std::scoped_lock l{ mtx };
            res += partial;
        }

    };

    for (unsigned t = 1; t < T; ++t) {
        workers.emplace_back(worker_proc, t);
    }

    worker_proc(0);
    for (auto& w : workers) {
        w.join();
    }

    return res / n;
}


double randomize_vector_par(std::vector<uint32_t>& V, uint32_t seed, uint32_t min_val = 0, uint32_t max_val = UINT32_MAX) {
    return randomize_vector_par(V.data(), V.size(), seed, min_val, max_val);
}


bool randomize_test() {
    std::vector<uint32_t> V1(100), V2(100);
    if (randomize_vector(V1, 0, A, B) != randomize_vector_par(V2, 0, A, B)) {
        return false;
    }
    auto pr = std::ranges::mismatch(V1, V2);
    return pr.in1 == V1.end() && pr.in2 == V2.end();
}


int main() {

    {
        const size_t N = 1u << 10;
        std::vector<uint32_t> arr(N);

        std::size_t T_max = get_num_threads();
        std::vector<profiling_results_t> res(T_max);

        for (unsigned T = 1; T <= T_max; ++T) {
            set_num_threads(T);

            res[T - 1].T = T;

            auto t1 = std::chrono::steady_clock::now();
            res[T - 1].result = randomize_vector_par(arr, 120);
            auto t2 = std::chrono::steady_clock::now();

            res[T - 1].time = duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
            res[T - 1].speedup = res[0].time / res[T - 1].time;
            res[T - 1].efficiency = res[T - 1].speedup / T;
        }

        std::cout << "Result,Time,Speedup,Efficiency(average_cpp_mtx)" << "\n";
        for (int i = 0; i < res.size(); i++) {
            std::cout << res[i].result << "," << res[i].time << "," << res[i].speedup << "," << res[i].efficiency << "\n";
        }
    }
    
    {
        const size_t N = 1u << 10;
        std::vector<uint32_t> arr(N);

        std::size_t T_max = get_num_threads();
        std::vector<profiling_results_t> res(T_max);

        for (unsigned T = 1; T <= T_max; ++T) {
            set_num_threads(T);

            res[T - 1].T = T;

            auto t1 = std::chrono::steady_clock::now();
            res[T - 1].result = randomize_vector(arr, N, 120);
            auto t2 = std::chrono::steady_clock::now();

            res[T - 1].time = duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
            res[T - 1].speedup = res[0].time / res[T - 1].time;
            res[T - 1].efficiency = res[T - 1].speedup / T;
        }

        std::cout << "Result,Time,Speedup,Efficiency(average_cpp_mtx)" << "\n";
        for (int i = 0; i < res.size(); i++) {
            std::cout << res[i].result << "," << res[i].time << "," << res[i].speedup << "," << res[i].efficiency << "\n";
        }

    }
}