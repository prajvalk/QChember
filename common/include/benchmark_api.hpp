#ifndef NEWSCF_BENCHMARK_API
#define NEWSCF_BENCHMARK_API

#include <string>
#include <chrono>
#include <iostream>
#include <fstream>

namespace newscf::benchmark {

    using namespace std::chrono;

    struct BenchmarkHandle {
        std::string testsuite;
        std::string current_test;
        bool isTestRunning = false;
        std::ofstream file;

        high_resolution_clock::time_point before;
        high_resolution_clock::time_point after;
    };

    void init_handle(BenchmarkHandle h) {
        std::cout << "==== new_scf Benchmarking Toolkit \n";
        std::cout << "==== testsuite: " << h.testsuite << " \n";
        std::cout << "==== Generates a CSV compatible report by default. \n";
        h.file = std::ofstream (h.testsuite+"_benchmark.csv");
        h.file << "testcase, time_ms \n";
    }

    void start_time(BenchmarkHandle h, std::string test) { 
        h.file << "==== Running " << test << " ... ";
        h.isTestRunning = true;
        h.before = high_resolution_clock::now();
    }

    void end_time(BenchmarkHandle h) {
        h.after = high_resolution_clock::now();
        high_resolution_clock::duration diff = h.after - h.before;
        auto diff_ms = duration_cast<milliseconds> (diff);
        std::cout << "Done. Took "+diff_ms.count();
        h.file << h.current_test << ", " << diff_ms.count() << "\n";
        h.isTestRunning = false; 
    }

    void close_handle(BenchmarkHandle h) {
        /*if(h.isTestRunning) 
            end_time(h);*/
        h.file.close();
        std::cout << "=== Benchmark complete. \n";
        std::cout << "=== Please send the report " << h.testsuite << "_benchmark.csv to the developers with your new_scf build tag if possible. \n";
    }

}

#endif