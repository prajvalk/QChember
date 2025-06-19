#ifndef NEWSCF_TESTINGAPI_HPP
#define NEWSCF_TESTINGAPI_HPP

#include <iostream>
#include <string>
#include <cmath>

namespace newscf::testing {

    struct TestHandle {
        std::string testsuite;
        int passed_tests = 0;
        int failed_tests = 0;
        int exitcode;

        int assertions = 0;
        std::string current_test;
        bool is_test_running = false;
        bool fail_flag = false;
        int failed_assertions = 0;
        std::string err_msg;
    };

    inline void init_tests() {
        std::cout << "==== new_scf Testing Toolkit \n"; 
    }

    inline void complete_test (TestHandle& handle) {
        if(!handle.fail_flag) {
            std::cout << "PASSED. " << handle.assertions << " assertions passed out of " << handle.assertions << "\n";
            handle.passed_tests++;
        } else {
            std::cout << "FAILED. " << handle.assertions << " assertions passed out of " << handle.assertions+handle.failed_assertions << "\n";
            std::cout << "==== Failed assertions: \n";
            std::cout << handle.err_msg << "\n";
            handle.failed_tests++;
        }

        handle.assertions = 0;
        handle.current_test = "";
        handle.is_test_running = false;
        handle.fail_flag = false;
        handle.failed_assertions = 0;
        handle.err_msg = "";
    }

    inline void add_test(std::string name, TestHandle& handle) {
        if(handle.is_test_running) complete_test(handle);
        handle.current_test = name;
        handle.is_test_running = true;

        std::cout << "==== Running " << name << " Test... ";
    }

    template <typename T>
    inline void test_assert_eq (const T a, const T b, const T tol, TestHandle& handle, std::string file, int lineno) {
        if (a == b) {
            handle.assertions++;
        } else {
            handle.err_msg += "\t "+file+":"+std::to_string(lineno)+" | Assertion failed at between " + std::to_string(a) + " and " + std::to_string(b) + " with tol=" + std::to_string(tol) + "\n";
            handle.failed_assertions++;
            handle.fail_flag = true;
        }
    }

    template <typename T>
    inline void test_assert_eq (const T a, const T b, TestHandle& handle, std::string file, int lineno) {
        test_assert_eq (a, b, T(0), handle, file, lineno);
    }

    #define TEST_ASSERT_EQ(a, b, h) \
        test_assert_eq (a, b, h, __FILE__, __LINE__);

    #define TEST_ASSERT_EQ_TOL(a, b, tol, h) \
        test_assert_eq (a, b, tol, h, __FILE__, __LINE__);

    inline void end_tests(TestHandle& handle) {
        if(handle.is_test_running) complete_test(handle);

        std::cout << "==== \n";
        std::cout << "==== TESTSUITE " << handle.testsuite << "\n";
        std::cout << "==== Passed Tests: " << handle.passed_tests << "\n";
        std::cout << "==== Failed Tests: " << handle.failed_tests << "\n";
        std::cout << "==== ";

        handle.passed_tests = 0;
        handle.failed_tests = 0;

        std::cout << "\n";

        if (handle.failed_tests >= 0) handle.exitcode = -1 - handle.failed_tests;
    }

}

#endif