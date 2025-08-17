/*
 * NewSCF
 * Copyright (C) 2025, Prajval K
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef NEWSCF_TESTING_COMMON_INTERFACE_SYSTEM_HPP
#define NEWSCF_TESTING_COMMON_INTERFACE_SYSTEM_HPP

#include <string>
#include <iostream>
#include <cmath>

namespace newscf::testing_cis {
     // Testing Toolkit

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
        std::cout << "==== NewSCF Testing Toolkit \n";
    }

    inline void complete_test (TestHandle& handle) {
        if(!handle.fail_flag) {
            std::cout << "\x1b[1;32mPASSED.\x1b[0m " << handle.assertions << " assertions passed out of " << handle.assertions << "\n";
            handle.passed_tests++;
        } else {
            std::cout << "\x1b[1;31mFAILED.\x1b[0m" << handle.assertions << " assertions passed out of " << handle.assertions+handle.failed_assertions << "\n";
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
        if (static_cast<T>(fabs(T(a) - T(b))) <= T(tol)) {
            handle.assertions++;
        } else {
            handle.err_msg += "\t "+file+":"+std::to_string(lineno)+" | Assertion failed at between " + std::to_string(a) + " and " + std::to_string(b) + " with tol=" + std::to_string(tol) + "\n";
            handle.failed_assertions++;
            handle.fail_flag = true;
        }
    }

    template <typename T>
    inline void test_assert_eq (const T a, const T b, TestHandle& handle, std::string file, int lineno) {
        if constexpr (std::is_same_v<T, double>) {
            test_assert_eq (a, b, T(1e-15), handle, file, lineno);
        } else if constexpr (std::is_same_v<T, float>) {
            test_assert_eq (a, b, T(1e-5), handle, file, lineno);
        } else {
            test_assert_eq (a, b, T(1e-16), handle, file, lineno);
        }
    }

    #define TEST_ASSERT_EQ(a, b, h) \
        test_assert_eq (a, b, h, __FILE__, __LINE__);

    #define TEST_ASSERT_EQ_TOL(a, b, tol, h) \
        test_assert_eq (a, b, tol, h, __FILE__, __LINE__);

    inline void end_tests(TestHandle& handle) {
        if(handle.is_test_running) complete_test(handle);

        std::cout << "==== \n";
        std::cout << "==== TESTSUITE " << handle.testsuite << "\n";
        std::cout << "==== \x1b[0;32mPassed Tests: " << handle.passed_tests << "\x1b[0m\n";
        std::cout << "==== \x1b[0;31mFailed Tests: " << handle.failed_tests << "\x1b[0m\n";
        std::cout << "==== ";

        if (handle.failed_tests > 0) handle.exitcode = -1 - handle.failed_tests;
        else if (handle.failed_tests == 0) handle.exitcode = 0;

        handle.passed_tests = 0;
        handle.failed_tests = 0;

        std::cout << "\n";
    }

}

#endif