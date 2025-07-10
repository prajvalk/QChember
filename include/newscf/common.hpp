/**
 * NewSCF Library | Commons Header
 * Mozilla Public License 2.0
 * Copyright (C) 2025, Prajval K
 */

#ifndef NEWSCF_COMMON_HPP
#define NEWSCF_COMMON_HPP

#ifndef NEWSCF_DISABLE_LOGGING
#include <string>
#include <cstring>
#include <iostream>
#endif

#include <unordered_map>
#include <stdexcept>
#include <cctype>
#include <algorithm>
#include <cmath>

namespace newscf {

// Logging Kit
#ifndef NEWSCF_DISABLE_LOGGING

#ifndef NEWSCF_LOG_LEVEL
#define NEWSCF_LOG_LEVEL WARN
#endif

#define LOG(LEVEL, MSG)            newscf_log(LEVEL, MSG, __FILE__, __LINE__)
#define LOG_ASIS(LEVEL, MSG, F, L) newscf_log(LEVEL, MSG, F, L)

    enum NEWSCF_LOGLEVELS {
        DEV_DUMP,
        DEV_INFO,
        DEV_WARN,
        DEV_ERROR,
        INFO,
        WARN,
        ERROR
    };

    inline void newscf_log (const NEWSCF_LOGLEVELS level, const std::string& msg, const std::string& file, const int line) {
        if       (level < NEWSCF_LOG_LEVEL)  return;

        std::string level_str;

        if       (level == DEV_DUMP)         level_str = "DEV_DUMP";
        else if  (level == DEV_INFO)         level_str = "DEV_INFO";
        else if  (level == DEV_WARN)         level_str = "DEV_WARN";
        else if  (level == DEV_ERROR)        level_str = "DEV_ERROR";
        else if  (level == INFO)             level_str = "INFO";
        else if  (level == WARN)             level_str = "WARN";
        else                                          level_str = "ERROR";

        std::cout << "==== "+level_str+" ==== ("+file+":"+std::to_string(line)+"): "+msg << std::endl;
    }

#else
#define LOG(LEVEL, MSG)
#define LOG_ASIS(LEVEL, MSG, F, L)
#endif

// Error Handling
#define HANDLE_ERROR(MSG, CODE) handle_error (MSG, CODE, __FILE__, __LINE__);

#define HANDLE_EXCEPTION(EXC, MSG, CODE) handle_error (EXC, MSG, CODE, __FILE__, __LINE__);

    inline void handle_error (const std::string &msg, int exitcode, const std::string& fn, const int no) {
        LOG_ASIS (ERROR, msg, fn, no);
        throw std::runtime_error(msg);
        exit(exitcode);
    }

    inline void handle_error (const std::exception& ex, const std::string &msg, int exitcode, const std::string& fn, const int no) {
        LOG_ASIS (ERROR, msg, fn, no);
        LOG_ASIS (DEV_ERROR, ex.what(), fn, no);
        throw std::runtime_error(msg);
        throw ex;
        exit(exitcode);
    }

// Memory Toolkit

    inline void* malloc(const size_t size) {
        LOG (DEV_INFO, "Allocating "+std::to_string(size)+" bytes of memory.");
        return std::malloc(size);
    }

    inline void* memset(void* ptr, const int value, const size_t size) {
        LOG(DEV_INFO, "Setting "+std::to_string(size)+" bytes of memory.");
        return std::memset(ptr, value, size);
    }

    inline void memcpy(void* dest, const void* src, const size_t size) {
        LOG(DEV_INFO, "Copying "+std::to_string(size)+" bytes of memory.");
        std::memcpy(dest, src, size);
    }

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

// Utility Toolkit

    inline std::string normalize_symbol(const std::string& symbol) {
        if (symbol.empty()) return symbol;
        std::string norm;
        norm += static_cast<char>(std::toupper(symbol[0]));
        for (size_t i = 1; i < symbol.size(); ++i)
            norm += symbol[i];
        return norm;
    }

    static const std::unordered_map<std::string, int> periodic_table = {
        {"H", 1}, {"He", 2},
        {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10},
        {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18},
        {"K", 19}, {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24}, {"Mn", 25}, {"Fe", 26},
        {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}, {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34},
        {"Br", 35}, {"Kr", 36}, {"Rb", 37}, {"Sr", 38}, {"Y", 39}, {"Zr", 40}, {"Nb", 41}, {"Mo", 42},
        {"Tc", 43}, {"Ru", 44}, {"Rh", 45}, {"Pd", 46}, {"Ag", 47}, {"Cd", 48}, {"In", 49}, {"Sn", 50},
        {"Sb", 51}, {"Te", 52}, {"I", 53}, {"Xe", 54}, {"Cs", 55}, {"Ba", 56}, {"La", 57}, {"Ce", 58},
        {"Pr", 59}, {"Nd", 60}, {"Pm", 61}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64}, {"Tb", 65}, {"Dy", 66},
        {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70}, {"Lu", 71}, {"Hf", 72}, {"Ta", 73}, {"W", 74},
        {"Re", 75}, {"Os", 76}, {"Ir", 77}, {"Pt", 78}, {"Au", 79}, {"Hg", 80}, {"Tl", 81}, {"Pb", 82},
        {"Bi", 83}, {"Po", 84}, {"At", 85}, {"Rn", 86}, {"Fr", 87}, {"Ra", 88}, {"Ac", 89}, {"Th", 90},
        {"Pa", 91}, {"U", 92}, {"Np", 93}, {"Pu", 94}, {"Am", 95}, {"Cm", 96}, {"Bk", 97}, {"Cf", 98},
        {"Es", 99}, {"Fm", 100}, {"Md", 101}, {"No", 102}, {"Lr", 103}, {"Rf", 104}, {"Db", 105}, {"Sg", 106},
        {"Bh", 107}, {"Hs", 108}, {"Mt", 109}, {"Ds", 110}, {"Rg", 111}, {"Cn", 112}, {"Nh", 113}, {"Fl", 114},
        {"Mc", 115}, {"Lv", 116}, {"Ts", 117}, {"Og", 118}
    };

    constexpr double angstrom_to_bohr = 1.8897259886;

    inline int get_atomic_number(const std::string& symbol_raw) {
        const std::string symbol = normalize_symbol(symbol_raw);
        const auto it = periodic_table.find(symbol);
        if (it == periodic_table.end()) {
            HANDLE_ERROR ("Unknown atomic symbol: " + symbol_raw, 4);
        }
        return it->second;
    }

    inline const std::string& get_atomic_symbol(const int atomic_number) {
        // Static reverse map cache: atomic number -> symbol
        static std::unordered_map<int, std::string> number_to_symbol;

        // Build reverse map once on first call
        if (number_to_symbol.empty()) {
            for (const auto& [sym, num] : periodic_table) {
                number_to_symbol[num] = sym;
            }
        }

        const auto it = number_to_symbol.find(atomic_number);
        if (it == number_to_symbol.end()) {
            HANDLE_ERROR("Unknown atomic number: " + std::to_string(atomic_number), 4);
            static const std::string empty;
            return empty; // Just to satisfy return type, this won't be reached if HANDLE_ERROR throws or exits.
        }
        return it->second;
    }

    inline int get_angular_momentum_from_symbol (const char ch) {
        switch (ch) {
            case 'S': return 0;
            case 'P': return 1;
            case 'D': return 2;
            case 'F': return 3;
            case 'G': return 4;
            case 'H': return 5;
            case 'I': return 6;
            default:
                LOG (WARN, "Unsupported angular momentum level "+std::to_string(ch));
                HANDLE_ERROR ("Angular momentum too large! Unsupported!", 99002);
                return 99;
        }
    }

    inline bool starts_with(const std::string& str, const std::string& prefix) {
        return str.size() >= prefix.size() &&
               str.compare(0, prefix.size(), prefix) == 0;
    }


    // --- FORTRAN-style float parser ---
    inline double parse_fortran_double(const std::string& token) {
        std::string fixed_token = token;
        std::replace(fixed_token.begin(), fixed_token.end(), 'D', 'E');
        std::replace(fixed_token.begin(), fixed_token.end(), 'd', 'E');
        return std::stod(fixed_token);
    }

    inline long long fact2 (const int i) {
        if (i <= 1) return 1;
        return i * fact2(i - 2);
    }

}

#endif //NEWSCF_COMMON_HPP