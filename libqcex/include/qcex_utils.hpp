#ifndef NEWSCF_QCEX_UTILS_HPP
#define NEWSCF_QCEX_UTILS_HPP

#include <algorithm>
#include <string>

#include "logging_api.hpp"

inline bool is_number(const std::string &s) {
  return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}

#define HANDLE_ERROR(MSG, CODE) \
    handle_error (MSG, CODE, __FILE__, __LINE__);

#define HANDLE_EXCEPTION(EXC, MSG, CODE) \
    handle_error (EXC, MSG, CODE, __FILE__, __LINE__);

inline void handle_error (const std::string &msg, int exitcode, std::string fn, int no) {
    LOG (ERROR, msg);
    throw std::runtime_error(msg);
    exit(exitcode);
}

inline void handle_error (const std::exception& ex, const std::string &msg, int exitcode, std::string fn, int no) {
    LOG (ERROR, msg);
    LOG (DEV_ERROR, ex.what());
    throw std::runtime_error(msg);
    throw ex;
    exit(exitcode);
}


#include <unordered_map>
#include <stdexcept>
#include <cctype>

inline std::string normalize_symbol(const std::string& symbol) {
    if (symbol.empty()) return symbol;
    std::string norm;
    norm += std::toupper(symbol[0]);
    for (size_t i = 1; i < symbol.size(); ++i)
        norm += std::tolower(symbol[i]);
    return norm;
}

inline int get_atomic_number(const std::string& symbol_raw) {
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

    std::string symbol = normalize_symbol(symbol_raw);
    auto it = periodic_table.find(symbol);
    if (it == periodic_table.end()) {
        HANDLE_ERROR ("Unknown atomic symbol: " + symbol_raw, 4);
    }
    return it->second;
}


#endif