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

#ifndef NEWSCF_UTILS_HPP
#define NEWSCF_UTILS_HPP

#include <newscf/cis.hpp>

#include <unordered_map>

namespace newscf::cis::utils {
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
                LOG (WARNING, "Unsupported angular momentum level "+std::to_string(ch));
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

#endif //NEWSCF_UTILS_HPP
