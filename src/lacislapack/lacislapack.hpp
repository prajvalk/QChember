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

#ifndef LACISLAPACK_HPP
#define LACISLAPACK_HPP

#ifndef LACIS_DISABLE_LAPACK_CHECKS
#define _LAPACK_INFO_CHECK(INFO) assert(INFO >= 0);
#elif
#define _LAPACK_INFO_CHECK(INFO)
#endif

#endif //LACISLAPACK_HPP
