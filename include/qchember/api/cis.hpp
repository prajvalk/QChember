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

#ifndef NEWSCF_COMMON_INTERFACE_SYSTEM_HPP
#define NEWSCF_COMMON_INTERFACE_SYSTEM_HPP

#include <iostream>
#include <string>
#include <cstring>
#include <random>
#include <vector>

namespace newscf::cis {

/* Logging API */

enum LOGGING_MODES           { INFO, WARNING };
enum DEVELOPER_LOGGING_MODES { DEV_DUMP, DEV_INFO, DEV_WARN, DEV_ERR };
enum MEMORY_LOGGING_MODES    { MEM_INFO, MEM_WARN, MEM_ERR };

inline void newscf_log (const LOGGING_MODES mode,
                        const std::string& message,
                        const std::string& file,
                        const int line) {
    if (mode == INFO)         std::cout << "(INFO "    << file << ":" << std::to_string(line) << ") " << message << std::endl;
    else if (mode == WARNING) std::cout << "(WARNING " << file << ":" << std::to_string(line) << ") " << message << std::endl;
}

#define LOG(MODE, MSG) \
        newscf_log(MODE, MSG, __FILE__, __LINE__)

#define LOG_ASIS(MODE, MSG, F, L) \
        newscf_log(MODE, MSG, F, L)

inline void newscf_devlog (const DEVELOPER_LOGGING_MODES mode,
                           const std::string& message,
                           const std::string& file,
                           const int line) {
    std::string prefix = "";
    if      (mode == DEV_DUMP) prefix = "DUMP";
    else if (mode == DEV_INFO) prefix = "INFO";
    else if (mode == DEV_WARN) prefix = "WARN";
    else if (mode == DEV_ERR)  prefix = "ERR ";

    std::cout << "[" << prefix << " " << file << ":" << line << "]: " << message << std::endl;
}

#define DEVLOG(MODE, MSG) \
        newscf_devlog (MODE, MSG, __FILE__, __LINE__)

#define DEVLOG_ASIS(MODE, MSG, F, L) \
        newscf_devlog (MODE, MSG, F, L)

inline void newscf_memlog (const MEMORY_LOGGING_MODES mode,
        const std::string& message,
        const std::string& file,
        const int line) {
        std::string prefix = "";
        if      (mode == MEM_INFO) prefix = "INFO";
        else if (mode == MEM_WARN) prefix = "WARN";
        else if (mode == MEM_ERR) prefix = "ERR ";
        std::cout << "[" << prefix << " " << file << ":" << line << "]: " << message << std::endl;
}
#ifndef DISABLE_MEMORY_LOGGING
#define MEMLOG(MODE, MSG) \
        newscf_memlog (MODE, MSG, __FILE__, __LINE__)

#define MEMLOG_ASIS(MODE, MSG, F, L) \
        newscf_memlog (MODE, MSG, F, L)
#else
#define MEMLOG(MODE, MSG)
#define MEMLOG_ASIS(MODE, MSG, F, L)
#endif

/* Memory API */

inline void* newscf_malloc (const size_t size) {
        return std::malloc(size);
}

#define MALLOC(SIZE) \
        newscf_malloc(SIZE); \
        MEMLOG(MEM_INFO, "Allocating "+std::to_string(SIZE)+" bytes");

inline void* newscf_memset (void* ptr, const int val, const size_t size) {
        return std::memset(ptr, val, size);
}

#define MEMSET(PTR, VAL, SIZE) \
        newscf_memset(PTR, VAL, SIZE); \
        MEMLOG(MEM_INFO, "Setting "+std::to_string(SIZE)+" bytes to "+std::to_string(VAL));

inline void newscf_memcpy (void* dest, const void* src, const size_t size) {
        std::memcpy(dest, src, size);
}

#define MEMCPY(DEST, SRC, SIZE) \
        newscf_memcpy(DEST, SRC, SIZE); \
        MEMLOG(MEM_INFO, "Copying "+std::to_string(SIZE)+" bytes");

inline void newscf_free (void* ptr) {
        std::free(ptr);
}

#define FREE(PTR) \
        newscf_free(PTR); \
        MEMLOG(MEM_INFO, "Freeing allocated memory");

inline void errorhandler (const std::string& message, const int code, const std::string& file, const int line) {
        DEVLOG_ASIS(DEV_ERR, message, file, line);
        throw std::runtime_error(message);
}

#define HANDLE_ERROR(MSG, CODE) \
        errorhandler (MSG, CODE, __FILE__, __LINE__);

/* CIS (Vector)/Matrix/Tensor System */

enum LACIS_TYPE {
        DENSE_LAPACK
};

enum LACIS_SAVE_TYPE {
        MATRIXMARKET,
        LZ4, CHNGLATER
};

 enum NormMetric { L1, L2 };

template <LACIS_TYPE BACKEND, typename T>
class LACIS {
public:
        void* data;
        size_t ndata;
        void* dims;
        size_t ndims;

        LACIS();

        ~LACIS();

        LACIS<BACKEND, T>& operator=(const LACIS<BACKEND, T>& rhs);

        void resizeToVector(size_t size);

        void resizeToMatrix(size_t r, size_t c);

        void resizeToMatrix(size_t i);

        void resizeToTensor3D(size_t i, size_t j, size_t k);

        void resizeToTensor3D(size_t i);

        void resizeToTensor4D(size_t i, size_t j, size_t k, size_t l);

        void resizeToTensor4D(size_t i);

        void fill(T t);

        T vectorGet (size_t i);

        void vectorSet (size_t i, T val);

        T matrixGet (size_t i, size_t j);

        void matrixSet (size_t i, size_t j, T val);

        T tensor3DGet (size_t i, size_t j, size_t k);

        void tensor3DSet (size_t i, size_t j, size_t k, T val);

        T tensor4DGet (size_t i, size_t j, size_t k, size_t l);

        void tensor4DSet (size_t i, size_t j, size_t k, size_t l, T val);

        template<LACIS_SAVE_TYPE>
        void save_to_file(const std::string& file);

        template<LACIS_SAVE_TYPE>
        void load_from_file(const std::string& file);

        template <NormMetric type>
        T compareTo(const LACIS<BACKEND, T>& lhs);

        void vectorPrint();
};

enum LinearSolverType {
        DSYSV, DSYSV_ROOK, DSYSV_RK,
        DSYSV_AA, DSYSV_AA_2
};

/*
 * Linear Solver
 * AX = B
 *
 * Supported backends: LAPACK/DSYSV
 */
template <LinearSolverType>
class LinearSolver {
public:
        /*
         * DSYSV: <double * LDA * N>[A] + <double * LDB * NRHS>[B] + <int * N>[IPIV] + <double * LWORK>[WORK]
         */
        void* work;

        /*
         * DSYSV: 4 Partitions (A | B | IPIV | WORK)
         */
        void* work_ptrs;

        /*
         * DSYSV: UPLO
         */
        char CH1;
        char CH2;

        /*
         * DSYSV: N
         */
        int  I1;

        /*
         * DSYSV: NRHS
         */
        int  I2;

        /*
         * DSYSV: LDA
         */
        int  I3;

        /*
         * DSYSV: LDB
         */
        int  I4;

        /*
         * DSYSV_AA_2: LTB
         */
        int  I5;

        int  LWORK;
        int  INFO;

        LinearSolver() {
                work = nullptr;
                work_ptrs = nullptr;
                CH1 = 'N';
                CH2 = 'N';
                I1 = 0;
                I2 = 0;
                I3 = 0;
                I4 = 0;
                I5 = -1;
                LWORK = -1;
                INFO = -1;
        }

        ~LinearSolver() {
                FREE(work);
                FREE(work_ptrs);
        }

        /* Set Symmetric Matrix Structuture
         * Used by: DSYSV
         */
        void set_uplo (char ch);

        /* Set Matrix A dimension */
        void set_n (int n);

        /* Set number of RHS to be solved */
        void set_nrhs (int nrhs);

        /* Set Matrix A, leading dimension */
        void set_lda (int lda);

        /* Set Matrix B, leading dimension */
        void set_ldb (int ldb);

        /* Initialized solver memory
         * Call after: set_n(), set_nrhs()
         */
        void init();

        /* Set Matrix A;
         * Call after: init()
         */
        template <LACIS_TYPE LT, typename LTtype>
        void set_a (const LACIS<LT, LTtype>& A);

        /* Set Matrix B;
         * Call after: init()
         */
        template <LACIS_TYPE LT, typename LTtype>
        void set_b (const LACIS<LT, LTtype>& B);

        /* Solve A x = {b1, b2, b3, ...}
         * Call after: init(), set_a(), set_b()
         */
        void solve();

        /* Copy answers {x1, x2, x3, ...}
         * Call after: solve()
         */
        template <LACIS_TYPE LT, typename LTtype>
        void copy_result (LACIS<LT, LTtype>& res);
};

enum GeneralizedEigenvalueSolverType {
        DSYGV, DSYGV_2
};

template <GeneralizedEigenvalueSolverType>
class GeneralizedEigenvalueSolver {
public:
        void* work;
        void* work_ptrs;
        char CH1; // JOBZ
        char CH2; // UPLO
        int  I1; // ITYPE
        int  I2; // N
        int  I3; // LDA
        int  I4; // LDB
        int  I5;
        int  LWORK1;
        int  LWORK2;
        int  INFO;

        GeneralizedEigenvalueSolver() {
                work = nullptr;
                work_ptrs = nullptr;
                CH1 = 'N';
                CH2 = 'N';
                I1 = 0;
                I2 = 0;
                I3 = 0;
                I4 = 0;
                I5 = 0;
                LWORK1 = -1;
                LWORK2 = -1;
                INFO = 0;
        }

        ~GeneralizedEigenvalueSolver() {
                FREE(work);
                FREE(work_ptrs);
        }

        void set_itype (int  itype);
        void set_jobz  (char jobz);
        void set_uplo  (char ch);
        void set_n     (int n);
        void set_lda   (int lda);
        void set_ldb   (int ldb);

        void init();

        template <LACIS_TYPE LT, typename LTtype>
        void set_a    (const LACIS<LT, LTtype>& A);

        template <LACIS_TYPE LT, typename LTtype>
        void set_b    (const LACIS<LT, LTtype>& B);

        void solve();

        template <LACIS_TYPE LT, typename LTtype>
        void copy_eigenvalues (LACIS<LT, LTtype>& res);

        template <LACIS_TYPE LT, typename LTtype>
        void copy_eigenvectors (LACIS<LT, LTtype>& res);
};

// LACIS-LAPACK Default Exports
typedef LACIS<DENSE_LAPACK, double> lacis_t;
typedef LinearSolver<DSYSV> linearsolver_t;

// Randomization API
/*
 * (Pseudo) (Batched Singular) Random Generator for Stack-like Use
 */
template<typename T, typename Engine = std::mt19937_64, typename Dist = std::uniform_real_distribution<T>>
class CISRandomGenerator {
public:
	CISRandomGenerator(size_t sz, T lower = T(0), T upper = T(1))
		: engine_(std::random_device{}()), dist_(lower, upper), workspace_(sz)
	{
		generateBatch();
	}

	const std::vector<T>& getBatch() const noexcept {
		return workspace_;
	}

	T next() {
		if (workptr_ >= workspace_.size()) {
			generateBatch();
			workptr_ = 0;
		}
		return workspace_[workptr_++];
	}

private:
	void generateBatch() {
		for (auto &val : workspace_)
			val = dist_(engine_);
	}

	Engine engine_;
	Dist dist_;
	std::vector<T> workspace_{};
	size_t workptr_ = 0;
};


/*
 * (Pseudo) (Non-Batched Multi-reference) Random Generator for Generator-like Use
 */
template<typename T, typename Engine = std::mt19937_64, typename Dist = std::uniform_real_distribution<T>>
class CISRandomGenerator2 {
private:
	Engine engine;
	std::vector<Dist> distributions; // One distribution per parameter
	size_t nparams;

public:
	CISRandomGenerator2(size_t params, const T* bounds)
		: engine(std::random_device{}()), nparams(params)
	{
		distributions.reserve(params);
		for (size_t i = 0; i < params; ++i) {
			T lb = bounds[2 * i];
			T ub = bounds[2 * i + 1];
			distributions.emplace_back(lb, ub);
		}
	}

	void next(T* result) {
		for (size_t i = 0; i < nparams; ++i) {
			result[i] = distributions[i](engine);
		}
	}
};



} // CIS Interface

#endif //NEWSCF_COMMON_INTERFACE_SYSTEM_HPP
