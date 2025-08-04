# Common Interfaces System (CIS) Reference Specification

Specifications marked with (*) are already implemented by the interface header and cannot be specialized by vendor modules.

_Specification Version: 1.0-alpha-aug4-2025_

## Logging Interface

* ```LOG(MODE, MSG)``` (*)
* ```LOG_ASIS(MODE, MSG, FILE, LINE)``` (*)
* ```DEVLOG(MODE, MSG)``` (*)
* ```DEVLOG(MODE, MSG, FILE, LINE)``` (*)
* ```MEMLOG(MODE, MSG)``` (*)
* ```MEMLOG(MODE, MSG, FILE, LINE)``` (*)

Supported normal logging modes: ```INFO, WARNING``` <br>
Supported developer logging modes: ```DEV_DUMP, DEV_INFO, DEV_WARN, DEV_ERR``` <br>
Supported memory logging modes: ```MEM_INFO, MEM_WARN, MEM_ERR```

## Memory Interface

* ```MALLOC(SIZE)``` (*)
* ```MEMSET(PTR, VAL, SIZE)``` (*)
* ```MEMCPY(DEST, SRC, SIZE)``` (*)
* ```FREE(PTR)``` (*)

## Error Handling

* ```HANDLE_ERROR(MSG, CODE)``` (*)

## Linear Algebra System

Linear Algebra Common Interfaces System (LACIS):

Supported Backends: ```DENSE_LAPACK``` <br>
Supported Save Types: ```MATRIXMARKET, LZ4``` <br>
Supported Norm Calculation Modes: ```L1, L2``` <br>

### ```LACIS<BACKEND, T>```
* ```void* data```
* ```void* dims```
* ```size_t ndata```
* ```size_t ndims```
* ```LACIS()```
* ```~LACIS()```
* ```LACIS<BACKEND, T>& operator= (const LACIS<BACKEND, T>&)```

Implements equidimensional data structures: 

* ```void resizeToVector(size_t)```
* ```void resizeToMatrix(size_t)```
* ```void resizeToTensor3D(size_t)```
* ```void resizeToTensor4D(size_t)```

Implements generic data structures:

* ```void resizeToVector(size_t)```
* ```void resizeToMatrix(size_t, size_t)```
* ```void resizeToTensor3D(size_t, size_t, size_t)```
* ```void resizeToTensor4D(size_t, size_t, size_t, size_t)```

Get/Set Operations:

* ```T vectorGet(size_t)```
* ```T matrixGet(size_t, size_t)```
* ```T tensor3DGet(size_t, size_t, size_t)```
* ```T tensor4DGet(size_t, size_t, size_t, size_t)```
* ```void vectorSet(size_t, T)```
* ```void matrixSet(size_t, size_t, T)```
* ```void tensor3DSet(size_t, size_t, size_t, T)```
* ```void tensor4DSet(size_t, size_t, size_t, size_t, T)```

Utility Functions:

* ```void fill(T)```
* ```T compareTo<NormMetric>(const LACIS<BACKEND, T>&)```
* ```void vectorPrint()``` 

IO Functions

* ```void save_to_file<LACIS_SAVE_TYPE>(const std::string&)```
* ```void load_from_file<LACIS_SAVE_TYPE>(const std::string&)```

### LinearSolver

Supported Solver Types: ```DSYSV, DSYSV_ROOK, DSYSV_RK, DSYSV_AA, DSYSV_AA_2```

* ```void* work```
* ```void* work_ptrs```
* ```char CH1, CH2```
* ```int I1, I2, I3, I4, I5```
* ```int LWORK, INFO```
* ```LinearSolver()``` (*)
* ```~LinearSolver()``` (*)

See ```cis.hpp``` for full interface.

### GeneralizedEigenvalueSolver

Supported Solver Types: ```DSYGV```

See ```cis.hpp``` for full interface.

## Convenience Export Types

* ```typedef LACIS<DENSE_LAPACK, double> lacis_t;```
* ```typedef LinearSolver<DSYSV> linearsolver_t;```