# Cholesky
<!-- TOC -->

- [Cholesky](#cholesky)
  - [Core library](#core-library)
  - [Linear algebra](#linear-algebra)
  - [Third party](#third-party)
  - [Tests](#tests)

<!-- /TOC -->
## Core library
Is a project independent library aimed at easing the build of new libraries .It is included with:

```
  #include "core/core.h"
  using namespace __core__;
```
For more information see [core/Readme.md](core/Readme.md).

## Linear algebra

Is the main set of functions aimed at computing the Cholesky decomposition. It is included with:

```
  #include "linear-algebra/linear-algebra.h"
  using namespace __core__;
```
For more information see [linear-algebra/Readme.md](linear-algebra/Readme.md).

## Third party

Is a set of wrappers to popular libraries like _Suite Sparse_ and _METIS_ needed in certain algorithms or test cases. It is included with:

```
  #include "third-party/third-party.h"
  using namespace __third_party__;
```
For more information see [third-party/Readme.md](third-party/Readme.md).

## Tests
Contains various test cases designed to verify the correctness of the results computed. In particular contains the tests:
- `int choleskyDenseTest(std::string filename)` which loads a sparse matrix and computes the Cholesky decomposition as if the matrix was dense.
- `int sparseCholesky(std::string filename,std::string outName)` it loads a sparse matrix, computes the sparse left looking Cholesky algorithm and stores the results in `outName`, storing the results of the original matrix `A`, the Cholesky decomposition `L` and the elimination tree `ETree` in `.csv` files.
- The test matrices which serve as input for the previous tests are located in the folder:
> test/matrices
