# Deltares
This is an adapted and stripped down version of f90tw. The original library did not work on windows, because the powershell commands failed. Furthermore, we only use googletest, so boost test was removed as an option.
The CMake code was simplified so that it is simpler for the user of the library.

# f90tw

f90tw project provides fortran wrappers for a limited subset of the [boost](https://www.boost.org/) and [google](https://github.com/google/googletest) test frameworks functionality. At the same time, offers a rather simple mechanism for setting up and managing test suites in fortran projects. All the way down, f90tw implementation follows a preprocessor-based approach. The motivation was to utilize already available test frameworks and assertions implementation with minimal effort.

## Build

You will need to obtain the repository from github and follow the usual cmake build steps. It is assumed that a) c/c++ and fortran compilers are available and, b) you have installed boost and gtest packages. In short, [clone](https://github.com/git-guides/git-clone) the project from github, enter the project's directory and issue the usual [cmake](https://cmake.org/documentation/) build commands. On linux boxes, the flow of commands looks like:

```bash
% git clone https://github.com/loudove/f90tw
% cd f90tw
% mkdir build
% cd build
% cmake -DBOOSTFWK=OFF -DGTESTFWK=ON ..
% make
```

Using the cmake options `BOOSTFWK` and `GTESTFWK` you can select the framework to be used. `BOOSTFWK` is `ON` by default. Moreover, options `TESTS` and `EXAMPLES` controls the build of f90tw tests and examples.

There is not specific installation directory. After a successful build, the module files and the libraries will be located under the `<build>/f90tw/` directory (or `<build>/f90tw/<conf>` on windows). The actual source files of the modules will be located in the `<build>/f90tw` directory with the names `assertions_boost.f90` and `assertion_gtest.f90` for the boost and gtest frameworks respectively. You can use [fprettify](https://github.com/pseewald/fprettify) to format them properly in order to inspect them.

If `TESTS` or `EXAMPLES` option are `ON`, the corresponding target executable will be placed in corresponding build folder (`<build>/tests/<conf>` and `<build>/examples/<conf>` respectively).

## Use

A fortran module with the interface of the fortran wrappers and the corresponding libraries with the actual functionality are provided for each test framework:

- Boost: **assertions_boost** module (implementable in [`assertions_boost.fpp`](./f90tw/assertions_boost.fpp)), `libf90tw_boost.a` library (implemented in [`boost_assertions.cpp`](./f90tw/boost_assertions.cpp)). After a successful build a "readable" code for the module can be found in `<build>/f90tw/assertions_boost.f90`.
- Gtest: **assertions_gtest** module (implemented in [`assertions_gtest.fpp`](./f90tw/assertions_gtest.fpp)), `libf90tw_gtest.a` library (implemented in `gtest_assertions.cpp`). Also in this case, after a successful build a "readable" code for the module will be located in `<build>/f90tw/ftest_boost.f90`.

If you include f90tw in your cmake project, the variables `f90tw_INCLUDE_DIR` (pointing to the heads directory) and `f90tw_MODULE_DIR` (pointing to the modules and libraries directory) are defined to help you access f90tw functionality. You can find hinds of f90tw use with cmake in [`examples/CMakeLists.txt`](./examples/CMakeLists.txt). The variables `FPP_EXE` (preprocessor), `FPP_ID` (preprocessor id), `FPP_OPTIONS` (preprocessor options) and `FPP_SWITCH` (preprocessor options's switch) are also exposed allowing improved flexibility.

To facilitate the setup of simple tests, additional functionality is provided and demonstrated in the accompanying examples found in the [`example`](./examples) directory. A detailed walk-through of the implementation of the examples can be found [here](./examples/README.md).

Assuming that you want to test the functionality of a fortran module named *example*, implemented in [`example.f90`](./examples/example.f90). The test can be implemented with a pair of files, using either preprocessor macros or directives.

- Boost: this is an example of preprocessor macros approach. A fortran module, ([`test_example_boost.fpp`](./examples/test_example_boost.fpp)), with the subroutines implementing the actual tests and its c/c++ pair, [`test_example_boost.cpp`](./examples/test_example_boost.cpp), implementing the boost tests wrapping their fortran counterparts.
- Gtest: this is an example of process directives approach.Also in this case, a fortran module [`test_example_gtest.fpp`](./examples/test_example_gtest.fpp) implements the tests and [`test_example_gtest.cpp`](./examples/test_example_gtest.cpp) their c/c++ wrappers. **PLEASE NOTE** that these files should be compiled/preprocessed by setting the USE GTEST preprocessor definition in order to avoid problems.

### Preprocessor macros

In the fortran code file, we first include the header [`"f90tw_test.h"`](./f90tw/f90tw_test.h) and then, we implemented the test module with the help of f90tw preprocessor macros:

- **TESTMODULE**( MODULENAME, MODULETOTEST): (*optional*) prepares the initial statements of the module(module declaration, use of the module to be tested and set `implicit none`).
  - MODULENAME: the name of the module
  - MODULETOTEST: the name of the module to test

- **TESTCONTAINS**: (*optional*) expands to fortran `contain` statement. It is
separated from TESTMODULE macro to allow for the declaration of module resources common to all tests.

- **TESTCODE**( *TESTTYPE*, *TESTSUITENAME*, *TESTNAME*, *SUBNAME*, **...** ): prepares the specific test.
  - *TESTTYPE*: the framework macro to be used for test declaration
  - *TESTSUITENAME*: the name of the test suite
  - *TESTNAME*: the name of the test
  - *SUBNAME*: the name of the fortran subroutine to be implemented
  - **...** : the rest of the arguments which are essentially, the fortran source code implementing the test.
    **PLEASE NOTE** that each code line should be terminated with `";"`. The use of pound character (`"#"`) in fortran string renders preprocessing impossible, while, in order to avoid similar catastrophic errors, the string concatenation operator (`//`) should be replaced with the `F90CONCAT` macro and the line continuation operator (`&`) with the `F90LNCONT` macro. Moreover, `F90SPOT` macro expands to the `filename:line` string, while `F90SPOTMSG(MSG)` macro appends the MSG string to it. You can simplify the implementation by just calling a fortran method implemented elsewhere. This approach is maybe preferable since you will avoid the drawbacks of preprocessing relevant to code's clarity and debugging.

- **ENDTESTMODULE**(MODULENAME)  : (*optional*) module end statement.
  - MODULENAME: the name of the module (the same with the one used with with **TESTMODULE** macro).

Using this approach, the c++ implementation becomes rather fast and easy since it is based on preprocessing the same fortran file, with different definitions of the macros resulting automatically in a) the declarations of the fortran test methods and b) the framework tests which essentially wrap these methods. A more detailed description of the fortran and c/c++ files is provided [here](./examples/README.md).

### Directives

An alternative approach is to implement the test using standard fortran (see [`test_example_gtest.f90`](./examples/test_example_gtest.f90)) and the `!$f90tw` directive (`[!cC]\$[fF]90[tT][wW]]` in regular expression) to indicate the test cases i.e. the subroutines to be wrapped. Essentially, the lines with the directives will be extracted and included in the c/c++ counterpart ((see [`test_example_gtest.cpp`](./examples/test_example_gtest.cpp))) in order to create the test suite. A more detailed description of the fortran and c/c++ files is provided [here](examples/README.md).

## Assertions tests

In addition to the assertion wrappers for [boost.test](./README.md#boost.test-assertions-tests) and [gtest](./README.md#gtest-assertions-tests), a method for accessing f90tw version is available:

- f90tw_version(*major*, *minor*, *patch*)
  - *major* : major version
  - *minor* : minor version
  - *patch* : patch version

### Boost.test assertions tests

The following subset of boost.test framework is supported with &lt;level&gt; = (WARN&#124;CHECK&#124;REQUIRE) :

Fortran (F90_&lt;C/C++&gt;) | operator | argumens type
--------------------------- | -------- | -------------
F90_BOOST_&lt;level&gt;_&lt;operator&gt; | MESSAGE <sup>(*)</sup> | logical(KIND=(C_BOOL&#124;4)),<br> integer(KIND=C_INT)
F90_BOOST_&lt;level&gt;_&lt;operator&gt; | BITWISE_EQUAL <sup>(*)</sup> | integer(KIND=C_INT)
F90_BOOST_&lt;level&gt;_&lt;operator&gt; | (EQUAL&#124;NE) | logical(KIND=(C_BOOL&#124;4)),<br> integer(KIND=C_INT),<br> real(KIND=C_FLOAT),<br> real(KIND=C_DOUBLE)
F90_BOOST_&lt;level&gt;_&lt;operator&gt; | (CLOSE&#124;CLOSE_FRACTION) <sup>(*)</sup> | real(KIND=C_FLOAT),<br> real(KIND=C_DOUBLE)
F90_BOOST_&lt;level&gt;_&lt;operator&gt; | SMALL | real(KIND=C_FLOAT),<br> real(KIND=C_DOUBLE)
F90_BOOST_&lt;level&gt;_&lt;operator&gt; | (GE&#124;GT&#124;LE&#124;LT) | integer(KIND=C_INT),<br> real(KIND=C_FLOAT),<br> real(KIND=C_DOUBLE)
F90_BOOST_TEST_MESSAGE | - | character(KIND=C_CHAR,LEN=*)
F90_BOOST_TEST_ERROR | - | character(KIND=C_CHAR,LEN=*)
F90_BOOST_TEST_FAIL | - | character(KIND=C_CHAR,LEN=*)
F90_BOOST_TEST_CHECKPOINT | - | character(KIND=C_CHAR,LEN=*)
F90_BOOST_TEST_INFO | - | character(KIND=C_CHAR,LEN=*)
F90_BOOST_TEST_INFO_SCOPE | - | character(KIND=C_CHAR,LEN=*)
F90_BOOST_TEST_CONTEXT | - | character(KIND=C_CHAR,LEN=*)
F90_BOOST_TEST_PASSPOINT | - | character(KIND=C_CHAR,LEN=*)

For more details please check the [Boost.Test](https://www.boost.org/doc/) documentation.

### Gtest assertions tests

The following subset of gtest framework is supported with &lt;level&gt; : (ASSERT&#124;EXPECT) :

Fortran (F90_&lt;C/C++&gt;) | operator | argumens type
--------------------------- | -------- | -------------
F90_&lt;level&gt;_&lt;operator&gt; |  (EQ&#124;NE) | logical(KIND=(C_BOOL&#124;4)),<br> integer(KIND=C_INT),<br> real(KIND=C_FLOAT),<br> real(KIND=C_DOUBLE)
F90_&lt;level&gt;_&lt;operator&gt; |  (GT&#124;GE&#124;LT&#124;LE) | integer(KIND=C_INT),<br> real(KIND=C_FLOAT),<br> real(KIND=C_DOUBLE)
F90_&lt;level&gt;_FLOAT_EQ <sup>(*)</sup> |  - | real(KIND=C_FLOAT)
F90_&lt;level&gt;_DOUBLE_EQ <sup>(*)</sup> |  - | real(KIND=C_DOUBLE)
F90_&lt;level&gt;_NEAR <sup>(*)</sup> |  - | real(KIND=C_FLOAT),<br> real(KIND=C_DOUBLE)
F90_&lt;level&gt;_&lt;operator&gt; |  (TRUE&#124;FALSE) <sup>(*)</sup> | logical(KIND=(C_BOOL&#124;4))
F90_&lt;level&gt;_&lt;operator&gt; |  (STREQ&#124;STRNE&#124;<br>STRCASEEQ&#124;STRCASENE) | character(KIND=C_CHAR,LEN=*)

In the gtest case, all the assertion methods can be used with a message as the last argument (of type character(KIND=C_CHAR,LEN=*) ) in order to overwrite the default assertion message.

For more details please check the [gtest](https://github.com/google/googletest) documentation.
