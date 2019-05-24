### General C++ style guidelines
The majority of the Athena++ source code consists of about 63,000 lines of C++ (*as of 6/11/18*). The most important guidelines to follow when editing or adding to this codebase are:
* All code should conform to the C++11 standard, and it must compile with `g++`, `icpc`, and `clang++` **without warnings** (sets of warning flags defined for each compiler in `tst/ci/set_warning_cflag.sh`). Avoid all language extensions.

* All code must be thoroughly commented. Use `//` (new/C++ style) (not `/* ... */` old/C style) for comments.  All `.cpp` files must begin with license boilerplate (as found in, e.g., `main.cpp`). <!-- Final implementation should include [doxygen](http://www.stack.nl/~dimitri/doxygen/index.html) annotation of classes and functions. -->

* Indent by 2 space characters, limit maximum line length to 90 characters, and do not use tabs for indenting.

* Use a 1 space indent for `public:`, `private:`, `protected:` class access specifiers in class definitions.

* Wrap long lines (arithmetic expressions, function declarations, function calls, conditionals, etc.) as follows:
  * If you are breaking the line "in the middle" of the RHS expression of the `=`, function argument/parameter list, etc., then align subsequent lines with the first argument / opening parenthesis. 
  * Or, break immediately after `=` or opening parenthesis, and start a new line with a +4 space indent (aligning subsequent lines with this +4 space indent).
  * See Google C++ Style Guide sections on [Function Calls](https://google.github.io/styleguide/cppguide.html#Function_Calls) and [Function Declarations and Definitions](https://google.github.io/styleguide/cppguide.html#Function_Declarations_and_Definitions) for more details and examples. 

* Strongly avoid use of C++ macros to define inline functions.  Avoid use of `#ifdef...#endif` blocks if at all possible.

* Use header (`.hpp`) files to define classes and function prototypes.  Function implementation should be in corresponding source code (`.cpp`) file.  Exceptions may be necessary (e.g. `athena_array.hpp`).

* File names: all lower case, words separated by underscore, e.g. `parameter_input.cpp`

* Type names (classes, structs, typedefs, enums, ...) ([PascalCase](https://en.wikipedia.org/wiki/Camel_case)): capitalize each word and avoid use of underscores, e.g. `InputParameters`

* Use [fixed-/precise-width integer types (available since C++11)](http://en.cppreference.com/w/cpp/types/integer) when 16-bit (`std::int16_t`, `std::int16_t`) or 64-bit (`std::uint64_t`, `std::int64_t`) integers are desired. Do not rely upon the [compiler and architecture-dependent](https://google.github.io/styleguide/cppguide.html#Integer_Types) `long` and `long long` types. See [C++ types](http://en.cppreference.com/w/cpp/language/types) for a table summarizing the variable width in bits of these types for several data models (LP32, ILP32, LLP64, LP64).
  * Use `ll` or `LL` literal suffix for `long long` (shorthand for `long long int`) for signed integer literals that are supposed to be at least 64-bits wide (e.g. for operations involving `std::int64_t` variables). Unlike `long int`, `long long int` is guaranteed to be at least 64 bits in size since the type was specified in the C99 standard. See [Wikipedia: C data types](https://en.wikipedia.org/wiki/C_data_types)
  * Similarly, use `ull` or `ULL` literal suffix for `unsigned long long int` for unsigned integer literals that must be compatible with `std::uint64_t` variables. 

* Only use unsigned integer types to store and manipulate bit patterns and `std::size_t` to store the result of the `sizeof` operator for portability. Do not mixed unsigned integer and signed integer types in arithmetic, since   the [usual arithmetic conversions](https://en.cppreference.com/w/cpp/language/operator_arithmetic) implicit rules may lead to challenging bugs. 

* [Do not define implicit conversions](https://google.github.io/styleguide/cppguide.html#Implicit_Conversions). Use explicit type conversions, especially from `int64_t` to (possibly) lower-precision `int` types. *Note:* only use these unsafe "narrowing" conversions when it is unavoidable, e.g. converting `std::size_t` unsigned integer type to signed `int` for MPI library calls. 

* Function names ([PascalCase](https://en.wikipedia.org/wiki/Camel_case)): treat like type names; use capitals for each word and avoid underscores, e.g. `MyFunction`.

* Variable names ([snake_case](https://en.wikipedia.org/wiki/Snake_case)): all lower case, words separated by underscore, e.g. `iso_csound`.

* Performance and clarity are the top priorities, with performance being the most important. No feature of C++ shall be used that sacrifices performance compared to an implementation in C. If clarity must be sacrificed for performance, the corresponding code segments must be thoroughly commented.

* Don't rely on external libraries (not even Boost, Blitz++, or I/O libraries), except when absolutely necessary (e.g. OpenMP, MPI, FFTW, HDF5).

* For each file, include **only** the C++ Standard Library headers that are direct dependencies for the contents of that file. For the C Standard Library subset of the library, include the equivalent C++ Standard Library name, e.g. `#include <cstring>` not the deprecated `#include <string.h>`. In accordance with the ISO C++11 standard, any references to identifiers from these headers should be qualified with the `std::` namespace, even though the library implementation may also (and usually does) place the symbols in the global namespace.
  * This rule excludes POSIX C extensions to the ISO C library, which are not in the C++ Standard Library: `signal.h, sys/stat.h, unistd.h` 

* Use [`nullptr` (available since C++11)](https://en.cppreference.com/w/cpp/language/nullptr), not `NULL`, to avoid rare issues whereby the latter is implicitly converted to an integral type

* Use [`override`](https://en.cppreference.com/w/cpp/language/override) or [`final`](https://en.cppreference.com/w/cpp/language/final) virtual specifier (both available since C++11) in class member function declarations when overriding a virtual function from a base class.

* Don't use `if` and `else if` without an `else` clause if it initializes a variable, even if the `else` branch is guaranteed to never execute at the time of writing. Such usages may cause bugs in the future when the guarantee is broken; they are typically caught by the `-Wsometimes-uninitialized` (Clang) or `-Wmaybe-uninitialized` (GCC)  compiler warning flags.

*Athena++ API-specific recommendations:*
* Use `il, iu`, `jl, ju`, `kl, ku` for local limits for loops along the x1, x2, x3 indices, respectively, to avoid confusion with the `is, ie`, `js, je`, `ks, ke` members of each `MeshBlock` class instance. The latter set of indices refer to the lower and upper boundaries of the real (non-ghost) zones in each direction. The former set of indices may extend into the ghost zones depending on the particular computation or stencil.

### Google C++ Style Guide
For topics not specified in this guide, follow the recommendations of the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html). In addition to the above guidelines, the following list contains some of the rules enforced in Athena++, organized by their [Cpplint](https://en.wikipedia.org/wiki/Cpplint) error category/code and severity level 1-5:
- `[build/]`
  - Use the following header include order: C headers, C++ headers, then project headers (exclude OpenMP, MPI, FFTW headers from this rule) `[build/include_order] [4]`
  - Alphabetizing header include statements within their category (C, C++, project headers) is required `[build/include_alpha] [4]`
  - All header files should have `#define` guards to prevent multiple inclusion. The format of the symbol name should be `<PROJECT>_<PATH>_<FILE>_HPP_`. To guarantee uniqueness, they should be based on the full path in a project's source tree. And, the `#endif` line should be `#endif  //<PROJECT>_<PATH>_<FILE>_HPP_` `[build/header_guard] [5]`
  - Do not use namespace using-directives.  Use using-declarations instead.  `[build/namespaces] [5]`
  - "Include what you use", e.g. always have `#include <cmath> ` in the file if you are using `std::sqrt()`, or `#include <string>` for string type `[build/include_what_you_use] [4]`
- `[legal/]`
  - You should have a line: "Copyright..." `[legal/copyright] [5]`
- `[readability/]`
  - You don't need a `;` after a `}`  `[readability/braces]` `[4]`
  - If an else has a brace on one side, it should have it on both  `[readability/braces] [5]`
  - Use C++-style casts like `static_cast<float>(double_value)`; do not use C-style cast formats like `int y = (int)x` `[readability/casting] [4]`
- `[runtime/]`
  - Use explicit keyword for one-parameter class constructors, so that they aren't confused with implicit conversions. `[runtime/explicit] [5]`
  - Never use `sprintf`. Use `snprintf` instead. `[runtime/printf] [5]` (*Note, the latter is safer if the destination buffer is too small*)
  - Static/global string variables are not permitted.  `[runtime/string] [4]`
- `[whitespace/]`
  - `{` should almost always be at the end of the previous line  `[whitespace/braces] [4]`
  - Use space before `else` keyword `[whitespace/braces] [5]`
  - If/else bodies with multiple statements require braces `[whitespace/braces] [5]`
  - Use space before and after inline comment starting characters ` // `, and at least two spaces is best between code and comments `[whitespace/comments] [4, 2]`
  - No trailing whitespace allowed `[whitespace/end_of_line] [4]`
  - Else clause should never be on same line as `else` keyword `[whitespace/newline] [4]`
  - Redundant blank line at the start (or end) of a code block should be deleted `[whitespace/blank_line] [2 (or 3)]`
  - No extra space before last semicolon on the line `[whitespace/semicolon] [5]`

The [following examples from Google Style Guide](http://google.github.io/styleguide/cppguide.html#Conditionals) are used to illustrate the proper use of curly braces, and horizontal and vertical whitespace around conditionals:
```c++
if (condition) { }  // Good - proper space after IF and before {.

if (condition) {  // no spaces inside parentheses
  ...  // 2 space indent.
} else if (...) {  // The else goes on the same line as the closing brace.
  ...
} else {
  ...
}
```
Do not use:
```c++
if(condition) {  } // Bad - space missing after IF.
if (condition){  } // Bad - space missing before {.
if(condition){  }   // Doubly bad.

// Not allowed - IF statement on one line when there is an ELSE clause
if (x) DoThis();
else DoThat();

// Not allowed - curly on IF but not ELSE
if (condition) {
  foo;
} else
  bar;

// Not allowed - curly on ELSE but not IF
if (condition)
  foo;
else {
  bar;
}

if (condition) // Bad-- { should always be at the end of the previous line
{
}
```

Paraphrasing Google's guide,
> The benefit of a style rule must be large enough to justify asking everyone to remember it.

In that spirit, the following are some of the Google C++ Style Guide rules that are explicitly **suppressed/ignored**:
- Use space after "," such as `f(i, j)` not `f(i,j)` `[whitespace/commas]`
- Use space after "=" unless in if-statement, for-statement, etc., e.g. `a = 2;` not `a=2;` `[whitespace/operators]` (*For readability, using whitespace around most operators is recommended*)
- All parameters passed by reference must be labeled const. `[runtime/references]`
- Small and focused functions are preferred `[readability/fn_size]`
  - This error is ignored only because of the developer effort necessary to correct these violations; reducing the size of these functions is a high priority.
  - As of writing, there are 4x functions that exceed this error's relative length threshold.
  - As of writing, Athena++ contains 10x files with over 1000 lines of code.

For the most up to date C++ style rules for Athena++, refer to the `CPPLINT.cfg` file in the root project directory.

*Last updated 1/22/19*.

### Checking C++ style
The [[Continuous Integration (CI)]] setup automatically checks changes to `master` to ensure that the style of C++ files remains consistent. Before executing the [[Regression Testing]] suite to validate the correctness of the solver, the testing pipeline uses [`cpplint.py`](https://en.wikipedia.org/wiki/Cpplint), an open source script developed by Google to identify style errors. The script is distributed [on GitHub](https://github.com/google/styleguide/tree/gh-pages/cpplint), and a copy of this script is stored in the Athena++ repository's `tst/style/` directory for convenience.

While `cpplint.py` accepts a variety of option flags at the command line, it also automatically searches the current, parent, and child directories for files named `CPPLINT.cfg` that contain configuration options for the style linter. One such file is version controlled in the Athena++ root directory and is primarily used to encode the filters for suppressing/ignoring style rules; any `CPPLINT.cfg` in subdirectories will override these project settings for the nested folder's contents.  

A similar setup is used to check [[Python Style]] in the CI framework.

The `cpplint_athena.sh` Bash script, also stored in the `tst/style/` directory, automatically pipes the appropriate C++ files in `src/` to both `cpplint.py` and additional style checks created specifically for Athena++ source code management, including:
* All `std::sqrt()` and `std::cbrt()` function calls reside in the `std::`, not global, namespace
* No tab characters are used

Note, `cpplint_athena.sh` skips files in `src/fft/plimpton/`, containing [Plimpton's library](http://www.sandia.gov/~sjplimp/docs/fft/README.html) source code packaged with Athena++ for parallel [[FFT]] capabilities.

Users are encouraged to use the linter wrapper script to check their modifications for any C++ style violations before sharing their local changes. The script must be invoked from the directory in which it resides:
```
cd tst/style/
./cpplint_athena.sh
```

Special comments can be added to the end of any line in the `.cpp`, `.hpp` files:
```
// NOLINT
// NOLINT(rule)
```
to ignore all style rules or just the selected `rule` category, respectively, when linting that particular line.
