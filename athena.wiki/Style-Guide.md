Conforming to a set of style guidelines is crucial for maintaining code clarity.  A consistent style makes it easier for others to understand the implementation when reading the source code, and it may prevent bugs from being introduced into the codebase. Style is very personal, and there is no one "correct" style for all codes.  However, in Athena++, all developers should make a best faith effort to conform to the following style guidelines:

* [[C++ Style]]
* [[Python Style]]

Athena++ is designed to be used with <a href="https://en.wikipedia.org/wiki/Lint_(software)">linters</a> for static analysis of the source code to catch style issues, syntax errors, and questionable usage of language constructs. User instructions for operating the recommended C++ and Python linter scripts are provided in each subpage. However, most editors can be configured to automatically invoke these linters.

<!-- Add sections on: 1) commit messages 2) Git hooks -->

#### Miscellaneous files
The Athena++ repository also includes Bash scripts, `athinput.` input files, Markdown documentation, YAML configurations, and other miscellaneous file formats. No specific style guidelines are enforced for these files, but it is recommended that users consult the style of existing examples of each file format in the repository before working with them.

 <!-- TODO: add Git hooks for checking style -->
