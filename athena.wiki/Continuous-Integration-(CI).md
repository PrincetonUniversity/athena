[Continuous integration](https://en.wikipedia.org/wiki/Continuous_integration) (CI) is the software development practice of frequently integrating small code changes into a central repository. CI best practices include the automation of the building, testing, and deploying the application. Apart from avoiding the "integration hell" and merge conflicts that can occur when infrequently adding large sets of changes, adopting CI offers several advantages:
- Bugs are detected quickly and fixed easily after they are introduced -> improving developer collaboration
- Lower testing costs -> saving developer time and effort performing manual testing
- Increased confidence in the stability and robustness of changes -> allowing developers to focus on writing higher quality code

Continuous integration may include many policies and tools, including the basic use of a version control system (VCS) / source control management (SCM) setup such as Git. Above all, it is a [cultural commitment from the developers](https://www.atlassian.com/continuous-delivery/how-to-get-to-continuous-integration) to commit changes frequently, fix broken code immediately, and write and maintain tests that cover all aspects of the project. See [Agile Alliance: Continuous Integration](https://www.agilealliance.org/glossary/continuous-integration) and [Travis CI: What is CI?][1] for more info.

Hereafter in this Wiki page, continuous integration will refer to the platforms and tools used for the **automated building and testing** of Athena++, also known as CI servers or build servers. While the initial setup and maintenance overhead of the build environment and automated testing suite are non-trivial, the costs are justified for collaborative software development on a project as large as Athena++.

### Setup summary
Continuous integration for the Athena++ private GitHub repository is currently provided by [[Travis CI]] and a [[Jenkins]] server. Both services independently check:
- All tests in the [[Regression Testing]] framework to check code functionality
- C++ and Python code <a href="https://en.wikipedia.org/wiki/Lint_(software)">linters</a> to confirm that conventions in the [[Style Guide]] are obeyed

In order to test as many potential user environments as possible, these test are repeated using a mix of compilers and library versions. All CI-specific files are located in `tst/ci/` folder and the `.travis.yml` hidden file in the root directory. All CI files are automatically excluded from the public version of the code when releases are created using the scripts in `pub/` directory (which is also excluded from public releases).

Builds are triggered on both platforms whenever:
1. A pull request targeting the `master` branch is opened or updated with new commits
2. Commits are pushed to `master`
3. A user manually triggers builds directly on the platforms (or by [comment commands](https://github.com/jenkinsci/ghprb-plugin/blob/master/README.md) in the Jenkins Pull Request builder)
<!-- the new GitHub Checks API also enables manually restarting builds for Travis CI on GitHub -->

Jenkins and Travis CI build statuses are reported through several channels:
* Pending, passing, or failing status is directly updated on the repository via the [GitHub status checks](https://help.github.com/articles/about-status-checks/) and the new Checks API tab (Travis CI builds of pull requests only)
* Finished test results are posted to the Slack workspace's `#development` channel
* Build failures (with attached logs) are emailed to certain CI server/repository administrators.

All of these settings are configurable.

#### Comparison of build environments
[[Jenkins]] and [[Travis CI]] serve complementary roles in the testing of Athena++:
- Jenkins runs tests on a single Linux cluster, but it has access to the extensive HPC libraries and tools that Princeton Research Computing offers, including [licensed software](https://researchcomputing.princeton.edu/software).
  - Therefore, the [[Regression Testing]] suite is run twice for each Jenkins build: once with open-source libraries and GCC, and once with Intel compilers and libraries.
- As a cloud-based service, Travis CI offers a myriad of [virtualization environments](https://docs.travis-ci.com/user/reference/overview/#Virtualization-environments), but these isolated environments are (mostly) built from scratch for each build. All dependencies must be re-installed in either a Docker container or a virtual machine (VM) image. User accounts are subject to resource and permissions constraints.
  - The Travis CI [Build Matrix](https://docs.travis-ci.com/user/customizing-the-build#Build-Matrix) tests both MPICH and OpenMPI each with both `g++` and `clang++` on Ubuntu containers. When combined with the 2x (MPICH, OpenMPI) + Apple `clang++` jobs, Travis CI executes a total of 6 independent jobs for each build.

The below table summarizes the rough division of build environments between the Travis CI and Jenkins services.

*Last updated on 6/8/2018.*

|                  |                  | Travis CI        | Jenkins          |
|------------------|------------------|------------------|------------------|
| **Hardware**     | Infrastructure   |  Docker container on Amazon EC2 (Linux) <br> VM image on Google GCE (macOS)| Compute node on PICSciE Perseus cluster |
|                  | CPU              | Unspecified x86_64             |  Intel Xeon CPU E5-2680 v4 @ 2.40GHz             |
|                  | Cores            | 2                | Dual-socket, 14 cores each |
|                  | Memory           |     4 GB         | 128 GB |
| **Operating systems** | Linux | Ubuntu 14.04 "Trusty" | Red Hat Enterprise Linux (RHEL) 7 |
|        | macOS | Xcode 9.3 macOS 10.13.2 High Sierra (MPICH) <br> Xcode 8.3 macOS 10.12 Sierra (OpenMPI) | - |
|        | Windows | -         | -
| **C++ compilers**| Intel `icc`     | - | 18.0.2 |
|                  | GNU `g++` | 4.8.4 (Ubuntu) | 6.3.1 |
|                  | `clang++` |  5.0.0 (Ubuntu) <br> Apple clang-902.0.39.1 (macOS) <br> Apple clang-802.0.38 (macOS) | - |
| **MPI libraries**    | Intel MPI | - | 2018 Release 2 (with `icc`) |
|                  | OpenMPI   | 3.0.2 | 1.10.2 (with `g++`)|
|                  | MPICH     | 3.2.1 | - |
| **Other libraries**  | Python  | 3.4.3 (Ubuntu) <br> 3.6 (macOS) | 2.7 |
|                  | FFTW | 3.3.7 | 3.3.4 |
|                  | HDF5 | 1.10.1  | 1.8.12 |
<!-- Maybe split the Travis CI column into 2x cols: Ubuntu Docker vs. macOS VM -->
<!-- Use "h5fc -showconfig" to get HDF5 library version on Perseus (not explicitly loading any HDF5 modules right now) -->

*Notes*:
- OpenMPI [changed their version numbering scheme](https://www.open-mpi.org/papers/versioning-update-2015/) in 2015 to something similar to Semantic Versioning (as Athena++ uses). OpenMPI 1.10.2 was released 1/21/2016, OpenMPI 2.1.1 was released 5/10/2017, and OpenMPI 3.0.2 was released 6/1/2018. See [OpenMPI Release Timeline](https://www.open-mpi.org/software/ompi/versions/timeline.php).
- FFTW and HDF5 libraries are built from the GNU Complier Collection (GCC) for all services and build jobs.
- After the first build, Travis CI caches the compiled libraries in order to reduce the startup time of the jobs. The cache can be manually deleted at any time.
  - Installed with Homebrew for macOS builds
  - Installed from source with GCC for Ubuntu builds
- Jenkins uses the [Environment Modules](http://modules.sourceforge.net/) package installed on the Linux clusters to load the requisite build environment libraries.
- Different Python libraries are used primarily to certify that the [[Regression Testing]] scripts remain compatible with both Python 3 and 2, especially regarding `print()` function usage, string/byte encoding, and list vs. iterator behavior.

[1]: https://docs.travis-ci.com/user/for-beginners/#What-is-Continuous-Integration-(CI)%3F
