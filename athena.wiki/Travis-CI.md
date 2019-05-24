[Travis CI](https://travis-ci.com/) is a cloud-based, hosted
["Software as a Service" (SaaS)](https://en.wikipedia.org/wiki/Software_as_a_service) that provides continuous integration for projects hosted on GitHub. Travis CI offers free builds for open source projects in public repositories; the Athena++ private development repository receives a free, single concurrent build under their Education program. Each job in a build is currently limited to 120 minutes under this allocation.
<!-- When did they downgrade the Education benefit from "two free concurrent builds"? -->
<!-- The API and core is open-source; Travis CI used to fully open-source but turned into a company in 2012 -->
<!-- TODO: add instructions about 30 min interactive debug sessions? -->

[Customizing the Build](https://docs.travis-ci.com/user/customizing-the-build) is defined completely in the `.travis.yml` hidden file in the root directory of the repository. The file is written in [YAML](https://en.wikipedia.org/wiki/YAML), a data serialization language (similar to JSON) that is commonly used for configuration files. This file is largely used to define the [phases and jobs](https://docs.travis-ci.com/user/for-beginners/#Builds%2C-Jobs%2C-Stages-and-Phases) in the build lifecycle from `install`-ing dependencies, to calling the test `script`, to defining the follow-up actions for `after_success` and `after_failure`.
<!-- TODO: add Build Stages (Beta): https://docs.travis-ci.com/user/build-stages/ -->

The actual regression tests are called during the `script` phase in `tst/ci/travis/run_tests_travis.sh`, a Bash script executed with `set -e` so that the script exits immediately after the first error. Additional helper scripts for installing library dependencies FFTW, MPICH, and OpenMPI are located in `tst/ci/travis/`. All Travis CI files are version controlled with Athena++ so that the build configuration details are never lost and can evolve transparently with the source code.

Any user who has access to the private Athena++ repository can view the Travis CI builds directly on https://travis-ci.com. From [Travis CI: Who has access to the builds?](https://docs.travis-ci.com/user/travis-ci-for-private/#Who-has-access-to-the-builds%3F):
> Access rights on Travis CI is based on the access rights on GitHub:
>
> - Users that can access a repository on GitHub can see the build status and logs on Travis CI.
> - Users that can push to a repository on GitHub can trigger, cancel and restart builds, and change its settings.
> - Users that have admin access to a repository on GitHub can enable/disable it on Travis CI.

As of 5/18/2018, the Travis CI GitHub App has replaced the deprecated GitHub Service + Webhook setup in the private Athena++ repository. The new interface supports the [GitHub Checks API](https://blog.github.com/2018-05-07-introducing-checks-api/), so the build environment, test stages, and build results for pull requests are also summarized directly on GitHub. <!-- will this be extended to commits? -->
An example for Athena++ can be seen [here](https://github.com/PrincetonUniversity/athena/runs/991733). Future customization will improve the annotations of failed builds. For example, it might provide a one line summary if the tests failed during the C++ style check for maximum line length <=90 versus during the compilation of a particular file in `pgen/`.

### [Travis CI: Building Specific Branches](https://docs.travis-ci.com/user/customizing-the-build/#building-specific-branches)
- When a commit is pushed to the Athena++ repository, Travis CI checks to see if the branch that contains the commit is safelisted for testing. 
- The `branches:` block in the `.travis.yml` file specifies the safelisted (`only:`) and blocklisted (`except:`) lists of branches for testing. 
- These parameters can use regular expressions to specify multiple branches. 

### [Travis CI: Building Pull Requests](https://docs.travis-ci.com/user/pull-requests/)
- This setting is toggled by the **Build pushed pull requests** switch on the [Travis CI athena settings]( https://travis-ci.com/PrincetonUniversity/athena/settings)
- Travis will only build a pull request if it can be merged without a merge conflict. Be sure to resolve any conflicts or rebase against the upstream (typically `master`) branch.
- The pull request builder will always clone the latest version of a PR, even if that changes between jobs in a build. A discrepancy in versions could arise between jobs if, for example, a new commit on the target branch or base branch is pushed immediately after starting a long running build.
- The testing scripts can refer to the environment variable `$TRAVIS_PULL_REQUEST` to create conditional test steps. Similarly, `$TRAVIS_BRANCH` can be used, but it should not be relied upon in pull request builds.
- Currently, Travis CI is setup with **Auto Cancellation** of both branch and pull request builds. When new commits trigger builds in either queue, any builds waiting in the queue are cancelled.
