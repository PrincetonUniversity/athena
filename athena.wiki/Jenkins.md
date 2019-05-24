[Jenkins](https://en.wikipedia.org/wiki/Jenkins_%28software%29) is an open source automation server written in Java. Plugins are used to manage functionalities and interfaces with a variety of VCS and other tools. Princeton University's [Research Computing (RC)](https://researchcomputing.princeton.edu/) department operates a Jenkins server as a shared resource for computational projects and users on campus.

The Jenkins setup for Athena++ is configured as three Projects in the `athena/` [Folder](https://www.cloudbees.com/products/cloudbees-jenkins-platform/team-edition/features/folders-plugin):
- `PrincetonUniversity_athena_jenkins_commit` automatically tracks, builds, and tests any commits pushed to the `master` branch on PrincetonUniversity/athena
- `PrincetonUniversity_athena_jenkins_PR` automatically builds any pull requests targeting `master` opened by whitelisted authors on PrincetonUniversity/athena
- `PrincetonUniversity_athena_jenkins_branch` automatically tracks, builds, and tests any commits pushed to any of the configured feature branches on PrincetonUniversity/athena. 
  - As branches are created, merged into `master`, and deleted, this configuration must be manually updated with the names and refspecs of branches for which testing with Jenkins CI is desired.

As of June 2018, there is no GitHub App for Jenkins integration; the server primarily interacts with GitHub via [Webhooks](https://help.github.com/articles/about-webhooks/) that deliver messages triggered by certain activities on the repository. They are configured by Athena++ repository administrators under [Settings > Webhooks](https://github.com/PrincetonUniversity/athena/settings/hooks).

All project configurations at https://jenkins.princeton.edu/ are version controlled via a Git repository that the Jenkins server administrator (David Luet) controls. The Jenkins server merely orchestrates the automation of the build and testing cycle; the actual tests run on slave nodes that can be configured on most of the PICSciE clusters.
The build consists of executing the following command through the [Slurm Workload Manager](https://en.wikipedia.org/wiki/Slurm_Workload_Manager) on the worker computer:
```bash
#!/usr/bin/env bash
export BASH_ENV=$HOME/.bashenv
salloc -N1 -n4 --time=3:50:00 ./tst/ci/jenkins/run_jenkins_perseus.sh
```

Currently, Jenkins is configured to use the Perseus cluster as the exclusive build environment for Athena++; four independent Build Executors are allocated to the projects under `perseus_kfelker`, but they are subject to the same queuing on Perseus as normal Slurm jobs. The Build Executor utilization of the two projects can be tracked at [Load Statistics](https://jenkins.princeton.edu/computer/perseus_kfelker/load-statistics).

### [Jenkins pull request builder plugin](https://wiki.jenkins.io/display/JENKINS/GitHub+pull+request+builder+plugin)
The [@buildbot-princeton](https://github.com/buildbot-princeton) user is a bot on GitHub that interfaces with the Jenkins server at Princeton. It must always have Write access as a repository Collaborator in order to function properly. This bot may post comments on pull requests, and it will listen to replies in the comments from Jenkins project administrators and whitelisted users. The default special command syntax is:
* "ok to test" to accept this pull request for testing
* "test this please" to launch a one-time test run of the pull request
* "add to whitelist" to add the author to the whitelist (all future pull requests by this GitHub user will automatically be built)
* "retest this please" to immediately start a new build
