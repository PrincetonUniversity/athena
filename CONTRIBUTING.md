# Introduction
Welcome! Thank you for considering contributing to Athena++. This project adheres to a [code of conduct](CODE_OF_CONDUCT.md), which you are expected to uphold by participating. 

The guidelines in this document are meant to help make the development of Athena++ straightforward and effective. They are a set of best practices, not strict rules, and this document may be modified at any time. Navigating the code can be daunting for new users, so if anything is unclear, please let us know!

<!-- ### Table of Contents -->


## Resources and quick links
* The latest development version of Athena++ is hosted in the public [PrincetonUniversity/athena](https://github.com/PrincetonUniversity/athena) GitHub repository.
* The public version of Athena++ was formerly distributed via [PrincetonUnviersity/athena-public-version](https://github.com/PrincetonUniversity/athena-public-version) GitHub repository. This repository is now archived. 
* The final version of the predecessor C code, Athena 4.2, has been migrated from its original SVN repository to GitHub at [PrincetonUniversity/Athena-Cversion](https://github.com/PrincetonUniversity/Athena-Cversion).
  * [Athena 4.2 Documentation](https://princetonuniversity.github.io/Athena-Cversion/AthenaDocs) is hosted on GitHub Pages.
  * The [Athena 4.2 Test Page](https://www.astro.princeton.edu/~jstone/Athena/tests/) contains useful algorithm test results.
* The [Athena++ website](https://www.athena-astro.app/) is hosted by [GitHub Pages](https://pages.github.com/).
* The latest version of the [Athena++ documentation](https://github.com/PrincetonUniversity/athena/wiki) is hosted as a GitHub Wiki attached to this repository.
* The Athena++ Slack workspace is located at [athena-pp.slack.com](https:://athena-pp.slack.com).
<!-- Could add links to New PR, New Issue, Issue Labels e.g. current "bugs" -->


# How to contribute
There are many ways to contribute! We welcome feedback, [documentation](#documentation), tutorials, scripts, [bug reports](#bug-reports), [feature requests](#suggesting-enhancements), and [quality pull requests](#pull-requests).

## Using the issue tracker
Both [bug reports](#bug-reports) and [feature requests](#suggesting-enhancements) should use the [GitHub issue tracker](https://github.com/PrincetonUniversity/athena/issues).

Please do not file an issue to ask a question on code usage.

### Bug reports
[Open a new Issue](https://github.com/PrincetonUniversity/athena/issues/new)

Fill out the relevant sections of the [`ISSUE_TEMPLATE.md`](https://github.com/PrincetonUniversity/athena/blob/master/.github/ISSUE_TEMPLATE.md) to the best of your ability when submitting a new issue.

### Suggesting enhancements
Feature requests are welcome, and are also tracked as [GitHub issues]( https://guides.github.com/features/issues/).

Please understand that we may not be able to respond to all of them because of limited resources.

## Submitting changes
Some requirements for code submissions:
- Athena++ is licensed under the BSD 3-Clause License; contributions must also use the BSD-3 license.
- The code must be commented and well documented, see [Documentation](#documentation).
- The Athena++ Wiki has a [Style Guide](https://github.com/PrincetonUniversity/athena/wiki/Style-Guide) section in the Programmer Guide. Please follow these conventions as closely as possible in order to promote consistency in the codebase.
- When implementing new functionality, add a regression test. See [Testing and continuous integration (CI)](#testing-and-continuous-integration-CI).
- If your submission fixes an issue in the [issue tracker](https://github.com/PrincetonUniversity/athena/issues), please reference the issue # in the pull request title or commit message, for example:
```
Fixes #42
```

The below instructions assume a basic understanding of the Git command line interface.
If you are new to Git or a need a refresher, the [Atlassian Bitbucket Git tutorial](https://www.atlassian.com/git/tutorials/learn-git-with-bitbucket-cloud) and the [Git documentation](https://git-scm.com/) are helpful resources.

The easiest way to contribute to Athena++ is to fork the repository to your GitHub account, create a branch on your fork, and make your changes there. When the changes are ready for submission, open a pull request (PR) on the Athena++ repository. The workflow could be summarized by the following commands:
1. Fork the repository to your GitHub account (only once) at https://github.com/PrincetonUniversity/athena/fork
2. Clone a local copy of your fork:
```
git clone https://github.com/<username>/athena ./athena-<username>
```
3. Create a descriptively-named feature branch on the fork:
```
cd athena-<username>
git checkout -b cool-new-feature
```
4. Commit often, and in logical groups of changes.
  * Use [interactive rebasing](https://help.github.com/articles/about-git-rebase/) to clean up your local commits before sharing them to GitHub.
  * Follow [commit message guidelines](https://www.git-scm.com/book/en/v2/Distributed-Git-Contributing-to-a-Project#_commit_guidelines); see also [How to Write a Git Commit Message](https://chris.beams.io/posts/git-commit/).
```
git add src/modified_file.cpp
# Use your editor to format the commit message
git commit -v
```
5. Push your changes to your remote GitHub fork:
```
git push -u origin cool-new-feature
```
6. When your branch is complete and you want to add it to Athena++, [open a new pull request to `master`](https://github.com/PrincetonUniversity/athena/pull/new/master).

### Forks and branches
The use of separate branches for both new features and bug fixes, no matter how small, is highly encouraged. Committing directly to `master` branch should be kept to a minimum. [Branches in Git are lightweight](https://git-scm.com/book/en/v2/Git-Branching-Branches-in-a-Nutshell), and merging small branches should be painless.

For the majority of development, users should use personal forks instead of branches on [PrincetonUniversity/athena](https://github.com/PrincetonUniversity/athena) (especially for larger development projects). The shared Athena++ repository should only contain a restricted set of main feature branches and temporary hotfix branches at any given time. <!-- consider reaching out to Athena++ developers before starting any significant PR/feature development to see if anyone is working on it or if we would consider merging it into Athena++-->

To update your fork with changes from [PrincetonUniversity/athena](https://github.com/PrincetonUniversity/athena), from the `master` branch on a cloned copy of the forked repo:
1. Add a remote named `upstream` for the original Athena++ repository:
```
git remote add upstream https://github.com/PrincetonUniversity/athena
```
2. Fetch the updates from the original Athena++ repository:
```
git fetch upstream
```
3. Merge the new commits into your forked `master`:
```
git merge --ff-only upstream/master
```
will work if you have not committed directly to your forked `master` branch.

If you have modified your forked `master` branch, the last two steps could be replaced by:
```
git pull --rebase upstream master
```
See [Developing on shared `branch`](#developing-on-shared-branch).

### Developing on shared `branch`
There are a few practices that should be followed when committing changes to a collaborative `branch` on [PrincetonUniversity/athena](https://github.com/PrincetonUniversity/athena) in order to avoid conflicts and headaches. These guidelines especially apply to developing on the fast changing `master` branch for those users with Admin permissions.

If you commit to an outdated local copy of `branch` (i.e. someone else has pushed changes to GitHub since you last checked), the `git push origin branch` command will be rejected by the server and prompt you to execute the `git pull` command. The default `git pull` behavior in this scenario is to create a merge-commit after you resolve any conflicts between your changes and the remote commits. However, these non-descriptive commit messages tend to clutter the repository history unnecessarily.
<!-- insert image of Network graph to compare linear and non-linear Git history -->

For example, searching the Athena++ repository history using the [GitHub website](https://github.com/PrincetonUniversity/athena/search?utf8=%E2%9C%93&q=merge+branch+%27master%27+of+https:&type=Commits) or the command line:
```
git log --oneline —grep="Merge branch ‘master' of https://github.com/PrincetonUniversity/athena$" | wc -l
```
returns many such commits. Most of them likely could have been avoided by either 1) doing local development on feature branches or 2) using `git pull --rebase` to perform a rebase instead of a merge when pulling conflicting updates.

If you frequently encounter such issues, it is recommended to enable the latter by default. In git versions >= 1.7.9, this can be accomplished with:
```
git config --global pull.rebase true
```

### Pull requests
When your changes are ready for submission, you may open a new pull request to `master` [from a branch on the main repository (Write access)](https://github.com/PrincetonUniversity/athena/pull/new/master) or from a branch on your forked repository. For the latter, go to the page for your fork on GitHub, select your development branch, and click the pull request button. Fill out the relevant sections of the [`PULL_REQUEST_TEMPLATE.md`](https://github.com/PrincetonUniversity/athena/blob/master/.github/PULL_REQUEST_TEMPLATE.md) to the best of your ability when submitting a new PR.

We will discuss the proposed changes and may request that you make modifications to your code before merging. To do so, simply commit to the feature branch and push your changes to GitHub, and your pull request will reflect these updates.

Before merging the PRs, you may be asked to squash and/or rebase some or all of your commits in order to preserve a clean, linear Git history. We will walk you through the interactive rebase procedure, i.e.
```
git rebase -i master
```

In general for Athena++, merging branches with `git merge —no-ff` is preferred in order to preserve the historical existence of the feature branch.

After the pull request is closed, you may optionally want to delete the feature branch on your local and remote fork via the GitHub PR webpage or the command line:
```
git branch d cool-new-feature
git push origin --delete cool-new-feature
```

### Code review policy
Currently, `master` is a GitHub [protected branch](https://help.github.com/articles/about-protected-branches/), which automatically:
* Disables force pushing on `master`
* Prevents `master` from being deleted

Additionally, we have enabled ["Require pull request reviews before merging"](https://help.github.com/articles/enabling-required-reviews-for-pull-requests/) to `master`. This setting ensures that all pull requests require at least 1 code review before the branch is merged to the `master` branch and effectively prohibits pushing **any** commit directly to `master`, even from users with Write access. Attempting to do so will result in an error such as:
```
Total 9 (delta 7), reused 0 (delta 0)
remote: Resolving deltas: 100% (7/7), completed with 7 local objects.
remote: error: GH006: Protected branch update failed for refs/heads/master.
remote: error: At least 1 approving review is required by reviewers with write access.
To git@github.com:PrincetonUniversity/athena.git
! [remote rejected] master -> master (protected branch hook declined)
error: failed to push some refs to 'git@github.com:PrincetonUniversity/athena.git'
```

Only collaborators with Admin permissions can bypass these restrictions. The decision to force the use of branches and pull requests for all changes, no matter how small, was made in order to:
1. Allow for isolated testing and human oversight/feedback/discussion of changes
2. Promote a [readable](https://fangpenlin.com/posts/2013/09/30/keep-a-readable-git-history/), [linear](http://www.bitsnbites.eu/a-tidy-linear-git-history/), and reversible Git history for computational reproducibility and maintainability
3. Most importantly, prevent any accidental pushes to `master`
<!-- Currently set # of required reviews to 1; other options to consider enabling in the future include: -->
<!-- "Dismiss stale PR approvals when new commits are pushed" -->
<!-- "Restrict who can push to this branch" (redundant with Require PR reviews)-->
<!-- "Require status checks to pass before merging" after separating CI build steps in new GitHub Checks API-->
<!-- "Include administrators" -->

When anyone opens a new pull request to `master`, GitHub will automatically request a code review from one or more users defined by the PR's modified files and the rules in the current [`.github/CODEOWNERS`](https://github.com/PrincetonUniversity/athena/blob/master/.github/CODEOWNERS) file. Only users with Admin permissions may modify this file to designate collaborators with at least Write access as "code owners". It is possible to use separate versions of this file on each branch to regulate PRs targeting those branches; see "[About CODEOWNERS](https://help.github.com/articles/about-codeowners/)" for more information.

## Testing and continuous integration (CI)
Automated testing is an essential part of any large software project. The [Regression Testing](https://github.com/PrincetonUniversity/athena/wiki/Regression-Testing) page in the Athena++ Wiki describes how to use and write new tests for the framework setup in the `tst/regression/` folder. Developers should run these tests to ensure that code changes did not break any existing functionalities.

Continuous integration is currently provided by both the Princeton Jenkins server and Travis CI service. These services automatically use the [Regression Testing](https://github.com/PrincetonUniversity/athena/wiki/Regression-Testing) framework to check code functionality and code <a href="https://en.wikipedia.org/wiki/Lint_(software)">linters</a> to ensure that conventions in the [Style Guide](https://github.com/PrincetonUniversity/athena/wiki/Style-Guide) are obeyed. The details of the infrastructure setup and instructions on how to use these services are covered in the [Continuous Integration (CI)](https://github.com/PrincetonUniversity/athena/wiki/Continuous-Integration-%28CI%29) Wiki page.

## Documentation
The development repository's [documentation](https://github.com/PrincetonUniversity/athena/wiki) is a [GitHub Wiki](https://help.github.com/articles/about-github-wikis/) and is written largely in Markdown. Limited math typesetting is supported via HTML. See existing Wiki source for examples, e.g. [Editing: Coordinate Systems and Meshes](https://github.com/PrincetonUniversity/athena/wiki/Coordinate-Systems-and-Meshes/_edit).

Any significant change or new feature requires accompanying documentation before being merged to `master`. While edits can be made directly using the online interface, the Wiki is a normal Git repository which can be cloned and modified. [However](https://help.github.com/articles/adding-and-editing-wiki-pages-locally/):
> You and your collaborators can create branches when working on wikis, but only changes pushed to the `master` branch will be made live and available to your readers.

## Community
The Athena++ private Slack workspace is located at [athena-pp.slack.com](https://athena-pp.slack.com). The default `#general` and `#random` channels are available for free-form discussion and user support, and topic-specific channels and private Direct Messages (DMs) with up to 8 other members can be started by anyone. Issues and pull requests on the GitHub repository should still be the main forum to discuss development details, but the Slack workspace is a useful centralized forum for general discussion, sharing new results, asking questions, and learning what others are working on. This Slack workspace was setup on the Free plan, which essentially limits the amount of file storage to 5GB and message history to 10k messages

At this time, the Slack workspace is closed to the general public. The workspace is configured such that anyone with a `@princeton.edu` email can join automatically at [this signup link](https://join.slack.com/t/athena-pp/signup). Any current member may invite new members.

## Versioning and releases
We intend to provide periodic releases, versioned according to CalVer, or [Calendar Versioninng](https://calver.org/). A detailed walkthrough of the steps a project maintainer must complete in order to mint a new release is provided in the following section.

We currently maintain Git tags and code versions in a one-to-one correspondence: all versions are tagged, and all tags have a version number. However, not all tags/versions are released (see below section discussing pre-release tagged versions). A release version is defined by drafting a [GitHub Release](https://help.github.com/articles/creating-releases/) along with release notes in the GitHub UI. 

Each release is accompanied by an Git annotated (not lightweight) tag. An annotated tag is a full Git object with its own tagger name, tagger email, and creation date. A lightweight Git tag is more appropriate for temporary or local/personal use than for publishing releases, since a lightweight tag is merely a pointer to a commit object (much like a branch that doesn't naturally move with commits and *shouldn't* be moved by users after it is shared). Therefore, the tag should be created from the Git CLI, not the GitHub UI which only supports creating lightweight tags as of 5/24/18.

### Example steps for tagging and releasing a new version
*Last updated April 2021*.

#### Prepare to increment version
Start in the root project directory of a fully up-to-date clone of the repository, with all changes committed and pushed. For convenience, we assume this is located at `~/athena/`.
1. Edit `src/main.cpp` to update the following line with the new version number and current month and year for proper `athena -h` output:
```c++
std::string athena_version = "version 21.0 - January 2021";
```
2. Use [`github_changelog_generator`](https://github.com/github-changelog-generator/github-changelog-generator) tool to help summarize changes since the previous version
  - The automatic process works best if the titles of the recently closed Issues and merged Pull Requests are first reviewed and retroactively edited to use proper active voice, formatting, and descriptiveness. E.g. "Add continuous integration with Travis CI and Jenkins" is a good Pull Request title.
  - This Ruby-based software can be installed locally via `gem install github_changelog_generator`. These instructions were written while using version `1.15.0-beta` of `github_changelog_generator`.
  - The GitHub user who is executing the tool should first [Create a personal access token for the command line](https://help.github.com/articles/creating-a-personal-access-token-for-the-command-line/) to avoid the API unauthenticated request limit.
  - Using this secret token value, execute `github_changelog_generator -t XXXX -o temp.md --since-tag v20.0`
  - Copy/paste the relevant sections from `temp.md` into `CONTRIBUTING.md` by replacing the old `Unreleased` section and creating sections for the new tag. The new sections can be created automatically by the script if it is instead executed after the tags exist.
  - Add `### Issues and Pull Requests:` subsection header above the new subsections. Clean up any redundant/irrelevant/confusing entries.
  - Write high-level summaries in new subsection(s) `### Added`, `### Removed`, `### Fixed/Changed` features using the automatically-generated entries as reference; these sections will be copied to the release notes.
  - Currently, the `CHANGLEOG.md` is completely excluded from the public repository; it is only for reference by Athena++ developers.
3. Add, commit, and push all of these bookkeeping changes (skipping any continuous integration) with `git add CHANGELOG.md src/main.cpp; git commit; git push` using a formulaic commit message such as:
> Bump version number and update CHANGLEOG
>
> [ci skip]

#### Create repository tag and GitHub Release

4. Push the 2x new tags to the GitHub remote with `git push --tags`
5. [Draft a new Release on the repository](https://github.com/PrincetonUniversity/athena/releases/new) using the `v21.0` tag.
  - Use a simple title, *Athena++ 21.0*, and header `## Release 21.0`
  - Copy/paste the newly-updated "Added", "Removed", "Fixed/Changed" sections from `CHANGLEOG.md` into the notes and write an introduction.

6. Announce the latest public release on the [Athena++ website](https://www.athena-astro.app/):
```
cd ~/athena; git checkout gh-pages
emacs download.html # Edit the latest version number and write description
git commit -am "Announce v21.0 public release"
git push
```

<!-- # Athena++ Code of Conduct
(add section here or store in external file, CODE_OF_CONDUCT.md) -->
