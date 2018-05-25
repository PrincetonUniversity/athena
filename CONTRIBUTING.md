# Introduction
Welcome! Thank you for considering contributing to Athena++.

The guidelines in this document are meant to help make the development of Athena++ straightforward and effective. They are a set of best practices, not strict rules, and this document may be modified at any time. Navigating the code can be daunting for new users, so if anything is unclear, please let us know!

<!-- ### Table of Contents -->

## Resources and quick links
* The latest development version of Athena++ is hosted in the private [PrincetonUniversity/athena](https://github.com/PrincetonUniversity/athena) GitHub repository.
* The public version of Athena++ is distributed via the public [PrincetonUnviersity/athena-public-version](https://github.com/PrincetonUniversity/athena-public-version) GitHub repository.
  * Periodic releases from the development repository are transferred by following the procedure described in [`pub/doc/pub_doc.tex`](https://github.com/PrincetonUniversity/athena/blob/master/pub/doc/pub_doc.tex) using the scripts in `pub/`.
<!-- Establish timeline for periodic releases? -->
* The final version of the predecessor C code, Athena 4.2, has been migrated from its original SVN repository to GitHub at [PrincetonUniversity/Athena-Cversion](https://github.com/PrincetonUniversity/Athena-Cversion).
  * [Athena 4.2 Documentation](https://princetonuniversity.github.io/Athena-Cversion/AthenaDocs) is hosted on GitHub Pages.
  * The [Athena 4.2 Test Page](https://www.astro.princeton.edu/~jstone/Athena/tests/) contains useful algorithm test results.
* The [Athena++ website](https://princetonuniversity.github.io/athena/index.html) is hosted by [GitHub Pages](https://pages.github.com/).
* The latest version of the [Athena++ documentation](https://github.com/PrincetonUniversity/athena/wiki) is hosted as a GitHub Wiki attached to the private repository.
* The Athena++ Slack workspace is located at [athena-pp.slack.com](https:://athena-pp.slack.com).
<!-- Could add links to New PR, New Issue, Issue Labels e.g. current "bugs" -->

## Organization

The [athena](https://github.com/PrincetonUniversity/athena) development repository is a private GitHub repo owned by the [PrincetonUniversity](https://github.com/PrincetonUniversity) organization (which is owned by the user [@cses](http://www.princeton.edu/researchcomputing/about/picscie/)). See "[About Organizations](https://help.github.com/articles/about-organizations/)" for more information.

There are three possible levels of [permissions for a repository belonging to an organization](https://help.github.com/articles/repository-permission-levels-for-an-organization/), and all are used in Athena++ development:
* **Admin**: the original/core Athena++ developers and a few users with extensive Athena++ development experience. <!-- They are responsible for maintaining the stability and quality of the overall codebase -->
* **Write**: a handful of experienced Athena++ developers who frequently contribute to significant portions of the codebase and/or work directly on a branch.
  * These users may push directly to the private repository's branches, with some caveats. See below section on the [Code review policy](#code-review-policy).
* **Read**: a larger group of users, developers, and students who have a compelling need for the latest development version of Athena++.
  * Read access allows the user to create a private forked repository, `<username>/athena`. The collaborator has complete freedom to experiment with changes within the fork without affecting the original project.
 * **Note**: the user's private fork will be deleted on GitHub if the individual is ever removed as a [collaborator](https://help.github.com/articles/adding-outside-collaborators-to-repositories-in-your-organization/).

[GitHub Teams](https://help.github.com/articles/organizing-members-into-teams/), `Athena++_PP` (Write) and `Athena++` (Admin), were originally used to control permission levels for groups of users; however, this is being phased out in favor of setting individual collaborator permissions.

The public version of Athena++ should serve as the primary resource for the majority of users. If you are reading this document, it means you have at least Read permissions to the private repository!

### Notification settings
By default, users with Read access to the private repository will only be subscribed to *participating notifications*, e.g. if the user makes or comments on an issue/PR/commit or if someone explicitly mentions their username. *Watching notifications* are a stronger notification setting that, typically, users have to explicitly opt-in for. See "[About notifications](https://help.github.com/articles/about-notifications/)" for more information. However, users with Write access are [automatically signed up for *watching notifications*](https://help.github.com/articles/watching-and-unwatching-repositories/) if their account has "Automatically watch repositories" enabled (default setting).

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
2. Clone a local copy of your private fork:
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

To update your private fork with changes from [PrincetonUniversity/athena](https://github.com/PrincetonUniversity/athena), from the `master` branch on a cloned copy of the forked repo:
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
When your changes are ready for submission, you may open a new pull request to `master` [from a branch on the main repository (Write access)](https://github.com/PrincetonUniversity/athena/pull/new/master) or from a branch on your private forked repository. For the latter, go to the page for your fork on GitHub, select your development branch, and click the pull request button. Fill out the relevant sections of the [`PULL_REQUEST_TEMPLATE.md`](https://github.com/PrincetonUniversity/athena/blob/master/.github/PULL_REQUEST_TEMPLATE.md) to the best of your ability when submitting a new PR.

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

### Slack Apps
The `#development` channel receives messages from the development GitHub repository when commits are made to `master` or an Issue/PR is opened or closed. This channel also receives messages from the Jenkins and TravisCI Slack Apps, which provide summaries of and links to every [Continuous Integration (CI)](https://github.com/PrincetonUniversity/athena/wiki/Continuous-Integration-%28CI%29) build test result.

Note, the GitHub + Slack integration was originally managed via a GitHub Webhook and "GitHub Notifications (Legacy)" Slack App until 5/19/2018. The current [GitHub + Slack App](https://slack.github.com/) interface enables embedded [rich link previews ](https://github.com/integrations/slack) when posting links from GitHub:
> When a user posts a GitHub link to issues and pull requests, directly linked comments, code blobs with line numbers, as well as organizations, repositories, and users in Slack, a preview of the link will be shown.

You may need to use the `/invite @github` command in private channels and Direct Messages to get rich link previews to work there. Posting links to the private [PrincetonUniversity/athena](https://github.com/PrincetonUniversity/athena) repository requires that you to link your Slack and GitHub accounts with the `/github signin` command. The GitHub + Slack App also allows you to open and close Issues and PRs from Slack with `/github close [issue link]`, for example.

Slack's simple file upload and sharing features are especially useful when compared to GitHub or email. Slack also integrates with cloud file storage apps such as Dropbox, Google Drive, and Box.

At this time, the Slack workspace is closed to the general public, but it is open to anyone who has Read access to the private repository and their associates. The workspace is configured such that anyone with a `@princeton.edu` email can join automatically at [this signup link](https://join.slack.com/t/athena-pp/signup). Any current member may invite new members. If all else fails, send your email address to [kfelker@math.princeton.edu](mailto:kfelker@math.princeton.edu) to request an invite.

## Versioning and public releases   
We intend to provide periodic releases to the [athena-public-version](https://github.com/PrincetonUniversity/athena-public-version)  repository based on the versions created in the private repository. See the header of [`CHANGELOG.md`](./CHANGELOG.md) for notes on Semantic Versioning and release practices for Athena++.

In the priave repository, each release is accompanied by an Git annotated (not lightweight) tag. Therefore, the tag should be created from the Git CLI, not the GitHub UI which only supports creating lightweight tags as of 5/24/18. The name of the tag corresponds to the same release number, prefixed with a `v` to distinguish them from other Git objects--- e.g. `vX.Y.Z`. See [Is "v1.2.3" a semantic version?](https://github.com/semver/semver/blob/master/semver.md#is-v123-a-semantic-version)

### Pre-release versions, release/support branches, tagging practices
While the release version is ideally drafted from a tagged commit on the `master` branch, it is possible that this development branch contains source code that is not intended to be released to the public. Or, the release may require temporary modifications that should not pollute the history of `master`.  In such cases, the content should be removed and the changes should be made in a separate **release/support branch**. Release branches should be named `release/X.Y.Z`. While these dead-end branches may not be merged back to master, they should be deleted after being tagged/released. The annotated tag at the tip of the branch should ensure that these commits are never garbage collected by Git.

For prerelease tags, `X.Y.Z-suffix` is allowed, with optional suffixes `-alpha`, `-beta`, `-rc`, `-dev`. Their intended meanings are:
- `-alpha`: version is feature-incomplete and/or unstable
- `-beta`: version is feature-complete, but may have bugs and/or may be significantly changed before release
- `-rc`: (release candidate) version will be released unless a last-minute bug is identified

The above version suffixes may all have an integer appended to them, e.g. `-rc.1`, `-rc.2`, ... . In contrast, the following tag suffix should not have additional numbers appended at the end:
- `-dev`: pre-release version that contains source code not meant for public release

<!-- Note, some projects don't use a period after numbered suffixes, e.g. -rc2-->
<!-- Although SemVer isn't so prescriptivist to define a set of suffixes, https://github.com/semver/semver/issues/114, it does define sytnax and how to use the suffixes in pt 9 , so that "pre-release versions have a lower precedence than the associated normal version" -->

<!-- TODO: add notes on automated CHANGELOG.md generation, drafting Release notes (hopefully automated via GitHub API in future), updating public Wiki, renaming version # in configure.py, updating gh-pages, importance of labeling PRs/Issues + using active voice grammar in their titles, the need for documenting a public API, how to create an annotated/signed tag -->

<!-- Need to add notes on modifications to pub/ scripts for releasing -->

<!-- # Athena++ Code of Conduct
(to add here or store in external file, CODE_OF_CONDUCT.md) -->
