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
 * **Note**: the user's private fork will be deleted on GitHub if the individual is ever removed as a [collaborator or from a Team and loses Read access](https://help.github.com/articles/adding-outside-collaborators-to-repositories-in-your-organization/).

[GitHub Teams](https://help.github.com/articles/organizing-members-into-teams/), `Athena++_PP` (Write) and `Athena++` (Admin), were originally used to control permission levels for groups of users; however, this is being phased out in favor of a newer [Nested Teams](https://blog.github.com/2017-06-13-nested-teams-add-depth-to-your-team-structure/) setup with:
* `Athena++` (Read)
   * `Athena++_core` (Write)
     * `Athena++_admin` (Admin)

These teams are currently setup with the same permissions hierarchy for the public version repository. Several Admin users also have [Team Maintainer](https://help.github.com/articles/giving-team-maintainer-permissions-to-an-organization-member/) roles for managing the membership and permissions of these teams. A GitHub user must be a member of the PrincetonUniversity Organization on GitHub. Princeton Research Computing has a [GitHub request form](https://forms.rc.princeton.edu/github/) that can be used to add usernames to the PrincetonUniversity GitHub organization (even if the user does not have a Princeton Netid). This is the preferred method for managing permissions; repository Collaborator status should be used only for short-term access.

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
We intend to provide periodic releases to the
[athena-public-version](https://github.com/PrincetonUniversity/athena-public-version)
repository based on the versions created in the private repository.

In the private repository, we currently maintain Git tags and code versions in a one-to-one correspondence: all versions are tagged, and all tags have a version number. However, not all tags/versions are released (see below section discussing pre-release tagged versions). A release version is defined by drafting a [GitHub Release](https://help.github.com/articles/creating-releases/) along with release notes in the GitHub UI. Once a release, say [Athena++ 1.1.1-dev](https://github.com/PrincetonUniversity/athena/releases/tag/v1.1.1-dev), is created in the private repository, a closely-related release with same version number, [Athena++ 1.1.1](https://github.com/PrincetonUniversity/athena-public-version/releases/tag/v1.1.1), should soon be drafted in the public repository. See the header of [`CHANGELOG.md`](./CHANGELOG.md) for notes on Semantic Versioning and GitHub Release note practices for Athena++.
<!-- the above bijection is true for now, but the latter (surjective) restriction may be loosened in the future, e.g. tags such as "last_public" -->

In both the private and public repositories, each release is accompanied by an Git annotated (not lightweight) tag. An annotated tag is a full Git object with its own tagger name, tagger email, and creation date. A lightweight Git tag is more appropriate for temporary or local/personal use than for publishing releases, since a lightweight tag is merely a pointer to a commit object (much like a branch that doesn't naturally move with commits and *shouldn't* be moved by users after it is shared).

Therefore, the tag should be created from the Git CLI, not the GitHub UI which only supports creating lightweight tags as of 5/24/18. The name of the tag corresponds to the same release number, prefixed with a `v` to distinguish them from other Git objects--- e.g. `vX.Y.Z`. See [Is "v1.2.3" a semantic version?](https://github.com/semver/semver/blob/master/semver.md#is-v123-a-semantic-version)

Before creating new versions and Git tags, the following manual actions should be taken:
* The version string output by the `athena -h` command must be manually
updated in `src/main.cpp`.  
* The user should also manually revert the `README.md` file to the simple (no CI status badges nor link to this file) format before tagging on the public repository.

After pushing new version changes to
[athena-public-version](https://github.com/PrincetonUniversity/athena-public-version),
the following updates should be performed manually:
- Duplicate the private version tag in the public repository
  - There is a manually-managed, ad-hoc relationship between the private repository's and the public repository's tags and releases; in the future, this could be improved by replacing the `git reset --soft <base>` squash procedure performed by `./build_pub_repo_v2.sh`.
- Draft public [Release
  Notes](https://github.com/PrincetonUniversity/athena-public-version/releases)
  based on private [`CHANGELOG.md`](./CHANGELOG.md) entries
- Update the [public Athena++
  Wiki](https://github.com/PrincetonUniversity/athena-public-version/wiki)
  based on the [private Athena++ Wiki](https://github.com/PrincetonUniversity/athena/wiki).
  - Since version 1.1.0, there is a `public` branch on the `athena.wiki` repository that should be continually rebased on top of `master` with a few commits that remove the irrelevant content. The `master` branch of `athena-public-version.wiki` should mirror that private wiki `public` branch via force-pushing.
  - Previously, this consisted of manually force-pushing from the private repository's `master` branch to the public repository's `master` branch *and then* removing or modifying any content that is sensitive/secret or irrelevant to the public repository (e.g. the pages on "Continuous Integration").
- Announce release on the [Athena++
  website](https://princetonuniversity.github.io/athena/index.html) by
  modifying the HTML files on the `gh-pages` branch of the private repository

A detailed example of these steps are illustrated in a below section.

### Pre-release versions, release/support branches, tagging practices
While the release version is ideally drafted from a tagged commit on the `master` branch, it is possible that this development branch contains source code that is not intended to be released to the public. Or, the release may require temporary modifications that should not pollute the history of `master`. In such cases, the content should be removed and the changes should be made in a separate **release/support branch**. A good formatting convention for naming release branches is `release/X.Y.Z`. While these dead-end branches will not be merged back to master, they should be deleted after being tagged/released. The annotated tag at the tip of the branch should ensure that these commits are never garbage collected by Git; all of the commits in the former `release/X.Y.Z` branch remain reachable by walking back from the `vX.Y.Z` tag.

For prerelease tags, `X.Y.Z-suffix` is allowed, with optional suffixes `-alpha`, `-beta`, `-rc`, `-dev`. Their intended meanings are:
- `-alpha`: version is feature-incomplete and/or unstable
- `-beta`: version is feature-complete, but may have bugs and/or may be significantly changed before release
- `-rc`: (release candidate) version will be released unless a last-minute bug is identified

The above version suffixes may all have an integer appended to them, e.g. `-rc.1`, `-rc.2`, ... . In contrast, the following tag suffix should not have additional numbers appended at the end:
- `-dev`: pre-release version that contains source code not meant for public release

### Example steps for tagging and releasing new public version
*Last updated 8/17/18*.

We now illustrate the many policies and practices discussed in this section using a concrete example. Suppose that the latest released version of Athena++ is 1.1.1 on both the private [`athena`](https://github.com/PrincetonUniversity/athena) and public [`athena-public-version`](https://github.com/PrincetonUniversity/athena-public-version) repositories (which is true at the time of writing), and suppose that the `master` branch on the private repository has changed significantly since the `v1.1.1` and/or `v1.1.1-dev` tags. If the changes merit a public release of a new MINOR revision 1.2.0, here are the actions that should be taken.

#### Prepare to increment version
Start in the root project directory of a fully up-to-date clone of the private repository, with all changes committed and pushed. For convenience, we assume this is located at `~/athena/`.
1. Edit `src/main.cpp` to update the following line with the new version number and current month and year for proper `athena -h` output:
```c++
std::string athena_version = "version 1.2.0 - August 2018";
```
2. Use [`github_changelog_generator`](https://github.com/github-changelog-generator/github-changelog-generator) tool to help summarize changes since 1.1.1
  - The automatic process works best if the titles of the recently closed Issues and merged Pull Requests are first reviewed and retroactively edited to use proper active voice, formatting, and descriptiveness. E.g. "Add continuous integration with Travis CI and Jenkins" is a good Pull Request title.
  - This Ruby-based software can be installed locally via `gem install github_changelog_generator`. These instructions were written while using version 1.15.0-beta.
  - The GitHub user who is executing the tool should first [Create a personal access token for the command line](https://help.github.com/articles/creating-a-personal-access-token-for-the-command-line/) to avoid the API unauthenticated request limit.
  - Using this secret token value, execute `github_changelog_generator -t XXXX -o temp.md --since-tag v1.1.1`
  - Copy/paste the relevant sections from `temp.md` into `CONTRIBUTING.md` by replacing the old `Unreleased` section and creating sections for the new `v1.2.0` (and possibly `v1.2.0-dev` tags (which do not yet exist). The new sections can be created automatically by the script if it is instead executed after the tags exist.
  - Add `### Issues and Pull Requests:` subsection header above the new subsections. Clean up any redundant/irrelevant/confusing entries.
  - Write high-level summaries in new subsection(s) `### Added`, `### Removed`, `### Fixed/Changed` features using the automatically-generated entries as reference; these sections will be copied to the release notes.
  - Currently, the `CHANGLEOG.md` is completely excluded from the public repository; it is only for reference by Athena++ developers.
3. Add, commit, and push all of these bookkeeping changes (skipping any continuous integration) with `git add CHANGELOG.md src/main.cpp; git commit; git push` using a formulaic commit message such as:
> Bump version number to v1.2.0 and update CHANGLEOG
>
> [ci skip]

#### Create private repository tag(s) and GitHub Release
If the `master` branch contains some feature that must be excluded from the public release (see earlier subsection), then here are the specific steps to follow:
1. Create separate annotated tag with `-dev` suffix to checkpoint the state of the private repository immediately before the upcoming public release with `git tag -s v1.2.0-dev`
  - Use `-s` to create an annotated tag object and sign it with your local GPG public key.
  - Use `-a` if no local GPG key exists or an unsigned, annotated tag object is desired.
  - Since an annotated tag is being created and the `-m` option is absent, Git will open an editor to create the tag message. Write something descriptive, such as:
  > Development repo master branch before removing Multigrid (again)

2. On tip of `master` branch, create a "dead-end" branch for polishing the public release: `git checkout -b release/1.2.0`
3. Manually modify/remove any source code that should remain private and cannot be removed by simply excluding entire files. Commit these changes to the dead-end branch `release/1.2.0` and `git push` to the GitHub remote.
4. Create annotated tag for the public release version by running `git tag -s v1.2.0` on the tip of this branch.
  - Follow the same advice as the above tag. For example, a good tag message might be:
  > Public release of Athena++ version 1.2.0
  >
  > Multigrid gravity removed from v1.2.0-dev tag. Created on dead-end branch release/1.2.0, which will be deleted.

5. Push the 2x new tags to the private repository's GitHub remote with `git push --tags`
6. [Draft a new Release on the private repository](https://github.com/PrincetonUniversity/athena/releases/new) using the `v1.2.0-dev` tag.
  - Even though the `v1.2.0` tag will have the latest creation date, the `v1.2.0-dev` should be considered the most recent version on the private repository in this scenario.
  - Use a simple title, *Athena++ 1.2.0-dev*, and header `## Release 1.2.0`
  - Copy/paste the newly-updated "Added", "Removed", "Fixed/Changed" sections from `CHANGLEOG.md` into the notes and write an introduction.

If nothing needs to be manually removed from `master` before the public release, then only the final 3 steps need to be executed. There will be no `v1.2.0-dev` tag, and the `v1.2.0` tag will be the latest release on the private repository.

#### Modify and use private -> public release script
1. Create a copy of the public release helper-script outside local working copy: `cp ./pub/build_public_repo-v2.sh ~/; cd ~`
2. Edit the copy of `build_public_repo-v2.sh` to change the variable indicating which branch of the private repository is the source for this public release:  `PRIV_REPO_PRIV_BRANCH="release/1.2.0"` (cannot be tag object `v1.2.0`!)
  - Also, possibly edit the Bash array variable `PRIVATE` to add/remove any files or subfolders that should be completely excluded from any public version. Examples include the `pub/` subdirectory containing the version-controlled script and instructions for creating such releases.
3. Execute the script by creating a new local clone `./build_pub_repo_v2.sh athena_working; cd athena_working`
  - The script may take several minutes to run, since it scans thousands of Git commits to history of files listed in the `PRIVATE` array.
4. Check that no private source code remains in the `athena_working/` files or filtered Git history.
5. Edit `README.md` to remove the private repository's CI status badges and link to `CONTRIBUTING.md` and stage this change.
5. Commit the squashed changes with `git commit -m "Public release of Athena++ version 1.2.0"`

#### Create public repository tag and GitHub Release
1. Create new annotated tag public repository counterpart of private repository annotated tag `v1.2.0` with `git tag -s v1.2.0`
  - Although the state of the public repository at this `v1.2.0` tag should be similar to the state of the private repository at its `v1.2.0` tag, the two tags objects have no direct Git relationship (different object hashes, tag times, etc.).
  - Typically, this tag message will be identical or similar to the previous commit message, e.g.:
  > Public release of Athena++ version 1.2.0

2. Share the new commit and tag with the public repository's GitHub remote with `git push --tags`
3. [Draft a new Release on the public repository](https://github.com/PrincetonUniversity/athena-public-version/releases/new) based on this `v1.2.0` tag. The notes can be copied from the earlier private release notes; remove any content that is irrelevant to the public version.

#### Publish Wiki and announce release
Now, we assume we have up-to-date, fresh clones of the private and public GitHub Wikis located at `~/athena.wiki/` and `~/athena-public-version.wiki/`, respectively.
1. For the private Wiki, rebase the `public` branch on the latest `master` branch with `cd ~/athena.wiki; git checkout public; git rebase master` and resolve any conflicts.
2. Make, stage, and commit any additional modifications to the `public` branch as desired and then force push the rebased branch with `git push -f`.
3. Force the public Wiki's `master` branch to mirror the private Wiki's `public` branch:
```
cd ~/athena-public-version.wiki
git remote add upstream git@github.com:PrincetonUniversity/athena.wiki.git`
git fetch upstream
git branch master --set-upstream-to upstream/public
git pull --rebase
git push -f origin master
```
4. Announce the latest public release on the [Athena++
  website](https://princetonuniversity.github.io/athena/index.html):
```
cd ~/athena; git checkout gh-pages
emacs download.html # Edit the latest version number and write description
git commit -am "Announce v1.2.0 public release"
git push
```

#### Cleanup
1. Delete the private repository's release branch, either via the online GitHub UI or `cd ~/athena/; git branch -D release/1.2.0; git push origin --delete release/1.2.0`
2. Delete the temporary build directory of the public release: `rm -rfd ~/athena-working`

<!-- Note, some projects don't use a period after numbered suffixes, e.g. -rc2-->
<!-- Although SemVer isn't so prescriptivist to define a set of suffixes, https://github.com/semver/semver/issues/114, it does define sytnax and how to use the suffixes in pt 9 , so that "pre-release versions have a lower precedence than the associated normal version" -->

<!-- TODO: add notes on automated CHANGELOG.md generation, drafting Release notes (hopefully automated via GitHub API in future, such as conventional-changelog tool), importance of labeling PRs/Issues + using active voice grammar in their titles, the need for documenting a public API -->

<!-- Need to add notes on modification from old procedure with pub/build_pub_repo.sh to new pub/build_pub_repo_v2.sh script; do we want to restart use of the last_public lightweight tag? -->

<!-- # Athena++ Code of Conduct
(add section here or store in external file, CODE_OF_CONDUCT.md) -->
