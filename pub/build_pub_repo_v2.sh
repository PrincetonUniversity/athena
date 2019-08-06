#!/usr/bin/env bash
# usage: build_pub_repo_v2.sh local_repo_name
# example: build_pub_repo_v2.sh athena_working

# fileglob of private files and directories
# the filename needs to include the path with respect
# to the top of the Git repo
# Here is an example for excluding multiple directories
# private_files="{pub,private,inputs}"
# Note: for a single directory, you should not use the curly
# brackets e.g.
# incorrect: private_files="{pub}"
#   correct: private_files="pub"
private_files="{pub,tst/ci,.travis.yml,.codecov.yml,.github,.github_changelog_generator,CONTRIBUTING.md,CHANGELOG.md}"

# Repo HTTPS URLs:
# priv_repo_url="https://github.com/PrincetonUniversity/athena.git"
# pub_repo_url="https://github.com/PrincetonUniversity/athena-public-version.git"
# Repo SSH URLs:
priv_repo_url="git@github.com:PrincetonUniversity/athena.git"
pub_repo_url="git@github.com:PrincetonUniversity/athena-public-version.git"

# data for this script
## branch on remote private repo from which to create public release
## (must be an actual branch and not a tag in order to rewrite local history/filter-branch)
priv_repo_priv_branch="release/1.1.1"
## local workspace branch for drafting the public release
priv_repo_pub_branch="public"

pub_remote_name="public_repo"
pub_repo_branch="master"

pre_push_hook="pub/pre-push"
post_commit_hook="pub/post-commit"

local_repo_name=$1

# clone private repo
git clone $priv_repo_url $local_repo_name
cd $local_repo_name
git pull --tags

# Copy the pre-push hook to the cloned repo
if [ -e $pre_push_hook ]
then
    cp $pre_push_hook .git/hooks/
    chmod 766 .git/hooks/pre-push
else
    echo the pre-push hook does not exists
    exit 1
fi

# get the private repo branch to be made public
# (if this is called after the "git filter-branch ... " command, the checkout
#  will return a non-filtered copy of the remote branch)
git checkout $priv_repo_priv_branch

# clean up the files and directories we don't want to make public
git filter-branch --force --index-filter \
    "git rm -r --cached --ignore-unmatch $private_files" \
    --prune-empty --tag-name-filter cat -- --all

# add the public_repo remote
git remote add $pub_remote_name $pub_repo_url
git fetch --all # will not overwrite the filtered tags

# Create the "public" branch by pulling the master branch of
# the public repo. The "public" branch of the private repo is
# now tracking the master branch of the public repo so that
# "git push" and "git pull" will be towards public_repo:master
git checkout -b $priv_repo_pub_branch
git branch -u $pub_remote_name/$pub_repo_branch $priv_repo_pub_branch
# squash commits on previous public version
git reset --soft $pub_remote_name/$pub_repo_branch

# configure git so it pushes from the public branch to public_repo:master
git config push.default upstream

# purge the deleted files from the packfile
git for-each-ref --format='delete %(refname)' refs/original | git update-ref --stdin
git reflog expire --expire=now --all
git gc --aggressive --prune=now

# git diff public..public_repo/master
# git commit -m "Second public release of Athena++"
# git branch -vv
