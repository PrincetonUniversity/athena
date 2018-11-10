#!/usr/bin/env bash
# usage: build_pub_repo.sh local_repo_name
# example: build_pub_repo.sh athena_working

# fileglob of private files and directories
# the filename needs to include the path with respect
# to the top of the Git repo
# Here is an example for excluding multiple directories
# private_files="{pub,private,inputs}"
# Note: for a single directory, you should not use the curly
# brackets e.g.
# incorrect: private_files="{pub}"
#   correct: private_files="pub"
private_files="{pub,tst/ci,.travis.yml,.github,CONTRIBUTING.md,CHANGELOG.md}"

# Repo HTTPS URLs
priv_repo_url="https://github.com/PrincetonUniversity/athena.git"
pub_repo_url="https://github.com/PrincetonUniversity/athena-public-version.git"
# Repo SSH URLs:
# priv_repo_url="git@github.com:PrincetonUniversity/athena.git"
# pub_repo_url="git@github.com:PrincetonUniversity/athena-public-version.git"

# data for this script: assumes we are rebasing private repo "master" on "last_public" tag
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

# Move the tag last_public
## copy current tag
if ! [ `git tag | grep -q last_public` ]
then
    git tag last_public_old last_public
fi

## remove local tag if it exists
if ! [ `git tag | grep -q last_public` ]
then
    git tag -d last_public
fi

## remove remote tag if it exists
if ! [ `git ls-remote --tags origin | grep -q last_public` ]
then
    git push --no-verify origin :refs/tags/last_public
fi

## Tag the latest master commit
git tag last_public `git rev-parse master`

## push the new tag to origin
git push --no-verify origin last_public

# rename last_public_old to last_public
git tag -d last_public
git tag last_public last_public_old
git tag -d last_public_old

# Copy the pre-push hook to the cloned repo
if [ -e $pre_push_hook ]
then
    cp $pre_push_hook .git/hooks/
    chmod 766 .git/hooks/pre-push
else
    echo the pre-push hook does not exists
    exit 1
fi

# clean up the files and directories we don't want to make public
git filter-branch --force --index-filter \
    "git rm -r --cached --ignore-unmatch $private_files" \
    --prune-empty --tag-name-filter cat -- --all

# rebase to eliminate merge commit since last_public
git rebase last_public

# add the public_repo remote
git remote add $pub_remote_name $pub_repo_url
git fetch --all

# Create the public branch by pulling the master branch of
# the public repo. The public branch of the private repo is
# now tracking the master branch of the public repo so that
# "git push" and "git pull" will be towards public_repo:master
git checkout -b $priv_repo_pub_branch $pub_remote_name/$pub_repo_branch

# configure git so it pushes from the public branch to public_repo:master
git config push.default upstream

# purge the deleted files from the packfile
git for-each-ref --format='delete %(refname)' refs/original | git update-ref --stdin
git reflog expire --expire=now --all
git gc --aggressive --prune=now
