#! /bin/bash
# usage: build_pub_repo.sh local_repo_name
# example: build_pub_repo.sh athena_working

# fileglob of private files and directories
# the filename needs to include the path with respect
# to the top of the Git repo
# Here is an example for excluding multiple directories
# PRIVATE="{pub,private,inputs}"
# Note: for a single directory, you should not use the curly
# brackets e.g.
# incorrect: PRIVATE="{pub}"
#   correct: PRIVATE="pub"
PRIVATE="{pub,tst/ci,.travis.yml}"

# Repo URLs
PRIV_REPO_URL="https://github.com/PrincetonUniversity/athena.git"
PUB_REPO_URL="https://github.com/PrincetonUniversity/athena-public-version.git"

# data for this script
PRIV_REPO_PUB_BRANCH="public"

PUB_REMOTE_NAME="public_repo"
PUB_REPO_BRANCH="master"

PRE_PUSH_HOOK="pub/pre-push"
POST_COMMIT_HOOK="pub/post-commit"

LOCAL_REPO_NAME=$1


# clone private repo
git clone $PRIV_REPO_URL $LOCAL_REPO_NAME
cd $LOCAL_REPO_NAME
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
if [ -e $PRE_PUSH_HOOK ]
then
    cp $PRE_PUSH_HOOK .git/hooks/
    chmod 766 .git/hooks/pre-push
else
    echo the pre-push hook does not exists
    exit 1
fi

# clean up the files and directories we don't want to make public
git filter-branch --force --index-filter \
    "git rm -r --cached --ignore-unmatch $PRIVATE" \
    --prune-empty --tag-name-filter cat -- --all

# rebase to eliminate merge commit since last_public
git rebase last_public

# add the public_repo remote
git remote add $PUB_REMOTE_NAME $PUB_REPO_URL
git fetch --all

# Create the public branch by pulling the master branch of
# the public repo. The public branch of the private repo is
# now tracking the master branch of the public repo so that
# "git push" and "git pull" will be towards public_repo:master
git checkout -b $PRIV_REPO_PUB_BRANCH $PUB_REMOTE_NAME/$PUB_REPO_BRANCH

# configure git so it pushes from the public branch to public_repo:master
git config push.default upstream
