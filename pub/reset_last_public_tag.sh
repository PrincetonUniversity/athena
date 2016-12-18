#! /bin/bash
# Resets the last_public Git tag to a given commit specified by SHA1.
# usage: reset_last_public_tag.sh sha1
# example: reset_last_public_tag.sh 586eb4a

# remove local tag if it exists
if ! [ `git tag | grep -q last_public` ]
then
    git tag -d last_public
fi

# remove remote tag if it exists
if ! [ `git ls-remote --tags origin | grep -q last_public` ]
then
    git push --no-verify origin :refs/tags/last_public
fi

# create new local tag
git tag last_public $1

# push to remote
git push origin last_public
