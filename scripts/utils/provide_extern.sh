#!/bin/bash
###############################################################################

# Ensure we have an 'extern' folder at the base and symlink to provided args

cd ${DIR_SCRIPTS}

mkdir -p ${DIR_SCRIPTS}/../extern/${1} > /dev/null 2>&1

#Make symlinks if target does not exist / fix if broken
function ens_sym_link {
    local src=$1
    local tar=$2

    if [ -L ${tar} ] ; then
        if [ -e ${tar} ] ; then
            #Good link
            :
        else
            #Broken link - recreate
            ln -sf ${src} ${tar}
        fi
    elif [ -e ${tar} ] ; then
        #Not a link - presumably extant data
        :
    else
        #Missing, make a link
        ln -s ${src} ${tar}
    fi
}


ens_sym_link ${3} ../extern/${1}/${2}

###############################################################################

# >:D
