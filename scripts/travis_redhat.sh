#!/bin/sh
# usage: ./travis_centos.sh [:7|:8|:latest|:stream]
#
## -e : Make sure all errors cause the script to fail
## -x be verbose; write what we are doing, as we do it
set -ex

# centos-stream should be pulled from quay.io
# note: a RedHad account is required after 2021.06,
#       see https://quay.io/repository/centos/centos?tab=tags
if [ $1 = ":stream" ]; then
    prefix="quay.io/centos/"
fi

if [ $1 = ":7" ]; then # for 7
    pm="yum"
    cmake="cmake3"
else # for 8+, stream
    pm="dnf"
    cmake="cmake"
    cmd_powertools="dnf config-manager --set-enabled powertools &&"
fi

sudo docker pull "${prefix}centos$1"                                     \
&&                                                                       \
sudo docker create --name mobydick ${prefix}centos$1 /bin/bash -c        \
"$pm install -y ${pm}-plugins-core epel-release                       && \
 $pm upgrade -y                                                       && \
 $cmd_powertools                                                         \
 $pm install -y git make gcc gcc-gfortran gcc-c++ environment-modules && \
 $pm install -y $cmake                                                && \
 $pm install -y mpich-devel                                           && \
 $pm --enablerepo=\"epel\" install -y openblas-devel lapack-devel     && \
 . /etc/profile.d/modules.sh                                          && \
 module avail && module load mpi && module list                       && \
 cd /tmp                                                              && \
 cd arpack-ng                                                         && \
 git status                                                           && \
 git log -2                                                           && \
 mkdir -p build && cd build                                           && \
 $cmake -DEXAMPLES=ON -DMPI=ON -DICB=ON ..                            && \
 make all && make test"                                                  \
&&                                                                       \
sudo docker cp -a ${GITHUB_WORKSPACE} mobydick:/tmp                      \
&&                                                                       \
sudo docker start -a mobydick

