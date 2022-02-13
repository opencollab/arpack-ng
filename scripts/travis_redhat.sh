#!/bin/sh
# usage: ./travis_centos.sh [:7|:8|:latest|:stream]
#
## -e : Make sure all errors cause the script to fail
## -x be verbose; write what we are doing, as we do it
set -ex

pm="dnf"
cmake="cmake"
cmd_powertools="dnf config-manager --set-enabled powertools &&"

podman pull registry.access.redhat.com/ubi8/ubi                          \
&&                                                                       \
podman create --name ubi ${prefix}centos$1 /bin/bash -c                  \
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
podman cp -a ${GITHUB_WORKSPACE} ubi:/tmp                                \
&&                                                                       \
podman start -a ubi

