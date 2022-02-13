#!/bin/sh
# usage: ./travis_centos.sh [:7|:8|:latest|:stream]
#
## -e : Make sure all errors cause the script to fail
## -x be verbose; write what we are doing, as we do it
set -ex

podman pull registry.access.redhat.com/ubi8/ubi                          \
&&                                                                       \
podman create --name ubi ${prefix}centos$1 /bin/bash -c                  \
"cat /etc/os-release                                                  && \
 dnf install -y dnf-plugins-core epel-release                         && \
 dnf upgrade -y                                                       && \
 dnf config-manager --set-enabled powertools                          && \
 dnf install -y git make gcc gcc-gfortran gcc-c++ environment-modules && \
 dnf install -y cmake                                                 && \
 dnf install -y mpich-devel                                           && \
 dnf --enablerepo=\"epel\" install -y openblas-devel lapack-devel     && \
 . /etc/profile.d/modules.sh                                          && \
 module avail && module load mpi && module list                       && \
 cd /tmp                                                              && \
 cd arpack-ng                                                         && \
 git status                                                           && \
 git log -2                                                           && \
 mkdir -p build && cd build                                           && \
 cmake -DEXAMPLES=ON -DMPI=ON -DICB=ON ..                             && \
 make all && make test"                                                  \
&&                                                                       \
podman cp -a ${GITHUB_WORKSPACE} ubi:/tmp                                \
&&                                                                       \
podman start -a ubi

