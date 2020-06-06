#!/bin/sh
## -e : Make sure all errors cause the script to fail
## -x be verbose; write what we are doing, as we do it
set -ex

# freebsd-like: systems based on sh/tcsh (not bash) with ICB
#   note: when you PR, docker-cp provides, in the container, the branch associated with the PR (not master where there's nothing new)
#         1. docker create --name mobydick IMAGE CMD        <=> create a container (= instance of image) but container is NOT yet started
#         2. docker cp -a ${TRAVIS_BUILD_DIR} mobydick:/tmp <=> copy git repository (CI worker, checkout-ed on PR branch) into the container
#                                                               note: docker-cp works only if copy from/to containers (not images)
#         3. docker start -a mobydick                       <=> start to run the container (initialized with docker-cp)
sudo docker pull ubuntu                                                                                           \
&&                                                                                                                \
sudo docker create --name mobydick ubuntu /bin/"$1" -c                                                            \
"apt-get -y update                                                                                             && \
 apt-get -y install build-essential diffutils findutils                                                        && \
 apt-get -y install git gfortran gcc g++ openmpi-bin libopenmpi-dev automake autoconf libtool pkg-config       && \
 apt-get -y install libblas-dev liblapack-dev                                                                  && \
 apt-get -y install libeigen3-dev                                                                              && \
 cd /tmp                                                                                                       && \
 cd arpack-ng                                                                                                  && \
 git status                                                                                                    && \
 git log -2                                                                                                    && \
 sed -e 's/LOG_FLAGS = /LOG_FLAGS = --allow-run-as-root --oversubscribe /' -i PARPACK/EXAMPLES/MPI/Makefile.am && \
 sed -e 's/LOG_FLAGS = /LOG_FLAGS = --allow-run-as-root --oversubscribe /' -i PARPACK/TESTS/MPI/Makefile.am    && \
 ./bootstrap                                                                                                   && \
 ./configure --enable-mpi --enable-icb-exmm --disable-dependency-tracking                                      && \
 export VERBOSE=1                                                                                              && \
 make all                                                                                                      && \
 make check                                                                                                    && \
 find . -name test-suite.log | xargs tail -n 300"                                                                 \
&&                                                                                                                \
sudo docker cp -a ${TRAVIS_BUILD_DIR} mobydick:/tmp                                                               \
&&                                                                                                                \
sudo docker start -a mobydick
