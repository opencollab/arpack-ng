#!/bin/sh
## -e : Make sure all errors cause the script to fail
## -x be verbose; write what we are doing, as we do it
set -ex

# for non-LTS && EOL'ed ubuntu, the mirror url must be modified to run `apt update`
cmd=ls
if [ $2 = ":eoan" ]; then
    cmd="sed -i 's/\(security\|archive\).ubuntu/old-releases.ubuntu/g' /etc/apt/sources.list"
fi 

sudo docker pull "$1$2"                                                                                           \
&&                                                                                                                \
sudo docker create --name mobydick "$1$2" /bin/bash -c                                                            \
"${cmd} && \
 apt-get    update                                                                                             && \
 ln -snf /usr/share/zoneinfo/Europe/Paris /etc/localtime && echo 'Europe/Paris' > /etc/timezone                && \
 apt-get -y install build-essential                                                                            && \
 apt-get -y install dialog apt-utils                                                                           && \
 apt-get -y install git gfortran gcc g++ openmpi-bin libopenmpi-dev automake autoconf libtool pkg-config cmake && \
 apt-get -y install libblas-dev liblapack-dev                                                                  && \
 cd /tmp                                                                                                       && \
 cd arpack-ng                                                                                                  && \
 git status                                                                                                    && \
 git log -2                                                                                                    && \
 sed -e 's/LOG_FLAGS = /LOG_FLAGS = --allow-run-as-root --oversubscribe /' -i PARPACK/EXAMPLES/MPI/Makefile.am && \
 sed -e 's/LOG_FLAGS = /LOG_FLAGS = --allow-run-as-root --oversubscribe /' -i PARPACK/TESTS/MPI/Makefile.am    && \
 ./bootstrap                                                                                                   && \
 ./configure --enable-icb --enable-mpi --disable-dependency-tracking                                           && \
 export VERBOSE=1                                                                                              && \
 make all                                                                                                      && \
 make check; find . -name test-suite.log | xargs tail -n 50                                                    && \
 make distcheck                                                                                                && \
 sed -e 's/mpirun /mpirun --allow-run-as-root --oversubscribe /' -i CMakeLists.txt                             && \
 mkdir -p build && cd build                                                                                    && \
 cmake -DEXAMPLES=ON -DMPI=ON -DICB=ON ..                                                                      && \
 make all                                                                                                      && \
 make test; tail -n 50 ./Testing/Temporary/LastTest.log"                                                          \
&&                                                                                                                \
sudo docker cp -a ${GITHUB_WORKSPACE} mobydick:/tmp                                                               \
&&                                                                                                                \
sudo docker start -a mobydick
