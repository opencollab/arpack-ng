#!/bin/sh
## -e : Make sure all errors cause the script to fail
## -x be verbose; write what we are doing, as we do it
set -ex

# Need to build boost-python from source for python 3 (repository package is built for python 2).

export BTVER="73"
INSTALL_BP=$(echo "\
 cd /tmp                                                                                                       && \
 apt-get -y install python3-minimal python3-pip python3-numpy                                                  && \
 pip3 install --system numpy                                                                                   && \
 apt-get -y install wget                                                                                       && \
 wget -q https://sourceforge.net/projects/boost/files/boost/1.${BTVER}.0/boost_1_${BTVER}_0.tar.gz             && \
 tar -xf boost_1_${BTVER}_0.tar.gz && cd boost_1_${BTVER}_0                                                    && \
 ./bootstrap.sh --with-libraries=python --with-python=/usr/bin/python3 --with-toolset=gcc >/dev/null 2>&1      && \
 ./b2 install                                                                             >/dev/null 2>&1      && \
 apt-get install locate                                                                                        && \
 updatedb                                                                                                      && \
 export BOOST_ROOT=/usr/local                                                                                  && \
 export Boost_NO_SYSTEM_PATHS=ON                                                                                ; \
 ")
export INSTALL_BP

# Testing boost-python specifying the correct python version.

export BPVER="38"
TEST_BP=$(echo "-DPYTHON3=ON -DBOOST_PYTHON_LIBSUFFIX=${BPVER}")
export TEST_BP

# Create command to run.

CMD=$(echo "\
 apt-get    update                                                                                             && \
 ln -snf /usr/share/zoneinfo/Europe/Paris /etc/localtime && echo 'Europe/Paris' > /etc/timezone                && \
 apt-get -y install build-essential                                                                            && \
 apt-get -y install dialog apt-utils                                                                           && \
 apt-get -y install git gfortran gcc g++ openmpi-bin libopenmpi-dev automake autoconf libtool pkg-config cmake && \
 apt-get -y install libblas-dev liblapack-dev libeigen3-dev                                                    && \
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
 make distcheck                                                                                                 ; \
 ${INSTALL_BP}                                                                                                    \
 cd /tmp                                                                                                       && \
 cd arpack-ng                                                                                                  && \
 git status                                                                                                    && \
 git log -2                                                                                                    && \
 sed -e 's/mpirun /mpirun --allow-run-as-root --oversubscribe /' -i CMakeLists.txt                             && \
 mkdir -p build && cd build                                                                                    && \
 cmake -DEXAMPLES=ON -DMPI=ON -DICB=ON -DICBEXMM=ON ${TEST_BP} ..                                              && \
 make all                                                                                                      && \
 make test; tail -n 50 ./Testing/Temporary/LastTest.log                                                         ; \
")
export CMD

# Run command over docker.

sudo docker pull "$1$2"                                                                                           \
&&                                                                                                                \
sudo docker create --name mobydick "$1$2" /bin/bash -c "${CMD}"                                                   \
&&                                                                                                                \
sudo docker cp -a "${TRAVIS_BUILD_DIR}" mobydick:/tmp                                                             \
&&                                                                                                                \
sudo docker start -a mobydick
