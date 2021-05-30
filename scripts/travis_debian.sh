#!/bin/sh
## -e : Make sure all errors cause the script to fail
## -x be verbose: write what we are doing, as we do it
set -ex

# This script is debian-based: enable to run debian, ubuntu (debian-based) and linuxmint (ubuntu-based).

# Need to build boost-python from source for python 3 (repository package is built for python 2).

export BTVER="67"
INSTALL_BP=$(echo "\
 cd /tmp                                                                                                       && \
 apt-get -y install python3-minimal python3-pip python3-numpy                                                  && \
 pip3 install --system numpy                                                                                   && \
 apt-get -y install wget                                                                                       && \
 wget -q https://sourceforge.net/projects/boost/files/boost/1.${BTVER}.0/boost_1_${BTVER}_0.tar.gz             && \
 tar -xf boost_1_${BTVER}_0.tar.gz && cd boost_1_${BTVER}_0                                                    && \
 ./bootstrap.sh --with-libraries=python --with-python=/usr/bin/python3 --with-toolset=gcc >/dev/null 2>&1      && \
 ./b2 install toolset=gcc                                                                 >/dev/null 2>&1      && \
 apt-get install locate                                                                                        && \
 updatedb                                                                                                      && \
 export BOOST_ROOT=/usr/local                                                                                  && \
 export Boost_NO_SYSTEM_PATHS=ON                                                                               && \
 ")
export INSTALL_BP

# Testing boost-python specifying the correct python version.

export BPVER="38"
TEST_BP=$(echo "-DPYTHON3=ON -DBOOST_PYTHON_LIBSUFFIX=${BPVER}")
export TEST_BP

# Testing ILP64.

export INSTALL_BLAS_LAPACK="libblas-dev liblapack-dev"

export INSTALL_MKL_ILP64="\
 cd /tmp                                                                                                       && \
 export DEBIAN_FRONTEND=noninteractive                                                                         && \
 apt-get -y install dialog                                                                                     && \
 echo yes | apt-get -y install intel-mkl libmkl-dev                                                            && \
 export FFLAGS='-DMKL_ILP64 -I/usr/include/mkl'                                                                && \
 export FCFLAGS='-DMKL_ILP64 -I/usr/include/mkl'                                                               && \
 export LIBS='-Wl,--no-as-needed -L/usr/lib/x86_64-linux-gnu -lmkl_sequential -lmkl_core -lpthread -lm -ldl'   && \
 export INTERFACE64=1                                                                                          && \
"

export TEST_MKL_ILP64="--with-blas=mkl_gf_ilp64 --with-lapack=mkl_gf_ilp64"

# Selecting what to test.

export TEST_DISTCHECK=""
if [ ".$3" = .REG ] # Regular: no boost, no ILP64.
then
    export INSTALL_BP=""
    export TEST_BP=""
    export INSTALL_MKL_ILP64=""
    export TEST_MKL_ILP64=""
    export TEST_DISTCHECK="make distcheck &&"
fi
if [ ".$3" = .BP ] # Boost Python only.
then
    export INSTALL_MKL_ILP64=""
    export TEST_MKL_ILP64=""
fi
if [ ".$3" = .ILP64 ] # ILP64 only.
then
    export INSTALL_BP=""
    export TEST_BP=""
    export INSTALL_BLAS_LAPACK=""
fi

# For non-LTS && EOL'ed ubuntu, the mirror url must be modified to run `apt update`
cmd=ls
if [ $2 = ":eoan" ]; then
    cmd="sed -i 's/\(security\|archive\).ubuntu/old-releases.ubuntu/g' /etc/apt/sources.list"
fi

# Create command to run.

CMD=$(echo "\
 ${cmd}                                                                                                        && \
 sed -e 's/main/main non-free contrib/' -i /etc/apt/sources.list                                               && \
 apt-get    update                                                                                             && \
 ln -snf /usr/share/zoneinfo/Europe/Paris /etc/localtime && echo 'Europe/Paris' > /etc/timezone                && \
 apt-get -y --allow-unauthenticated -o Dpkg::Options::=--force-confdef      upgrade                            && \
 apt-get -y --allow-unauthenticated -o Dpkg::Options::=--force-confdef dist-upgrade                            && \
 apt-get -y install build-essential                                                                            && \
 apt-get -y install dialog apt-utils                                                                           && \
 apt-get -y install git gfortran gcc g++ openmpi-bin libopenmpi-dev automake autoconf libtool pkg-config cmake && \
 apt-get -y install ${INSTALL_BLAS_LAPACK} libeigen3-dev                                                       && \
 ${INSTALL_MKL_ILP64}                                                                                             \
 cd /tmp                                                                                                       && \
 cd arpack-ng                                                                                                  && \
 git status                                                                                                    && \
 git log -2                                                                                                    && \
 sed -e 's/LOG_FLAGS = /LOG_FLAGS = --allow-run-as-root --oversubscribe /' -i PARPACK/EXAMPLES/MPI/Makefile.am && \
 sed -e 's/LOG_FLAGS = /LOG_FLAGS = --allow-run-as-root --oversubscribe /' -i PARPACK/TESTS/MPI/Makefile.am    && \
 ./bootstrap                                                                                                   && \
 ./configure ${TEST_MKL_ILP64} --enable-icb --enable-mpi --disable-dependency-tracking                         && \
 export VERBOSE=1                                                                                              && \
 make all                                                                                                      && \
 make check && find . -name test-suite.log | xargs tail -n 50                                                  && \
 ${TEST_DISTCHECK}                                                                                                \
 ${INSTALL_BP}                                                                                                    \
 cd /tmp                                                                                                       && \
 cd arpack-ng                                                                                                  && \
 git status                                                                                                    && \
 git log -2                                                                                                    && \
 sed -e 's/mpirun /mpirun --allow-run-as-root --oversubscribe --mca btl self /' -i CMakeLists.txt              && \
 mkdir -p build && cd build                                                                                    && \
 cmake -DEXAMPLES=ON -DMPI=ON -DICB=ON -DICBEXMM=ON ${TEST_BP} ..                                              && \
 make all                                                                                                      && \
 make test && tail -n 50 ./Testing/Temporary/LastTest.log                                                         \
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
