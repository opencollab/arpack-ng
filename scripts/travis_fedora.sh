#!/bin/sh
## -e : Make sure all errors cause the script to fail
## -x be verbose; write what we are doing, as we do it
set -ex
## Should we init a container?
if [ ".$1" = .setup ]
then
  # fedora
  #   note: when you PR, docker-cp provides, in the container, the branch associated with the PR (not master where there's nothing new)
  #         1. docker create --name mobydick IMAGE CMD        <=> create a container (= instance of image) but container is NOT yet started
  #         2. docker cp -a ${GITHUB_WORKSPACE} mobydick:/tmp <=> copy git repository (CI worker, checkout-ed on PR branch) into the container
  #                                                               note: docker-cp works only if copy from/to containers (not images)
  #         3. docker start -a mobydick                       <=> start to run the container (initialized with docker-cp)
    test . != ".$2" && mpi="$2" || mpi=openmpi
    time sudo docker pull fedora
    time sudo docker create --name mobydick fedora /tmp/arpack-ng/scripts/travis_fedora.sh $mpi
    time sudo docker cp -a ${GITHUB_WORKSPACE} mobydick:/tmp
    time sudo docker start -a mobydick ; e=$?
    exit $e
fi

test . != ".$1" && mpi="$1" || mpi=openmpi

## If we are called as root, setup everything
if [ $UID -eq 0 ]
then
    cat /etc/os-release
    # Ignore weak depencies
    echo "install_weak_deps=False" >> /etc/dnf/dnf.conf
    time dnf -y upgrade
    time dnf -y install environment-modules git gfortran openblas-devel cmake ${mpi}-devel make gcc-c++
    useradd test
    chown -R test /tmp
    sudo -u test $0 $mpi
## If we are called as normal user, run test
else
    . /etc/profile.d/modules.sh
    module load mpi
    export OMPI_MCA_rmaps_base_oversubscribe=yes
    cd /tmp
    cd arpack-ng
    git status
    git log -2
    mkdir -p build && cd build
    time cmake -DEXAMPLES=ON -DMPI=ON -DICB=ON ..
    export VERBOSE=1
    time make all
    time make test
    tail -n 300 ./Testing/Temporary/LastTest.log
fi
