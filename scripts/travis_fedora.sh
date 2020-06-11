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
  #         2. docker cp -a ${TRAVIS_BUILD_DIR} mobydick:/tmp <=> copy git repository (CI worker, checkout-ed on PR branch) into the container
  #                                                               note: docker-cp works only if copy from/to containers (not images)
  #         3. docker start -a mobydick                       <=> start to run the container (initialized with docker-cp)
    test . != ".$2" && mpi="$2" || mpi=openmpi
    test . != ".$3" && os="$3" || os=fedora:latest
    base=${version%:*}
    release=${os#*:}
    case $base in
	fedora)
	    reg=registry.fedoraproject.org
	    ;;
	centos)
	    reg=registry.centos.org
	    ;;
	*)
	    echo "Unknonw OS"
	    exit 1
    esac
    time sudo docker pull $reg/$os ||
	sudo docker pull $os
    time sudo docker create --name mobydick $os \
	/tmp/arpack-ng/scripts/travis_fedora.sh $mpi
    time sudo docker cp -a ${TRAVIS_BUILD_DIR} mobydick:/tmp
    time sudo docker start -a mobydick ; e=$?
    exit $e
fi

test . != ".$1" && mpi="$1" || mpi=openmpi

## If we are called as root, setup everything
if [ $UID -eq 0 ]
then
    if grep centos -i /etc/os-release ; then
	dnf install -y dnf-plugins-core #epel-release
	dnf config-manager --set-enabled PowerTools
    fi
    time dnf -y upgrade
    time dnf -y install environment-modules git \
        gfortran openblas-devel cmake ${mpi}-devel make gcc-c++
    useradd test
    chown -R test /tmp
    sudo -u test $0 $mpi
## If we are called as normal user, run test
else
    on_err() {
	env || :
	tail -n 300 config.log || :
	tail -n 300 ./Testing/Temporary/LastTest.log || :
    }
    trap on_err ERR
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
