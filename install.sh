build_dir=~/temp/easifem-extpkgs/arpack-ng/build

cmake -S ./ -B ${build_dir} -D CMAKE_INSTALL_PREFIX:PATH=${EASIFEM_EXTPKGS} -D BUILD_SHARED_LIBS=ON -D MPI=ON -D CMAKE_BUILD_TYPE:STRING=Release

cmake --build ${build_dir} --target install