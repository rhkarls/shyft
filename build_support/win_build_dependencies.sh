#!/bin/bash
export WORKSPACE=$(readlink --canonicalize --no-newline `dirname ${0}`/../..)
# to align the cmake support:
export SHYFT_DEPENDENCIES_DIR=${WORKSPACE}/shyft_dependencies
armadillo_name=armadillo-8.400.0
dlib_name=dlib-19.10
boost_ver=1_66_0
numpy_ver=1.13
cmake_common="-DCMAKE_INSTALL_MESSAGE=NEVER"
echo ---------------
echo Windows Update/build shyft-dependencies
echo WORKSPACE..............: ${WORKSPACE}
echo SHYFT_DEPENDENCIES_DIR.: ${SHYFT_DEPENDENCIES_DIR}
echo PACKAGES...............: miniconda w/shyft_env, doctest, boost_${boost_ver}, ${armadillo_name}, ${dlib_name}, numpy=${numpy_ver} 
WGET='curl -L -O'

# A helper function to compare versions
function version { echo "$@" | awk -F. '{ printf("%d%03d%03d%03d\n", $1,$2,$3,$4); }'; }

# the current versions we are building
mkdir -p ${SHYFT_DEPENDENCIES_DIR}
cd ${SHYFT_DEPENDENCIES_DIR}

if [ ! -d ${armadillo_name} ]; then 
    echo Building ${armadillo_name}
    if [ ! -f ${armadillo_name}.tar.xz ]; then 
        ${WGET}  http://sourceforge.net/projects/arma/files/${armadillo_name}.tar.xz
    fi;
    7z -y e ${armadillo_name}.tar.xz
    7z -y x ${armadillo_name}.tar
    pushd ${armadillo_name}
    mkdir -p build
    cd build && cmake -G"Visual Studio 15 2017 Win64" -DCMAKE_INSTALL_PREFIX=${SHYFT_DEPENDENCIES_DIR} -DCMAKE_INSTALL_LIBDIR=lib -DDETECT_HDF5=0 -DARMA_USE_WRAPPER=FALSE -DARMA_USE_LAPACK=TRUE -DARMA_USE_BLAS=TRUE ${cmake_common} ..
    cmake --build . --config Release 
    cmake -P cmake_install.cmake
	cp ../examples/lib_win64/* ${SHYFT_DEPENDENCIES_DIR}/lib
    popd
fi;
echo Done ${armadillo_name}


if [ ! -d ${dlib_name} ]; then
    echo Building ${dlib_name}
    if [ ! -f ${dlib_name}.tar.bz2 ]; then
        ${WGET} http://dlib.net/files/${dlib_name}.tar.bz2
    fi;
    7z -y e ${dlib_name}.tar.bz2
    7z -y x ${dlib_name}.tar
    pushd ${dlib_name}
    mkdir -p build
    dlib_cfg="-DDLIB_PNG_SUPPORT=0 -DDLIB_GIF_SUPPORT=0 -DDLIB_LINK_WITH_SQLITE3=0 -DDLIB_NO_GUI_SUPPORT=1 -DDLIB_DISABLE_ASSERTS=1 -DDLIB_JPEG_SUPPORT=0 -DDLIB_USE_BLAS=0 -DDLIB_USE_LAPACK=0 -DBUILD_SHARED_LIBS=0"
    cd build
	cmake -G"Visual Studio 15 2017 Win64" .. -DCMAKE_INSTALL_PREFIX=${SHYFT_DEPENDENCIES_DIR} -DCMAKE_INSTALL_LIBDIR=lib ${cmake_common} ${dlib_cfg} 
	cmake --build . --config Release --target install
	cmake --build . --config Debug --target install
    popd
fi;
echo Done ${dlib_name}

if [ ! -d doctest ]; then
    echo Building doctest
    git clone https://github.com/onqtam/doctest
    pushd doctest
    cmake -G"Visual Studio 15 2017 Win64" -DCMAKE_INSTALL_PREFIX=${SHYFT_DEPENDENCIES_DIR} ${cmake_common} && cmake -P cmake_install.cmake
    popd
fi;
echo Done doctest

cd ${WORKSPACE}
type python
python_tst=$?
if [ ! ${python_tst} -eq 0 ]; then
	if [ ! -d miniconda/Scripts ]; then
		echo Missing python install. try to make one using miniconda
		if [ -d miniconda ]; then
			rm -rf miniconda
		fi;
		if [ ! -f Miniconda3-latest-Windows-x86_64.exe ]; then
			${WGET}  http://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86_64.exe
		fi;
		echo 'start /wait "" .\Miniconda3-latest-Windows-x86_64.exe /InstallationType=JustMe /S /D=%cd%\miniconda' >install_miniconda.cmd
		./install_miniconda.cmd
		# Update conda to latest version, assume we start with 4.3 which
		# requires PATH to be set
		OLDPATH=${PATH}
		export PATH="${WORKSPACE}/miniconda:${WORKSPACE}/miniconda/Scripts:${WORKSPACE}/miniconda/Dlls:$PATH"

		old_conda_version=$(conda --version | sed "s/conda \(.*\)/\1/")
		echo "Old conda version is ${old_conda_version}"
		#activate
		conda config --set always_yes yes --set changeps1 no
		conda update conda
		new_conda_version=$(conda --version | sed "s/conda \(.*\)/\1/")
		echo "New conda version is ${new_conda_version}"
		conda install numpy=${numpy_ver}
		conda create -n shyft_env python=3.6 pyyaml numpy=${numpy_ver} netcdf4 gdal matplotlib requests nose coverage pip shapely  pyproj
	else
		export PATH="${WORKSPACE}/miniconda:${WORKSPACE}/miniconda/Scripts:${WORKSPACE}/miniconda/Dlls:$PATH"
	fi;
else
	echo using pre-installed `type -p python`
fi;
echo Done minconda/python
export BOOST_PYTHONHOME=`type -p python |  sed  -e 's_/python__' -e 's/^\///' -e 's_/_\\\\_g' -e 's/^./\0:/'`
echo Setting BOOST_PYTHONHOME to  ${BOOST_PYTHONHOME}
cd ${SHYFT_DEPENDENCIES_DIR}
if [ ! -d boost_${boost_ver} ]; then
    echo Building boost_${boost_ver}
    if [ ! -f boost_${boost_ver}.tar.gz ]; then
        ${WGET} http://sourceforge.net/projects/boost/files/boost/${boost_ver//_/.}/boost_${boost_ver}.tar.gz
    fi;
    7z -y e boost_${boost_ver}.tar.gz
    7z -y x boost_${boost_ver}.tar
    pushd boost_${boost_ver}
    echo "set PYTHONHOME=%BOOST_PYTHONHOME%" >bboost.cmd
    echo "call bootstrap.bat" >> bboost.cmd
    echo "b2 -d0 -j 6 define=BOOST_CONFIG_SUPPRESS_OUTDATED_MESSAGE link=shared variant=release,debug threading=multi runtime-link=shared address-model=64 --with-system --with-filesystem --with-date_time --with-serialization --with-python --prefix=%cd%\.. install" >> bboost.cmd
    ./bboost.cmd
    popd
fi;
echo  Done boost_${boost_ver}

cd ${SHYFT_DEPENDENCIES_DIR}
if [ ! -d pybind11 ]; then
    git clone https://github.com/pybind/pybind11.git
	pushd pybind11
	mkdir -p build
	cd build
	echo set PYTHONHOME=%BOOST_PYTHONHOME% >bdlib.cmd
	echo cmake -G\"Visual Studio 15 2017 Win64\" -DCMAKE_INSTALL_PREFIX=${SHYFT_DEPENDENCIES_DIR} -DPYBIND11_TEST=0 ${cmake_common} .. >>bdlib.cmd
	echo cmake -P cmake_install.cmake >>bdlib.cmd
	./bdlib.cmd
	popd
fi;
echo Done pybind11

cd ${WORKSPACE}
if [ -d shyft-data ]; then 
    pushd shyft-data
    git pull >/dev/null
    popd
else 
    git clone https://github.com/statkraft/shyft-data
fi;
echo Done shyft-data

