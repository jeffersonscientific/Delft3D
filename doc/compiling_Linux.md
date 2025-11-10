# Compiling Delft3D on Linux
For the latest overview of build steps see the Docker files located in `ci\dockerfiles\linux\`.

## Prerequisites
- Various basic Linux utilities
```
dnf update --assumeyes
dnf install --assumeyes epel-release
dnf config-manager --set-enabled powertools

dnf install --assumeyes \
    which binutils patchelf diffutils procps m4 make gcc gcc-c++ \
    openssl openssl-devel wget perl python3 xz curl-devel
```

- [Intel oneAPI Fortran Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler-download.html)
- [Intel oneAPI MPI Library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html)
- [Intel oneAPI Math Kernal Library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html)
```
INTEL_ONEAPI_VERSION=2024
COMMON_VARS_VERSION="2024.2.1"
COMPILER_DPCPP_CPP_VERSION="2024.2.1"
COMPILER_FORTRAN_VERSION="2024.2.1"
MKL_DEVEL_VERSION="2024.2.2"
MPI_DEVEL_VERSION="2021.13.1"

cat <<EOT > /etc/yum.repos.d/oneAPI.repo
[oneAPI]
name=Intel® oneAPI repository
baseurl=https://yum.repos.intel.com/oneapi
enabled=1
gpgcheck=1
repo_gpgcheck=1
gpgkey=https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
EOT

dnf install --assumeyes \
    intel-oneapi-common-vars-${COMMON_VARS_VERSION} \
    intel-oneapi-compiler-dpcpp-cpp-${COMPILER_DPCPP_CPP_VERSION} \
    intel-oneapi-compiler-fortran-${COMPILER_FORTRAN_VERSION} \
    intel-oneapi-mkl-devel-${MKL_DEVEL_VERSION} \
    intel-oneapi-mpi-devel-${MPI_DEVEL_VERSION}


```

- autoconf
- automake
- libtool
```
source /opt/intel/oneapi/setvars.sh

for URL in \
    'https://mirrors.kernel.org/gnu/autoconf/autoconf-2.72.tar.xz' \
    'https://mirrors.kernel.org/gnu/automake/automake-1.17.tar.xz' \
    'https://mirrors.kernel.org/gnu/libtool/libtool-2.4.7.tar.xz'
do
    BASEDIR=$(basename -s '.tar.xz' "$URL")
    if [[ -d "/var/cache/src/${BASEDIR}" ]]; then
        echo "CACHED ${BASEDIR}"
    else
        echo "Fetching ${BASEDIR}.tar.xz..."
        wget --quiet --output-document=- "$URL" | tar --extract --xz --file=- --directory='/var/cache/src/'
    fi

    pushd "/var/cache/src/${BASEDIR}"
    ./configure CC=icx CXX=icpx FC=ifx CFLAGS="-O3" CXXFLAGS="-O3" FCFLAGS="-O3"
    make --jobs=$(nproc)
    make install
    popd
done
```

- Ninja is required for buiding cmake
```
wget https://github.com/ninja-build/ninja/releases/download/v1.12.1/ninja-linux.zip
unzip ninja-linux.zip
chmod +x ninja
mv ninja /usr/bin/
rm ninja-linux.zip

echo "Installed ninja version:" $(ninja --version)
```

- [CMake](https://cmake.org/download/) version 3.30 or later
```
source /opt/intel/oneapi/setvars.sh

URL='https://github.com/Kitware/CMake/releases/download/v3.30.3/cmake-3.30.3.tar.gz'
BASEDIR=$(basename -s '.tar.gz' "$URL")
if [[ -d "/var/cache/src/${BASEDIR}" ]]; then
    echo "CACHED ${BASEDIR}"
else
    echo "Fetching ${BASEDIR}.tar.gz..."
    wget --quiet --output-document=- "$URL" | tar --extract --gzip --file=- --directory='/var/cache/src'
fi

export CC=icx CXX=icpx CFLAGS="-O3" CXXFLAGS="-O3"

pushd /var/cache/src/cmake-3.30.3
./bootstrap --parallel=$(nproc)
make --jobs=$(nproc)
make install
popd
```

## Build steps for third party libraries
- 

## Build steps
- build.sh
  Execute "./build.sh --help" to show the usage
  Currently used as default build process: "./build.sh fm-suite --compiler intel21"
  This will execute "src/setenv.sh" on Deltares systems. On other systems, the environment must be prepared upfront.
  For instructions, see [Setup your own Linux environment](Linux_setup.md).

## Alternative: without build-script
Refer to the [README](src\cmake\README) file in the CMake folder.
WARNING: When building without our build-script, the collection of the resulting binaries will need attention
