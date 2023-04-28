#!/bin/bash

set -e

print_msg () {
  echo "[*]" "$1"
}

###################### parallel configuration ######################

parallel=false

while [[ -n $# ]]; do
  case $1 in
    -p|--parallel)
      parallel=true
      ;;
    -*)
      echo "Unknown parameter: $1"
      exit 1
      ;;
    *)
      break
      ;;
  esac;
  shift
done

###################### building static libs ######################

print_msg "Compiling libmcl_rust"
if [[ "$parallel" = true ]]; then
  print_msg "Using parallel version"
  cargo rustc --release --crate-type=staticlib --features=parallel
else
  print_msg "Using non-parallel version"
  cargo rustc --release --crate-type=staticlib
fi

###################### cloning c-kzg-4844 ######################

print_msg "Cloning c-kzg-4844"
git clone https://github.com/ethereum/c-kzg-4844.git
cd c-kzg-4844 || exit
git -c advice.detachedHead=false checkout "$C_KZG_4844_GIT_HASH"
git submodule update --init

print_msg "Applying patches and building blst"
cd src
export CFLAGS="-Ofast -fno-builtin-memcpy -fPIC -Wall -Wextra -Werror"
make blst
unset CFLAGS
cd ..

###################### detecting os ######################

case $(uname -s) in
  "Linux")
    CSHARP_PLATFORM=linux-x64
    CLANG_PLATFORM=x86_64-linux
    ;;
  "Darwin")
    CSHARP_PLATFORM=osx-x64
    CLANG_PLATFORM=x86_64-darwin
    ;;
  *)
    echo "FAIL: unsupported OS"
    exit 1
    ;;
esac

###################### dotnet tests ######################

print_msg "Patching dotnet binding"
git apply < ../csharp.patch
cd bindings/csharp || exit

print_msg "Building dotnet"
make -B ckzg CSHARP_PLATFORM=$CSHARP_PLATFORM CLANG_PLATFORM=$CLANG_PLATFORM
dotnet restore

print_msg "Running dotnet tests"
dotnet test -c release --no-restore
cd ../..

###################### rust tests ######################

#print_msg "Patching rust binding"
#git apply < ../rust.patch
#cd bindings/rust || exit

#print_msg "Running rust tests"
#cargo test --release
#cd ../..

###################### python tests ######################

print_msg "Patching python binding"
git apply < ../python.patch
cd bindings/python || exit

print_msg "Running python tests"
make
cd ../..

###################### java tests ######################

print_msg "Patching java binding"
git apply < ../java.patch
cd bindings/java || exit

print_msg "Running java tests"
make CC_FLAGS=-lstdc++ build test
cd ../..

###################### nodejs tests ######################

print_msg "Patching nodejs binding"
git apply < ../nodejs.patch
cd bindings/node.js || exit

print_msg "Running nodejs tests"
make
cd ../../..

###################### go tests ######################

#print_msg "Patching go binding"
#git apply < ../go.patch
#cd bindings/go || exit

#print_msg "Running go tests"
#CGO_CFLAGS="-O2 -D__BLST_PORTABLE__" go test
#cd ../../..

###################### cleaning up ######################

print_msg "Cleaning up"
rm -rf c-kzg-4844
