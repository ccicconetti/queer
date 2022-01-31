# quantum-routing
Quantum routing simulations

## Building instructions

Dependencies:

- cmake >= 3.2
- recent C++ compiler (GNU gcc >= 7 or LLVM clang >= 10)
- Google's glog library
- non-ancient Boost libraries

Clone repository:

```
git clone https://github.com/ccicconetti/quantum-routing.git
cd quantum-routing
git submodule update --init --recursive
```

To build without optimizations and with debug symbols:

```
mkdir debug
cd debug
cmake -DCMAKE_BUILD_TYPE=debug ../
```

To build with optimizations and without debug symbols (does not compile unit tests):

```
mkdir release
cd release
cmake -DCMAKE_BUILD_TYPE=release ../
```

To run the unit tests:

```
build/Test/testqr
```
