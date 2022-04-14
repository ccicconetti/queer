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

## Execute experiments

1. build the software (see instructions), e.g., assume you build in `release/`
2. enter the experiment you want to run, e.g., `Experiments/001_Constant_Rate_PPP`
3. create a symbolic link from the corresponding executable in your build directory (e.g., `release/Experiments/001_Constant_Rate_PPP/main-001`) to the experiment directory
4. enter into the sub-experiment directory, e.g., `var-flows`
5. execute the script `./pre.sh` (if present): this will satisfy pre-run requirements, such as downloading external datasets
6. execute the script `./run.sh`: this will populate a directory called `data` with CSV output, one per batch of experiments; the meaning of the column can be retrieved by running the executable with `--explain-output`
7. execute the script `./post.sh` (if present): this will do some post-processing analysis on the results, whose output will be stored in `post`; note that you will need a valid Python2 interpreter and Internet access, which is required to download the utility Python script [percentile.py](https://raw.githubusercontent.com/ccicconetti/serverlessonedge/master/scripts/percentile.py) from GitHub. If the machine running the post-processing does not have Internet access, you may download the file and copy it into the sub-experiment directory with exec permissions
8. there may be some [Gnuplot](http://www.gnuplot.info/) scripts in the directory `graph`, you can run them by calling `gnuplot -persists SCRIPT.plt`

Full example, assuming you build in `release` and you have a working Gnuplot:

```
cd Experiments/001_Constant_Rate_PPP
ln -s ../../release/Experiments/001_Constant_Rate_PPP/main-001 .
cd var-flows
./run.sh
./post.sh
cd graph
gnuplot -persist admission-cdf.plt
```

The last command will show a plot similar to this one:

![](Docs/var-flows-admission-rate.png)
