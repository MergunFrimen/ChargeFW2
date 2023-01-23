<h1 align="center">
  <!-- <a href=""> -->
  <img src="./doc/logo.png" alt="ChargeFW2" width="500">
  <!-- </a> -->
</h1>


<p align="center">
  <a href="#compilation-requirements">Compilation requirements</a> •
  <a href="#installation">Installation</a> •
  <a href="#docker">Docker</a> •
  <a href="#usage">Usage</a> •
  <a href="#how-to-cite">How To Cite</a> 
</p>

<h4 align="center">Application for computing partial atomic charges using selected empirical methods. ChargeFW2 is the computational core of <a href="https://acc2.ncbr.muni.cz" target="_blank">ACC2</a>.</h4>

See the [short description](https://acc2.ncbr.muni.cz/static/methods.pdf) of implemented methods. 

## Compilation requirements
- [CMake](https://cmake.org/) 3.17
- [GCC](https://gcc.gnu.org/) 10 or [Clang](https://clang.llvm.org/) 10
- [Boost](https://www.boost.org/) 1.69
- [Eigen](http://eigen.tuxfamily.org) 3.3
- [fmt](https://fmt.dev) 6.2.1
- [nanoflann](https://github.com/jlblancoc/nanoflann) 1.4.3
- [JSON for Modern C++](https://github.com/nlohmann/json) 3.7.3
- [GEMMI](https://github.com/project-gemmi/gemmi) 0.4.7
- [pybind11](https://github.com/pybind/pybind11) 2.5.0

> **Note:**
> Tested on Fedora 32-36 and Ubuntu 20.04-22.04. Other version of the libraries might work too however this was not tested.

## Installation
After downloading and unpacking the sources, run the following in the ChargeFW2 directory:

```shell script
$ mkdir build
$ cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=<WHERE-TO-INSTALL> -DCMAKE_BUILD_TYPE=Release
$ make
$ make install
```

## Docker
Rather than installing all dependencies, you can run ChargeFW2 directly in a Docker container:

```shell script
$ docker build -t chargefw2 .
$ docker run -it --rm --entrypoint bash chargefw2
```

A prebuild image is available on [Docker Hub](https://hub.docker.com/r/frimen/chargefw2) which can be used directly:

```shell script
$ docker run -it --rm --entrypoint bash docker.io/frimen/chargefw2
```

The Docker container is setup for use in CLI workflows. Example using relative paths to files:

```shell script
$ docker build -t chargefw2 .
$ docker run -it --rm -v $PWD:$PWD chargefw2 --mode charges --input-file $PWD/doc/molecules.sdf --chg-out-dir $PWD/
```

There is also a devcontainer available for this project. You can read more about it [here](.devcontainer/README.md).

## Usage

The [documentation](doc/documentation.md) for the application and its [Python bindings](doc/ChargeFW2%20-%20tutorial.pdf) is located in the [doc](doc) folder.

## How to cite
If you found ChargeFW2 or Atomic Charge Calculator II helpful, please cite: [Raček, T., Schindler, O., Toušek, D., Horský, V., Berka, K., Koča, J., & Svobodová, R. (2020). Atomic Charge Calculator II: web-based tool for the calculation of partial atomic charges. Nucleic Acids Research](https://doi.org/10.1093/nar/gkaa367).
