<div style='-webkit-text-size-adjust:100%;-webkit-tap-highlight-color:transparent;--blue:#007bff;--indigo:#6610f2;--purple:#6f42c1;--pink:#e83e8c;--red:#dc3545;--orange:#fd7e14;--yellow:#ffc107;--green:#28a745;--teal:#20c997;--cyan:#17a2b8;--white:#fff;--gray:#6c757d;--gray-dark:#343a40;--primary:#007bff;--secondary:#6c757d;--success:#28a745;--info:#17a2b8;--warning:#ffc107;--danger:#dc3545;--light:#f8f9fa;--dark:#343a40;--breakpoint-xs:0;--breakpoint-sm:576px;--breakpoint-md:768px;--breakpoint-lg:992px;--breakpoint-xl:1200px;--font-family-sans-serif:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Helvetica Neue",Arial,"Noto Sans",sans-serif,"Apple Color Emoji","Segoe UI Emoji","Segoe UI Symbol","Noto Color Emoji";--font-family-monospace:SFMono-Regular,Menlo,Monaco,Consolas,"Liberation Mono","Courier New",monospace;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Helvetica Neue",Arial,"Noto Sans",sans-serif,"Apple Color Emoji","Segoe UI Emoji","Segoe UI Symbol","Noto Color Emoji";color:#212529;text-align:left;box-sizing:border-box;margin-top:0;line-height:1.2;font-size:7rem;font-weight:700!important;margin-bottom:30px;text-align:center'><span style="margin-bottom:30px">ChargeFW<span style='color:#28a745'>2</span></span></div>

<p align="center">
  <a href="#compilation-requirements">Compilation requirements</a> •
  <a href="#installation">Installation</a> •
  <a href="#docker">Docker</a> •
  <a href="#usage">Usage</a> •
  <a href="#how-to-cite">How To Cite</a> •
</p>

<h4 align="center">Application for computing partial atomic charges using selected empirical methods.</h4>
<h4 align="center">ChargeFW2 is the computational core of <a href="https://acc2.ncbr.muni.cz" target="_blank">ACC2II</a>.</h4>

<div style='height: 50px'></div>

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
