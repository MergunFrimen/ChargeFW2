# Dev Container

## Requirements
https://code.visualstudio.com/docs/devcontainers/containers#_system-requirements

## Installation
https://code.visualstudio.com/docs/devcontainers/containers#_installation

## How to use

VSCode should automatically prompt you to reopen the project in a container.
If not, you can manually open the project in a container by running the command `>Dev Containers: Open Folder in Container...` from the command palette.

After the container is built and VSCode finishes loading the devcontainer environment, don't forget to run:

```bash
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=.
```

to update *CMakeCache.txt* (the default `CMAKE_INSTALL_PREFIX=/usr/lib` won't work because the devcontainer is running unpriviledged).