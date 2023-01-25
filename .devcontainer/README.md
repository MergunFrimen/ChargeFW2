# Dev Container

## Requirements & Installation
[See this guide](https://code.visualstudio.com/docs/devcontainers/containers#_system-requirements)

## How To Use

VSCode should automatically prompt you to reopen the project in a container.
If not, you can manually open the project in a container by running the command
`>Dev Containers: Open Folder in Container...`
from the Command Palette. VSCode might then further prompt you to rebuild the container.

After the container is done building, run the following commands to build the project:

```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=.
make
make install
```

> The default `CMAKE_INSTALL_PREFIX=/usr/lib` won't work because the devcontainer is running unprivileged and won't be able to install files into that directory.