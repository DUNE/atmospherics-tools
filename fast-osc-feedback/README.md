# FastOdcFeedback Project

## Overview

FastOdcFeedback is a project designed to process and analyze oscillogram data using various scientific libraries and tools. The project leverages ROOT, Uproot, and other scientific libraries to handle data processing and visualization.

## Project Structure

```
CMakeLists.txt
fast-feedback.ipynb
fastfeedback.py
README.md
src/
    CMakeLists.txt
    LinkDef.h
    Oscillogram.cxx
    Oscillogram.h
```

### Key Files and Directories

- **[`CMakeLists.txt`](CMakeLists.txt)**: The main CMake configuration file for the project.
- **[`fast-feedback.ipynb`](fast-feedback.ipynb)**: A Jupyter notebook for interactive data analysis.
- **[`fastfeedback.py`](fastfeedback.py)**: The main Python file with tools for data processing and analysis.
- **[`README.md`](README.md)**: This documentation file.
- **[`src/`](src/)**: Contains source files for the Oscillogram library.
  - **[`CMakeLists.txt`](src/CMakeLists.txt)**: CMake configuration for the Oscillogram library.
  - **[`LinkDef.h`](src/LinkDef.h)**: ROOT dictionary generation file.
  - **[`Oscillogram.cxx`](src/Oscillogram.cxx)**: Implementation of the Oscillogram library.
  - **[`Oscillogram.h`](src/Oscillogram.h)**: Header file for the Oscillogram library.

## Dependencies

The project relies on several external libraries and tools:

- **ROOT**: A framework for data processing.
- **Uproot**: A library for reading ROOT files in Python.
- **Awkward Array**: A library for manipulating complex data structures.
- **Pandas**: A data analysis and manipulation library.
- **Matplotlib**: A plotting library.
- **Sparse**: A library for sparse matrices.
- **Polars**: A fast DataFrame library.

## Installation

To build and install the project, follow these steps:

1. **Clone the repository**:
    ```sh
    git clone <repository-url>
    cd fast-osc-feedback
    ```

2. **Build the project using CMake**:
    ```sh
    mkdir build
    cd build
    cmake ..
    make
    ```

## Usage

### Tools

The main classes for data processing are in [`fastfeedback.py`](fastfeedback.py). You can run it using Python:

```sh
python fastfeedback.py
```

### Jupyter Notebook

For interactive data analysis, you can use the Jupyter notebook [`fast-feedback.ipynb`](fast-feedback.ipynb). Launch Jupyter and open the notebook:

```sh
jupyter notebook fast-feedback.ipynb
```

## Configuration

The project uses CMake for configuration. The main configuration file is [`CMakeLists.txt`](CMakeLists.txt). It includes external projects and sets up the build environment.

### External Projects

The project includes an external project, OscProb, which is configured in the main [`CMakeLists.txt`](CMakeLists.txt):

```txt
ExternalProject_Add(OscProb
    DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}
    GIT_REPOSITORY https://github.com/joaoabcoelho/OscProb
    SOURCE_DIR ${OSCPROB_SOURCE_DIR}
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${MAKE}
    INSTALL_COMMAND "${oscprob_lib_install_commands}"
)
```

### Oscillogram Library

The Oscillogram library is defined in [`src/CMakeLists.txt`](src/CMakeLists.txt):

```txt
add_library(Oscillogram SHARED Oscillogram.cxx)
target_include_directories(Oscillogram PUBLIC ${OSCPROB_INCLUDE_DIR})
ROOT_GENERATE_DICTIONARY(Oscillogram_dict Oscillogram.h MODULE Oscillogram LINKDEF LinkDef.h)
target_link_libraries(Oscillogram PRIVATE ROOT::RIO)
```

## Contributing

Contributions are welcome! Please fork the repository and submit pull requests.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Contact

For any questions or issues, please open an issue on the GitHub repository or contact the project maintainers.

---