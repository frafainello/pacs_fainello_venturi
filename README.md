# PACS Project - Eikonal Equation

### Author: Francesco Fainello, Francesca Venturi

## Project Overview

## Problem Statement

## Usage Instructions

To set up and run the solver, follow these steps:

1. **Create a Project Directory:**
    -  First, create a directory where all the necessary files will be created and dowloaded:
        ```bash
        mkdir project
        ```

2. **Clone the Repositories:**
   - Then, clone the necessary repositories with the following command:
     ```bash
     cd project
     git clone --recursive git@github.com:pacs-course/pacs-examples.git
     git clone --recursive git@github.com:frafainello/pacs_fainello_venturi.git
     ```
   - This will load the `pacs-examples` and `pacs_fainello_venturi` folders into your project directory.

3. **Modify Makefile:**
   - After cloning, modify the `PACS_ROOT` variable in the `Makefile` to point to the `pacs-examples` folder. This step ensures that the project is correctly configured with the necessary paths. Then navigate into the `pacs_fainello_venturi` folder:
    ```bash
     cd pacs_fainello_venturi
     ```

4. **Load Necessary Modules:**
   - Before proceeding with the installation, load the required modules by running:
     ```bash
     chmod +x scripts/load_modules.sh
     scripts/load_modules.sh
     ```

5. **Installation:**
   - Run:
     ```bash
     make install
     ```
   - This command installs the necessary dependencies and prepares the environment.

6. **Compilation:**
   - Compile the code by running:
     ```bash
     make
     ```
   - This will generate the executable files needed to run simulations.

<!-- 7. **Set Library Path:**
   - Set the library path to include the `lib` directory:
     ```bash
     LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(pwd)/lib
     ``` -->

7. **Execution:**
   - Execute the main simulation by running:
     ```bash
     ./main_executable
     ```

8. **Generating Documentation:**
   - To generate documentation, run:
     ```bash
     make docs
     ```

9. **Cleaning the Build:**
   - To clean up the build files, use:
     ```bash
     make clean
     ```

These steps will guide you through setting up, compiling, and running the solver, along with generating documentation and cleaning up your build environment.


## Folder Structure
