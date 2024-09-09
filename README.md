# PACS Project - Eikonal Equation

### Author: Francesco Fainello, Francesca Venturi

### Project Overview

This project is centred on developing efficient, implicit numerical methods for solving the **Eikonal equation**, which plays a key role in modelling wave propagation and computing distance functions. Applications of the Eikonal equation span across diverse fields such as seismology, geospatial physics, and computer vision. 

The **Eikonal Equation** reads:
   \[
   \|\nabla u(x)\| = F(x)
   \]
where \( u(x) \) is the distance function and \( F(x) \) represents the wavefront propagation speed. This leads to the following variational formulation:
\[
   E(u(x)) = \int_{\Omega} \left( \|\nabla u(x)\| - 1 \right)^2 \, dx
   \]
   This energy functional is minimized to solve the Eikonal equation.


The focus of this project is on implicit techniques based on **variational formulations**, which are more robust and flexible when dealing with complex boundary conditions and anisotropic media compared to traditional explicit methods like the fast marching method.

### Problem Statement

Accurate solutions of the Eikonal equation are crucial for simulations in complex domains, such as subsurface geological structures. However, explicit methods often fail to handle the nonlinearities and discontinuities in real-world models. This project addresses these limitations by applying **energy minimization principles** and employing **iterative methods** such as the **Alternating Direction Method of Multipliers (ADMM)** and **relaxation of Laplacian iterations**.

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
