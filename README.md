# PACS Project - Eikonal Equation Solver

### Authors: Francesco Fainello, Francesca Venturi

## Project Overview

This project focuses on developing efficient, implicit numerical methods for solving the **Eikonal equation**, a fundamental equation in modeling wave propagation and computing distance functions. The Eikonal equation has broad applications across fields such as: seismology, geospatial physics and computer vision.

The **Eikonal Equation** is given by:

```math
\|\nabla u(x)\| = F(x)
```

Where:
- \( u(x) \) is the distance function.
- \( F(x) \) represents the wavefront propagation speed.

To solve this, we derive the following **variational formulation**:

```math
E(u(x)) = \int_{\Omega} \left( \|\nabla u(x)\| - 1 \right)^2 \, dx
```

This energy functional is minimized to compute the solution for the Eikonal equation.

## Problem Statement

This project addresses the numerical solution of this equation by employing **energy minimization principles** and applying **iterative methods** such as:
- **Method of Laplacian Iterations**
- **Relaxation of Laplacian Iterations**
- **Alternating Direction Method of Lagrange Multipliers (ADMM)**

These methods allow for efficient, scalable solutions applicable to a variety of complex scenarios.

---
## Usage Instructions

To set up and run the solver, follow these steps:

1. **Create a Project Directory:**
    -  First, create a directory where all the necessary files will be created and downloaded:
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

4. **Compilation:**
   - Compile the code by running:
     ```bash
     make
     ```
   - This will generate the executable files needed to run simulations.

5. **Execution**

   a. **Prepare the Environment**:  
      Make sure that all the necessary libraries and dependencies are installed. The simulation is designed to be run in a parallel computing environment using **MPI**.

      Ubuntu/Linux:
      ```bash
      sudo apt-get update
      sudo apt-get install openmpi-bin libopenmpi-dev
      ```

      MacOS:
      ```bash
      brew install open-mpi
      ```

   b. **Run the Simulation**:  
      To execute the main simulation, use the following command:

      ```bash
      mpirun -np <num_processes> ./main_executable
      ```

   #### Input Configuration

   During execution, the user will be prompted to provide input options for the simulation setup. The steps are outlined below:

   c. **Mesh Input**:  
      You will be asked to input the name of the `.vtk` file containing the mesh data. If no input is provided, the default file `mesh_10.vtk` will be used. This default mesh represents a **unitary cube** divided into 10 uniform intervals along each edge.

      Example prompt:
      ```
      Enter the name of the .vtk file containing the mesh (default: mesh_10.vtk):
      ```

   d. **Boundary Conditions**:  
      You will be asked to choose the boundary conditions to impose on the problem. The available options are:
      - **1**: Dirichlet Null BCs on the bottom face (z=0)  
      - **2**: Dirichlet Null BCs on the central point of the bottom face (z=0, x=0.5, y=0.5)  
      - **3**: Dirichlet Null BCs on one vertex (x=0, y=0, z=0)  
      - **4**: Dirichlet Null BCs on the whole boundary  

      If no input is provided, **Option 1** is used by default.

   e. **Method Selection**:  
      Next, you will be asked to select the method for solving the Eikonal equation. The possible options are:
      - **1**: Standard Eikonal  (Method of Laplacian Iterations)
      - **2**: Penalty Eikonal  (Relaxation of Laplacian Iterations)
      - **3**: Lagrangian Eikonal  (Alternating Direction Method of Lagrange Multipliers)

      If no input is provided, **Option 1 (Standard Eikonal)** is selected by default.

6. **Generating Plots:**
   - To generate the plots you will find in the report, run:
     ```bash
     cd plots
     ./plots.sh
     ```    

7. **Generating Documentation:**
   - To generate documentation, run:
     ```bash
     make docs
     ```

8. **Cleaning the Build:**
   - To clean up the build files, use:
     ```bash
     make clean
     ```

These steps will guide you through setting up, compiling, and running the solver, along with generating documentation and cleaning up your build environment.
