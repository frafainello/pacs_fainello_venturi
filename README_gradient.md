# Tools for computing gradients and basic finite different matrices #
In `gradient.hpp` you find several tools for computing gradients and basic linear finite element matrices. 

- Elemental Mass Matrix Either in the form
  $$
    M = f \int_{\K} \phi_i(x) \phi_j(x) \, dx
  $$
where $f$ is a factor representing the constant value of a function over the finite element $K$,
or in the form
$$
    M = \int_{\K} \phi_i(x) \phi_j(x) s(x) \, dx
$$
where $s(x)$ is a linear scalar function represented by its values at the element nodes: $s(x) = \sum_{k=1}^n s_k \phi_k(x)$.
- Elemental Stiffness Matrix In the form
$$
    K = f \int_{\K} \nabla \phi_i(x) \nabla \phi_j(x) \, dx
$$
where $f$ is a factor representing the constant value of a function over the finite element $K$.

There are also other utilities to compute the gradient of a function over a linear finite element.

All code is doxygenated and you can produce the documentation with doxygen
## Note ##
To avoid useless reallocation of matrices the functions for computing the various element matrices may be transformed in methods of a class that stores
the matices as member variables. This is not done here to keep the code as simple as possible, but it maybe an useful improvement for the future.