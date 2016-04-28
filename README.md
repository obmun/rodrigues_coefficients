Evaluation of different numerical techniques for the calculation of Ritto-Correa's coefficients
===============================================================================================

Ritto-Correa's a<sub>i</sub>, b<sub>i</sub> and c<sub>i</sub> coefficients are a set of values
introduced by Ritto-Corrêa and Camotim in the paper ["On the differentiation of the Rodrigues
formula and its significance for the vector-like parameterization of Reissner–Simo beam
theory"][1]. These coefficients are helpful in  the calculation of the [directional
derivatives][2] of the _rotation vector parametrization_ of the _orthogonal tensor_ ([orthogonal
matrix][3]) representing a finite rotation. They appear in the expression that relates the rotation
vector to the rotation tensor. In Ritto-Corrêa's paper the [Rodrigue's formula is used][6]. This
derivatives and the tensor parametrization are useful in continuum mechanics beam theory.

The expressions for a<sub>i</sub>, b<sub>i</sub> and c<sub>i</sub> can be found in sections 2.1 and
2.2 of Ritto-Correa's paper. b<sub>i</sub> and c<sub>i</sub> are related to the 1<sup>st</sup> and
2<sup>nd</sup> order derivatives of a<sub>i</sub>.

In this code we implement 3 different ways of calculating these coefficients:
  - By direct application of the corresponding _symbolic expressions_ based on elemental functions
  and its direct implementation using the corresponding C++ library math functions.
  - Using the series expansion of the coefficients, as found in [section 2.3 of the paper][1].
  - Using [Fike's and Alonso's hyper-dual numbers][6].

This code is written in C++11 and only makes uses of std library functions and the included
hyper-dual class as implemented by Fike and slightly modified by me. Build system is CMake.

Hyper-dual numbers
------------------

The hyper-dual numbers are an improvement over the concept of the **complex-step derivative
approximation**. For a brief but very clear introduction to complex-step derivative
approximations, you can take a look at section 2.1 of [Martins, Sturdza and Alonso's
paper][5]. By means of complex numbers, we can find approximations to the
derivative of a function without having to use finite differences and, consequently, **avoiding the
 substractive cancellation error** typical of finite differences approximations.

Fike and Alonso introduced the hyper-dual numbers in order to facilitate a method for
**second-derivative calculation** that was comparable to the complex-step approximation for the
first derivative in terms of accuracy, efficiency and ease of implementation. For more details
please check [Fike's paper][7].

Why this code
-------------

The idea of this work was not to obtain an approximate value of the coefficients, as we know their
analytic expressions, but to analyse its evaluation accuracy using different methods. The key issue
with the a<sub>i</sub>, b<sub>i</sub> and c<sub>i</sub> expressions is that they are prone to
cancellation errors.

Building and execution
----------------------

Nothing special. You just need a C++11 compliant compiler and CMake. Then:

    mkdir build
    cd build
    cmake ../
    make

The output of the program, for the moment, is a _big_ table of numbers with the evaluation of some
of the coefficients at multiple points around 0. Redirect it to a file or pipe it to `less` in the
following fashion:

    ./derivatives | less -S

Some additional notes
---------------------

### Rotations as tensors ###

Effectively, a rotation matrix is not only a matrix but a tensor. In continuum mechanics the
usage of tensors is common, so technical papers frequently speak about the rotation tensors. For
more information, check for example [section 1.3 on this Brown Uni web][4].

[1]: http://onlinelibrary.wiley.com/doi/10.1002/nme.532/full
[2]: https://en.wikipedia.org/wiki/Directional_derivative
[3]: https://en.wikipedia.org/wiki/Rotation_matrix
[4]: http://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Tensors/Tensors.htm
[5]: http://www.math.u-psud.fr/~maury/paps/NUM_CompDiff.pdf
[6]: https://en.wikipedia.org/wiki/Rodrigues'_rotation_formula
[7]: http://adl.stanford.edu/hyperdual/Fike_AIAA-2011-886.pdf