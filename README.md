# 1D-ADE
A FEM solver for 1D steady advection-diffusion equation. The equation solved here is 

$\frac{Pe}{h}\frac{du}{dx} = \frac{d^2u}{dx^2}  ,x\in(0,L)$

where $Pe$ is cell Peclet number, $h$ is the length of elements.

In this code, the equation is defined by a class `Problem`, where 6 parameters `Pe`, `h`, `L`, `bc1`, `bc2`, `base_function` are required.

2 kinds of basis function is avalible here. `base_function=1` means linear basis function, and `base_function = 2` mean quadratic basis function.
