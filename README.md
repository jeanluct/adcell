# adcell

![demos][demos]

Matlab code for solving the advection-diffusion equation for a
two-dimensional incompressible autonomous cellular flow.

## Installation

- Clone the *adcell* project:
  ```
  git clone https://github.com/jeanluct/adcell.git
  ```
  or with SSH
  ```
  git clone git@github.com:jeanluct/adcell.git
  ```

- After starting Matlab, `cd` to the `extern` subfolder of the
  *adcell* project, and run `mex fft2udotgrad_helper.c` from within
  Matlab.  You might need to use the `-largeArrayDims` option if
  Matlab complains.  This will compile a helper function for filling
  the sparse matrix of the Fourier-space operator `u.grad` in the
  advection-diffusion equation.

- If you cannot compile the MEX file, you'll be limited to `N=31`
  Fourier modes in each spatial direction.

- Now `cd ..` to come back to the root folder of the *adcell* project.

You're now ready to run the program!

## Running

- Run `adcell_demo('demoX')`, where `X` is a number from 1 to 6.  The
  argument specifies which flow and set of parameter values are used.
  The banner above is the output of `'demo1'`, `'demo2'`, and
  `'demo3'` at the final time.

- You can add your own run parameters by adding a case in the `switch`
  statement at the top of `adcell_demo`: you can specify the diffusion
  coefficient (effectively the inverse Peclet number) and the number
  of Fourier coefficients in each dimension.  `N` should be made
  larger as `Diff` is decreased.  Another useful parameter is `ks`,
  which determines how many cells fill the domain.  For `ks=4`, there
  will seem to be 8 cells in each direction, since they come in pairs.

- See the `examples` folder for more sample code.

- You can see a demo of the program with commentary in this [lecture][1].

## License

*adcell* is released under the [MIT License][2].  It uses [*jlt
lib*][3] and the Fourier differentiation matrix function
[fourdif.m][4] written by S. C. Reddy and J. A. C. Weideman.

[1]: https://youtu.be/-uqkG1i3UHU
[2]: https://github.com/jeanluct/adcell/raw/master/LICENSE
[3]: https://github.com/jeanluct/jlt
[4]: http://appliedmaths.sun.ac.za/~weideman/research/differ.html
[demos]: https://github.com/jeanluct/adcell/raw/master/img/demos_small.png
[demo1]: https://github.com/jeanluct/adcell/raw/master/img/demo1_small.png
[demo2]: https://github.com/jeanluct/adcell/raw/master/img/demo2_small.png
[demo3]: https://github.com/jeanluct/adcell/raw/master/img/demo3_small.png
