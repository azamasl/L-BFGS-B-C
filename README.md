# L-BFGS-B-C

The C version of the well-known [L-BFGS-B code](http://users.eecs.northwestern.edu/~nocedal/lbfgsb.html), version 3.0.

TO BE CONTINUED

# Installation

To use in C, go to the `src/` subdirectory and type `make`. This is unnecessary if you just want to use the Matlab wrapper. For an example of how to call the library from C, see the `driver1.c` file. The included Makefile includes a test of the installation using the problem defined in `driver1.c`.

# License

L-BFGS-B is released under the BSD 3-clause license, and I am releasing this software under the same license. See LICENSE for details

The L-BFGS-B website requests that you cite them. From their website:
"Condition for Use: This software is freely available, but we expect that all publications describing  work using this software , or all commercial products using it, quote at least one of the references given below. This software is released under the "New BSD License" (aka "Modified BSD License" or "3-clause license"). "

It would be nice to cite this website as well since it took a significant amount of work...

# Authors
This is a fork of
C version and Matlab wrapper are written by [Stephen Becker](http://amath.colorado.edu/faculty/becker/), stephen.becker@colorado.edu

The L-BFGS-B algorithm was written in the 1990s (mainly 1994, some revisions 1996) by Ciyou Zhu (in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal)

L-BFGS-B Version 3.0 is an algorithmic update from 2011, with coding changes by J. L. Morales.

