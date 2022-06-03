CLUSTER DEComposition suite
---------------------------

This repository contains software for re-expressing configuration
interaction (CI) wave functions as coupled cluster (CC)
expansions. The software is composed of two parts

1) a program that generates the necessary recursion relations and
stores the instructions on disk

2) a program that reads in a CI wave function and translates it into a
CC expansion.

If you use this program, please cite the following reference:

S. Lehtola, N. M. Tubman, M. Head-Gordon, and K. B. Whaley, "Cluster
decomposition of full configuration interaction wave functions: a tool
for chemical interpretation of systems with strong correlation",
J. Chem. Phys. 147, 154105 (2017). arXiv:1707.04376.

For further information on the program, read the above reference.

To compile the code, run make.

In contrast to the article, the program comes in a number of variants
that are meant for different input formats for the CI wave function.
# clusterdec_bit.x reads wave functions given as bit strings: coefficient alphastring betastring, where the strings are 0s and 1s indicating whether there is an electron in the orbital
# clusterdec_ubit.x reads wave functions given as unified bit stings: coefficient string, where the orbital occupations are 0, u, d, or 2, meaning no occupation, occupation by a spin-up electron, a spin-down electron, or two electrons, respectively
# clusterdec_int.x (built by make all_gmp, requires the GNU MultiPrecision Arithmetic Library to build) reads wave functions given as integers: coefficient aint bint, where aint and bint are the integer representations of the bitstring

In addition, the program also includes tools to convert wave functions from the abovementioned formats to each other: the programs convert_X_Y.x convert from the X format (bit, ubit, int) to the Y format (bit, ubit, int).


2022-06-03 Susi Lehtola
