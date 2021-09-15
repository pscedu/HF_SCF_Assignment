# Assignment #2 Hartree-Fock Self-Consistent Field Program

This assignment is meant to give you some experience using NumPy and SciPy
in a scientific code.

## Introduction
This assignment has been adapted from the ["Programming Tutorial in Chemistry
by Python" originally developed by Daniel Crawford](https://pycrawfordprogproj.readthedocs.io/en/latest/index.html)

The Hartree-Fock (HF) method is an approximation that allows one to solve the
molecular Schrodinger equation by assuming that the potential felt by
an electron is determined as a mean-field of all of the other electrons in
the system. Essentially this allows us to represent the wavefunction of a
molecule as a Slater determinant. We can then use the variational method to
show that solving equations is a minimization process that can be implemented
iteratively, which is termed the Self-Consistent Procedure (SCF)

To learn more about HF theory, it is recommended that you start with the notes
from [C. David Sherrill at the Georgia Institute of Technology](http://vergil.chemistry.gatech.edu/notes/hf-intro/hf-intro.html) or Chapter 3 of Szabo and Ostlund's [Modern Quantum Chemistry: Introduction to Advanced Electronic
Structure Theory](https://cmu.primo.exlibrisgroup.com/discovery/fulldisplay?docid=alma991002558099704436&context=L&vid=01CMU_INST:01CMU&lang=en&search_scope=MyInst_and_CI&adaptor=Local%20Search%20Engine&tab=Everything&query=any,contains,Szabo%20Ostlund&offset=0)

This project uses the [PySCF Module](https://pyscf.org/), see References at the bottom for citation.

## The Assignment

*Assumptions*: You have Python installed on your system

To complete the assignment:

1. Fork this repository into your own GitHub account.
2. Clone your forked repository
3. Go into your new directory and create a virtualenv with he command

>`python -m venv ./venv`

4. Activate your new virtual environment

>`source ./venv/bin/activate`

5. Install the necessary modules

>`pip install -r requirements.txt`

6. Now you are ready to start programming. There are two python files.

- main.py - this is the main controller program, it has the basic SCF
procedure outlined with all of the function calls that will be needed. You should not need to edit this file, but you should review it as it will be needed to run the final program.

- SCF.py - this is a file with stubs for all of the functions that you will need to implement. The goal will be to fill in each function. Notes are given below about each step of the procedure.

7. Submit your assignment by creating a pull request for your fork to the original repository and assign review to Shawn Brown.

**Note**: This process will allow you to see others assignments, while we encourage you to share thoughts and work together, but the purpose of this is not to get a grade, but to learn programming, so we are trusting you to do the assignment without cheating.

### The SCF Procedure



## References

1. [Programming Tutorial in Chemistry by Python](https://pycrawfordprogproj.readthedocs.io/en/latest/index.html), Daniel Crawford.

1. [An Introduction to Hartree-Fock Molecular Orbital Theory](http://vergil.chemistry.gatech.edu/notes/hf-intro/hf-intro.html), C. David Sherrill.

1. Szabo, A., Ostlund, N. S. (1996). [Modern Quantum Chemistry: Introduction to Advanced Electronic
Structure Theory](https://cmu.primo.exlibrisgroup.com/discovery/fulldisplay?docid=alma991002558099704436&context=L&vid=01CMU_INST:01CMU&lang=en&search_scope=MyInst_and_CI&adaptor=Local%20Search%20Engine&tab=Everything&query=any,contains,Szabo%20Ostlund&offset=0). Mineola: Dover Publications, Inc., available at CMU Library.

1. [Recent developments in the PySCF program package](https://doi.org/10.1063/5.0006074), Q. Sun, X. Zhang, S. Banerjee, P. Bao, M. Barbry, N. S. Blunt, N. A. Bogdanov, G. H. Booth, J. Chen, Z.-H. Cui, J. J. Eriksen, Y. Gao, S. Guo, J. Hermann, M. R. Hermes, K. Koh, P. Koval, S. Lehtola, Z. Li, J. Liu, N. Mardirossian, J. D. McClain, M. Motta, B. Mussard, H. Q. Pham, A. Pulkin, W. Purwanto, P. J. Robinson, E. Ronca, E. R. Sayfutyarova, M. Scheurer, H. F. Schurkus, J. E. T. Smith, C. Sun, S.-N. Sun, S. Upadhyay, L. K. Wagner, X. Wang, A. White, J. Daniel Whitfield, M. J. Williamson, S. Wouters, J. Yang, J. M. Yu, T. Zhu, T. C. Berkelbach, S. Sharma, A. Yu. Sokolov, and G. K.-L. Chan, J. Chem. Phys. 153, 024109 (2020)

1. [PySCF: the Python-based simulations of chemistry framework](https://wires.onlinelibrary.wiley.com/doi/10.1002/wcms.1340), Q. Sun, T. C. Berkelbach, N. S. Blunt, G. H. Booth, S. Guo, Z. Li, J. Liu, J. McClain, S. Sharma, S. Wouters, and G. K.-L. Chan, WIREs Comput. Mol. Sci. 8, e1340 (2018).

1. [Libcint: An efficient general integral library for Gaussian basis functions](https://doi.org/10.1002/jcc.23981), Q. Sun, J. Comp. Chem. 36, 1664 (2015).
