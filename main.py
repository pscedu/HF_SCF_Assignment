"""
This is the Hartree-Fock SCF Assignmentthat is for the class assignment for
Modern Programming in Data Analytics

It is intended to give students a basic exposure to the following concepts

1. Working with novel libraries
2. Converting math to python programming
3. Basic exposure to using NumPy and SciPy
4. Proper documentation of code and unit testing

The instructions for the assignement come in the Readme from the repository
at:
"""
from pyscf import gto #PySCF is a quantum chemistry python module


def main():
    """
    The main function call, meant to gather input and invoke the
    Hartree Fock Self Consistent Field Procedure

    Output:
        The total Hartee-Fock Energy of the molecule
    """

    ### Input Data
    mol = mol_h2o = gto.M(unit="Bohr",
                          atom = "O 0.000000000000  -0.143225816552   0.000000000000;"
                               + "H 1.638036840407   1.136548822547  -0.000000000000;"
                               + "H -1.638036840407   1.136548822547  -0.000000000000",
                          basis='STO-3g')
    mol_h2o.build()

    ### Convergence Criteria
    E_conv_threshold = 1.0E-10
    D_conv_threshold = 1.0E-8
    max_iterations = 1000

    ### get the integrals

    Suv = mol_h2o.intor('int1e_ovlp')
    Tuv = mol_h2o.intor('int1e_kin')
    Vuv = mol_h2o.intor('int1e_nuc')
    eri = mol_h2o.intor("int2e")


if __name__ == "__main__":
    main()
