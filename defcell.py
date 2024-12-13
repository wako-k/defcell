#!/usr/bin/env python3
####################################################################################################
#
# defcell.py: deform the unit cell using PyMOL API
#
# Copyright (C) 2024 Kei Wako
#
####################################################################################################
import sys
import argparse
 
# Tell PyMOL to launch quiet (-q), without script (-k) and without internal GUI (-i)
import __main__
__main__.pymol_argv = [ 'pymol', '-qki' ]
 
# Call the function below before using any PyMOL modules.
#import pymol
#pymol.finish_launching()  # not supported on macOS

from pymol import cmd

####################################################################################################
def main():

    # parse the command line options
    parser = argparse.ArgumentParser()
    parser.add_argument('pse_file', metavar='FILENAME',
                        help='file name of the PyMOL pse file')
    parser.add_argument('-s', '--strain', metavar='STRAIN', type=float,
                        help='specify a strain (%%)')
    parser.add_argument('-o', '--output', metavar='FILENAME', default='out.pdb',
                        help='output filename for pdb file (default: %(default)s)')
    args = parser.parse_args()

    # file name of the PyMOL pse file
    pse_file = args.pse_file

    # strain (%)
    strain_percent = args.strain

    # ratio of the strain
    strain = strain_percent/100.0

    # file name of the output pdb file
    out_file = args.output

    # open the PyMOL pse file
    cmd.load(pse_file)

    # get the names of the molecules
    mol_list = cmd.get_object_list()
#    print(mol_list)

    # deform the unit cell
    deform_cell(strain, mol_list)

    # move the molecules
    move_molecules(strain, mol_list)

    # save the file
    cmd.save(out_file)

####################################################################################################
def deform_cell(strain, mol_list):
    """deform the unit cell"""

    # get the first molecule
    first_mol = mol_list[0]
#    print("first_mol:", first_mol)

    # get the unit cell lengths
    symmetry = cmd.get_symmetry(first_mol)
#    print("symmetry:", symmetry)

    # the unit cell length
    cell_a = symmetry[0]
    cell_b = symmetry[1]
    cell_c = symmetry[2]
#    print("cell lengths:", cell_a, cell_b, cell_c)

    # the unit cell angles
    cell_alpha = symmetry[3]
    cell_beta  = symmetry[4]
    cell_gamma = symmetry[5]
#    print("cell angles:", cell_alpha, cell_beta, cell_gamma)

    # the space group name
    space_group = symmetry[6]
#    print("space_group:", space_group)

    # deform the cell length
    new_a = cell_a + (cell_a * strain)
#    print("new_a:", new_a)

    # change the cell size
    cmd.set_symmetry(first_mol, new_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma, space_group)

    # conform
#    new_symmetry = cmd.get_symmetry(first_mol)
#    print("new_symmetry:", new_symmetry)

    return

####################################################################################################
def move_molecules(strain, mol_list):
    """move the molecules"""

    # loop for the all molecules
    for mol in mol_list:
#        print("mol:", mol)

        # get the center of mass
        com = cmd.centerofmass(mol)
#        print("com:", com)

        # the x coordinate
        x = com[0]
#        print("x:", x)

        # calculate the shift vector
        shift = strain * x
#        print("shift:", shift)
        
        # move the molecule
        cmd.translate([shift, 0, 0], mol)

        # conform
#        new_com = cmd.centerofmass(mol)
#        print("new_com:", new_com)

    return

####################################################################################################
if __name__ == "__main__":
    # execute only if run as a script
    main()

