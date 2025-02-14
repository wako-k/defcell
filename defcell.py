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
    parser.add_argument('-s', '--strain', metavar='COMPONENT',
                        help='specify a strain component (e.g. exx, eyy, ezz, eyz, ezx, exy)')
    parser.add_argument('-a', '--amount', metavar='STRAIN', type=float,
                        help='specify an amount of strain (%%)')
    parser.add_argument('-o', '--output', metavar='FILENAME', default='out.pdb',
                        help='output filename for pdb file (default: %(default)s)')
    args = parser.parse_args()

    # file name of the PyMOL pse file
    pse_file = args.pse_file

    # strain component
    component_strain = args.strain

    # strain (%)
    amount_strain_percent = args.amount

    # ratio of the strain
    amount_strain = amount_strain_percent/100.0

    # file name of the output pdb file
    out_file = args.output

    # open the PyMOL pse file
    cmd.load(pse_file)

    # get the names of the molecules
    mol_list = cmd.get_object_list()

    # deform the unit cell
    deform_cell(component_strain, amount_strain, mol_list)

    # move the molecules
    move_molecules(component_strain, amount_strain, mol_list)

    # save the file
    cmd.save(out_file)


####################################################################################################
def deform_cell(component_strain, amount_strain, mol_list):
    """deform the unit cell"""

    # get the first molecule
    first_mol = mol_list[0]

    # get the unit cell lengths
    symmetry = cmd.get_symmetry(first_mol)

    # calculate the unit cell parameters
    new_symmetry = calc_cell_parameter(component_strain, amount_strain, symmetry)

    a = new_symmetry[0]
    b = new_symmetry[1]
    c = new_symmetry[2]
    alpha = new_symmetry[3]
    beta  = new_symmetry[4]
    gamma = new_symmetry[5]
    space_group = new_symmetry[6]

    # change the cell size
    cmd.set_symmetry(first_mol, a, b, c, alpha, beta, gamma, space_group)

    return


####################################################################################################
def calc_cell_parameter(component_strain, amount_strain, symmetry):
    """calculate the unit cell parameters"""

    new_symmetry = []

#    print("symmetry:", symmetry)
#    print("component_strain:", component_strain)
#    print("amount_strain:", amount_strain)

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

    match component_strain:
        case 'exx':
#            print("exx")
            # deform the cell length
            new_a = cell_a + (cell_a * amount_strain)
#            print("new_a:", new_a)
            # no changes
            new_b = cell_b
            new_c = cell_c
            new_alpha = cell_alpha
            new_beta  = cell_beta
            new_gamma = cell_gamma


    new_symmetry.append(new_a)
    new_symmetry.append(cell_b)
    new_symmetry.append(cell_c)
    new_symmetry.append(cell_alpha)
    new_symmetry.append(cell_beta)
    new_symmetry.append(cell_gamma)
    new_symmetry.append(space_group)
#    print("new_symmetry:", new_symmetry)

    return new_symmetry


####################################################################################################
def move_molecules(component_strain, amount_strain, mol_list):
    """move the molecules"""

    # loop for the all molecules
    for mol in mol_list:
#        print("mol:", mol)

        # calculate a shift vector for single molecule
        shift_vector = calc_shift_vector(component_strain, amount_strain, mol)
        
        # move the molecule
        cmd.translate(shift_vector, mol)

        # conform
#        new_com = cmd.centerofmass(mol)
#        print("new_com:", new_com)

    return


####################################################################################################
def calc_shift_vector(component_strain, amount_strain, mol):
    """calculate a shift vector for single molecule"""

    # get the center of mass
    com = cmd.centerofmass(mol)
#    print("com:", com)

    # the x coordinate
    x = com[0]
#    print("x:", x)

    # calculate the shift vector
    shift = amount_strain * x
#    print("shift:", shift)

    shift_vector = [shift, 0, 0]

    return shift_vector


####################################################################################################
if __name__ == "__main__":
    # execute only if run as a script
    main()

