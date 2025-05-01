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
import numpy as np
 
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

    # dummy data
#    symmetry = [1.0, 1.0, 1.0, 90, 90, 90, 'foo']
#    print("symmetry", symmetry)

    space_group = symmetry[6]

    # calculate the lattice vectors
    lattice_vectors = calc_lattice_vector(symmetry)
#    print("lattice_vectors", lattice_vectors)

    # calculate the strain tensor
    strain_tensor = calc_strain_tensor(component_strain, amount_strain)

    # calcuate the deformation gradient tensor (F = I + e)
    deformation_gradient = np.eye(3) + strain_tensor
#    print("deformation_gradient", deformation_gradient)

    # applay deformation to the lattice vectors
    deformed_lattice_vectors = np.dot(deformation_gradient, lattice_vectors)
#    print("deformed_lattice_vectors", deformed_lattice_vectors)

    # calculate the lattice constant
    new_symmetry = calc_lattice_constant(deformed_lattice_vectors, space_group)
#    print("new_symmetry", new_symmetry)

    return new_symmetry


####################################################################################################
def calc_lattice_constant(lattice_vectors, space_group):
    """calculate the lattice constant"""

    a1 = lattice_vectors[:, 0]
    a2 = lattice_vectors[:, 1]
    a3 = lattice_vectors[:, 2]

    a = np.linalg.norm(a1)
    b = np.linalg.norm(a2)
    c = np.linalg.norm(a3)

    cos_a = np.dot(a2, a3) / (b * c)
    cos_b = np.dot(a3, a1) / (c * a)
    cos_c = np.dot(a1, a2) / (a * b)

    alpha = np.degrees(np.arccos(cos_a))
    beta  = np.degrees(np.arccos(cos_b))
    gamma = np.degrees(np.arccos(cos_c))

    symmetry = [a, b, c, alpha, beta, gamma, space_group]

    return symmetry


####################################################################################################
def calc_strain_tensor(component_strain, amount_strain):
    """calculate the vector of the strain"""

    # create the strain vector
    # (exx, eyy, ezz, 2eyz, 2ezx, 2exy)

    match component_strain:
        case 'exx':
            strain_vector = np.array([amount_strain, 0.0, 0.0, 0.0, 0.0, 0.0])
        case 'eyy':
            strain_vector = np.array([0.0, amount_strain, 0.0, 0.0, 0.0, 0.0])
        case 'ezz':
            strain_vector = np.array([0.0, 0.0, amount_strain, 0.0, 0.0, 0.0])
        case 'eyz':
            strain_vector = np.array([0.0, 0.0, 0.0, amount_strain, 0.0, 0.0])
        case 'ezx':
            strain_vector = np.array([0.0, 0.0, 0.0, 0.0, amount_strain, 0.0])
        case 'exy':
            strain_vector = np.array([0.0, 0.0, 0.0, 0.0, 0.0, amount_strain])
        case _:
            print("Error. The component of strain is ", component_strain)
            sys.exit()

    # convert to strain tensor
    strain_tensor = np.array([
        [strain_vector[0],     strain_vector[5] / 2, strain_vector[4] / 2],
        [strain_vector[5] / 2, strain_vector[1],     strain_vector[3] / 2],
        [strain_vector[4] / 2, strain_vector[3] / 2, strain_vector[2]]
    ])

    return strain_tensor


####################################################################################################
def calc_lattice_vector(symmetry):
    """calclate lattice vectors from lattice constants"""

    #
    # get the lattice constant
    #

    # the unit cell length
    a = symmetry[0]
    b = symmetry[1]
    c = symmetry[2]

    # the unit cell angles
    alpha = symmetry[3]
    beta  = symmetry[4]
    gamma = symmetry[5]

    # the space group name
    space_group = symmetry[6]

    #
    # transform to the lattice vector
    #

    # convert to radian
    alpha_rad, beta_rad, gamma_rad = np.radians([alpha, beta, gamma])

    # the vector of a1
    a1 = np.array([a, 0, 0])

    # the vector of a2
    a2 = np.array([b * np.cos(gamma_rad), b * np.sin(gamma_rad), 0])

    # the vector of a3
    cx = c * np.cos(beta_rad)
    cy = c * (np.cos(alpha_rad) - np.cos(beta_rad) * np.cos(gamma_rad)) / np.sin(gamma_rad)
    cz = np.sqrt(c**2 - cx**2 - cy**2)
    
    a3 = np.array([cx, cy, cz])

    # create the matrix as a column vector
    lattice_vectors = np.column_stack((a1, a2, a3))

    return lattice_vectors


####################################################################################################
def move_molecules(component_strain, amount_strain, mol_list):
    """move the molecules"""

    # calculate the strain tensor
    strain_tensor = calc_strain_tensor(component_strain, amount_strain)

    # calcuate the deformation gradient tensor (F = I + e)
    deformation_gradient = np.eye(3) + strain_tensor

    # loop for the all molecules
    for mol in mol_list:

        # get the position vector for the center of mass
        r = np.array(cmd.centerofmass(mol))

        # calculate the deformed potition
        new_r = np.dot(deformation_gradient, r)

        # calculate the sift vector
        shift_vector = new_r - r
        
        # move the molecule
        cmd.translate(shift_vector.tolist(), mol)

        # conform
#        new_com = cmd.centerofmass(mol)
#        print("new_com:", new_com)

    return


####################################################################################################
if __name__ == "__main__":
    # execute only if run as a script
    main()

