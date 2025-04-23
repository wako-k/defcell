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
    space_group = symmetry[6]
#    print("space_group", space_group)

#    print("symmetry:", symmetry)
#    print("component_strain:", component_strain)
#    print("amount_strain:", amount_strain)

    # calculate the lattice vectors
    lattice_vectors = calc_lattice_vector(symmetry)
#    print("lattice_vectors", lattice_vectors)

    # calculate the strain tensor
    strain_tensor = calc_strain_tensor(component_strain, amount_strain)
#    print("strain_tensor", strain_tensor)

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

    a1 = lattice_vectors[0]
    a2 = lattice_vectors[1]
    a3 = lattice_vectors[2]
#    print("a1", a1)
#    print("a2", a2)
#    print("a3", a3)

    a = np.linalg.norm(a1)
    b = np.linalg.norm(a2)
    c = np.linalg.norm(a3)
#    print("a", a)
#    print("b", b)
#    print("c", c)

    cos_a = np.dot(a2, a3) / (b * c)
    cos_b = np.dot(a3, a1) / (c * a)
    cos_c = np.dot(a1, a2) / (a * b)
#    print("cos_a", cos_a)
#    print("cos_b", cos_b)
#    print("cos_c", cos_c)

    alpha = np.degrees(np.arccos(cos_a))
    beta  = np.degrees(np.arccos(cos_b))
    gamma = np.degrees(np.arccos(cos_c))

#    print("alpha", alpha)
#    print("beta", beta)
#    print("gamma", gamma)

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

#    print("strain_vector", strain_vector)

    # convert to strain tensor
    strain_tensor = np.array([
        [strain_vector[0],     strain_vector[5] / 2, strain_vector[4] / 2],
        [strain_vector[5] / 2, strain_vector[1],     strain_vector[3] / 2],
        [strain_vector[4] / 2, strain_vector[3] / 2, strain_vector[2]]
    ])
#    print("strain_tensor", strain_tensor)

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
#    print("cell lengths:", a, b, c)

    # the unit cell angles
    alpha = symmetry[3]
    beta  = symmetry[4]
    gamma = symmetry[5]
#    print("cell angles:", alpha, beta, gamma)

    # the space group name
    space_group = symmetry[6]
#    print("space_group:", space_group)

    #
    # transform to the lattice vector
    #

    # convert to radian
    alpha_rad, beta_rad, gamma_rad = np.radians([alpha, beta, gamma])
#    print("alpha_rad, beta_rad, gamma_rad:", alpha_rad, beta_rad, gamma_rad)

    # the vector of a1
    a1 = np.array([a, 0, 0])
#    print("a1:", a1)

    # the vector of a2
    a2 = np.array([b * np.cos(gamma_rad), b * np.sin(gamma_rad), 0])
#    print("a2:", a2)

    # the vector of a3
    cx = c * np.cos(beta_rad)
    cy = c * (np.cos(alpha_rad) - np.cos(beta_rad) * np.cos(gamma_rad)) / np.sin(gamma_rad)
    cz = np.sqrt(c**2 - cx**2 - cy**2)
    
    a3 = np.array([cx, cy, cz])
#    print("a3:", a3)

    lattice_vectors = np.array([a1, a2, a3])
#    print("lattice_vectors:", lattice_vectors)

    return lattice_vectors


####################################################################################################
def move_molecules(component_strain, amount_strain, mol_list):
    """move the molecules"""

#    print("mol_list", mol_list)

    # calculate the strain tensor
    strain_tensor = calc_strain_tensor(component_strain, amount_strain)
#    print("strain_tensor", strain_tensor)

    # calcuate the deformation gradient tensor (F = I + e)
    deformation_gradient = np.eye(3) + strain_tensor
#    print("deformation_gradient", deformation_gradient)

    # loop for the all molecules
    for mol in mol_list:
#        print("mol:", mol)

        # get the position vector for the center of mass
        r = np.array(cmd.centerofmass(mol))
#        print("r", r)

        # calculate the deformed potition
        new_r = np.dot(deformation_gradient, r)
#        print("new_r", new_r)

        # calculate the sift vector
        shift_vector = new_r - r
#        print("shift_vector", shift_vector)

        # calculate a shift vector for single molecule
#        shift_vector = calc_shift_vector(component_strain, amount_strain, mol)
        
        # move the molecule
        cmd.translate(shift_vector.tolist(), mol)

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

