import numpy as np
from math import sqrt
from math import cos
from math import sin
from math import pi
from random import random
dmin=2.0
dmax=12.0
step_size=3.0
layer=6.0
big_start_value=1000

# define the function for the rotation matrix
def rotation(theta):
    tx,ty,tz = theta
    Rx = np.array([[1,0,0], [0, cos(tx), -sin(tx)], [0, sin(tx), cos(tx)]])
    Ry = np.array([[cos(ty), 0, -sin(ty)], [0, 1, 0], [sin(ty), 0, cos(ty)]])
    Rz = np.array([[cos(tz), -sin(tz), 0], [sin(tz), cos(tz), 0], [0,0,1]])
    return np.dot(Rx, np.dot(Ry, Rz))

# read xyz coordinate file for the molecule to be solvated and move center to origin
filename = "molecule.xyz"
molecule_atoms = []
molecule_coords = []
with open(filename) as xyz:
    molecule_n_atoms = int(xyz.readline())
    title = xyz.readline()
    for line in xyz:
        atom,x,y,z = line.split()
        molecule_atoms.append(atom)
        molecule_coords.append((float(x), float(y), float(z)))

maxval = np.amax(molecule_coords,0)
minval = np.amin(molecule_coords,0)
delx = -(maxval[0] - (maxval[0]-minval[0])/2)
dely = -(maxval[1] - (maxval[1]-minval[1])/2)
delz = -(maxval[2] - (maxval[2]-minval[2])/2)

molecule_coords = [ (c[0] + delx, c[1] + dely, c[2] + delz) for c in molecule_coords ]
maxval = np.amax(molecule_coords,0)
minval = np.amin(molecule_coords,0)
print('maxval with molecule at origin',maxval)
print('minval with molecule at origin',minval)

# decide on the size of the simulation cell
xparam = 2*(maxval[0]+layer)
yparam = 2*(maxval[1]+layer)
zparam = 2*(maxval[2]+layer)
xcenter = xparam/2
ycenter = yparam/2
zcenter = zparam/2
print('center of the cell',xcenter,ycenter,zcenter)

# move the molecule to the center of the cell
molecule_coords = [ (c[0] + xcenter, c[1] + ycenter, c[2] + zcenter) for c in molecule_coords ]
maxval = np.amax(molecule_coords,0)
minval = np.amin(molecule_coords,0)
print('maxval with molecule at center of the cell',maxval)
print('minval with molecule at center of the cell',minval)

print('box dimensions',xparam,yparam,zparam)
print('center of cell',xcenter,ycenter,zcenter)

# read xyz coordinate file for the solvent
filename = "solvent.xyz"
solvent_atoms = []
solvent_coords = []
with open(filename) as xyz:
    solv_n_atoms = int(xyz.readline())
    title = xyz.readline()
    for line in xyz:
        atom,x,y,z = line.split()
        solvent_atoms.append(atom)
        solvent_coords.append((float(x), float(y), float(z)))

# create an array of possible points inside the box on which to place solvent molecules
xdots = np.arange(0, xparam, step_size)
ydots = np.arange(0, yparam, step_size)
zdots = np.arange(0, zparam, step_size)

maxvaldots = np.amax(xdots,0)
minvaldots = np.amin(xdots,0)
print('maxval dots in x dimension',maxvaldots)
print('minval dots in x dimension',minvaldots)


dots = []
for x in xdots:
    for y in ydots:
        for z in zdots:
            dots.append((float(x), float(y), float(z)))
            print('dot',x,y,z)


good_dots = []
for d in dots:
    smallest_d = big_start_value
    for c in molecule_coords:
        dist = sqrt((c[0] - d[0])**2 + (c[1] - d[1])**2 + (c[2] - d[2])**2)
        smallest_d = min(smallest_d, dist)
    if smallest_d > dmin:
        good_dots.append(d)
        print('adding this dot',d)

# create new array to hold the solvated molecule as it is created
solvated_molecule_atoms = []
solvated_molecule_coords = []
for i, c in enumerate(molecule_coords):
    solvated_molecule_atoms.append(molecule_atoms[i])
    solvated_molecule_coords.append(c)

# for each good point, rotate a solvent molecule and place its center on the point
# must evaluate for collisions
for c in good_dots:
    # calculate the rotation matrix and rotate a water molecule
    rotated_solvent_atoms = []
    rotated_solvent_coords = []
    theta = 2*pi*random(), 2*pi*random(), 2*pi*random()
    matr = rotation(theta)
    for i,d in enumerate(solvent_coords):
        rotated_solvent_atoms.append(solvent_atoms[i])
        rotated_solvent_coords.append(np.dot(matr,d))
    # place the water molecule on the good coord
    rotated_solvent_coords = [ (c[0] - d[0], c[1] - d[1], c[2] - d[2] ) for d in rotated_solvent_coords ]
    # test for collisions
    # s_dist = smallest_dance(rotated_solvent_coords, c)
    smallest_d = big_start_value
    for a in solvated_molecule_coords:
        for b in rotated_solvent_coords:
            dist = sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2)
            smallest_d = min(smallest_d, dist)
    if smallest_d > dmin:
        for i,f in enumerate(rotated_solvent_coords):
            solvated_molecule_atoms.append(rotated_solvent_atoms[i])
            solvated_molecule_coords.append(f)
# write output xyz coordinate file
n_atoms = len(solvated_molecule_coords)
filename = "out.xyz"
with open(filename, 'w') as xyz:
    xyz.write("%s\n" % n_atoms)
    xyz.write("energy\n")
    for i, c in enumerate(solvated_molecule_coords):
        xyz.write("%s %s %s %s\n" %(solvated_molecule_atoms[i], c[0], c[1], c[2]))
print('finished the program')
