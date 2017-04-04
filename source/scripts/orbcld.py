#!/usr/bin/env python3

# Copyright 2012-2017 Jan von Cosel
#
# This file is part of astronx.
#
# astronx if free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# astronx is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have recieved a copy of the GNU General Public License
# along with astronx. If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################
#

import os
import math
import sys
import random
import numpy as np

random.seed()

argsProvided = len(sys.argv)

if argsProvided <= 1:
    print("Mass of central object [kg]")
    Cmass = float(input())
else:
    Cmass = float(sys.argv[1])

if argsProvided <= 4:
    print("X, Y and Z position of cloud center [m]")
    posString = input()
    cldCen = []
    cldCen.append(float(posString.split()[0]))
    cldCen.append(float(posString.split()[1]))
    cldCen.append(float(posString.split()[2]))
else:
    cldCen = []
    cldCen.append(float(sys.argv[2]))
    cldCen.append(float(sys.argv[3]))
    cldCen.append(float(sys.argv[4]))

if argsProvided <= 5:
    print("Maximum cloud object mass [kg]:")
    maxMass = float(input())
else:
    maxMass = float(sys.argv[5])

if argsProvided <= 6:
    print("Cloud radius [m]:")
    cldRad = float(input())
else:
    cldRad = float(sys.argv[6])

if argsProvided <= 7:
    print("Number of objects in cloud:")
    Nobj = int(input())
else:
    Nobj = int(sys.argv[7])

if argsProvided <= 8:
    print("Output file:")
    outputFileName = input()
else:
    outputFileName = sys.argv[8]

# create the masses:
masses = []
for i in range(Nobj):
    randNum = 0.0
    while (randNum == 0.0):
        randNum = random.random()
    # take the log of the uniform random number to
    # get an exponential probability distribution:
    randNum = -math.log(randNum)
    masses.append(maxMass * randNum)
# scale the masses so that the largest mass corresponds
# to the maximum allowed mass:
largestMass = 0.0
for i in range(Nobj):
    if (masses[i] > largestMass):
        largestMass = masses[i]
for i in range(Nobj):
    masses[i] = masses[i] * maxMass / largestMass

x = []
y = []
z = []
vx = []
vy = []
vz = []
for i in range(Nobj):
    r = cldRad * random.random()
    theta = math.pi * random.random()
    phi = 2.0 * math.pi * random.random()
    xobj = cldCen[0] + (r * math.sin(theta) * math.cos(phi))
    yobj = cldCen[1] + (r * math.sin(theta) * math.sin(phi))
    zobj = cldCen[2] + (r * math.cos(theta))
    x.append(xobj)
    y.append(yobj)
    z.append(zobj)
    dist = math.sqrt((xobj - cldCen[0])**2 + (yobj - cldCen[1])**2 + (zobj - cldCen[2])**2)
    vabs = math.sqrt(6.6726e-11 * Cmass / dist)
    vDirec = np.cross([xobj, yobj, zobj], [0, 0, 1])
    vDirec /= np.linalg.norm(vDirec)
    vx.append(vDirec[0] * vabs)
    vy.append(vDirec[1] * vabs)
    vz.append(vDirec[2] * vabs)

outputFile = open(outputFileName, "w")
outputFile.write("    central_obj")
outputFile.write("{0:13.4e}".format(Cmass))
outputFile.write("{0:14.4e} {1:11.4e} {2:11.4e}".format(0.0, 0.0, 0.0))
outputFile.write("{0:14.4e} {1:11.4e} {2:11.4e}\n".format(0.0, 0.0, 0.0))
for i in range(Nobj):
    outputFile.write("    obj_{num:07d}".format(num=i+1))
    outputFile.write("{0:13.4e}".format(masses[i]))
    outputFile.write("{0:14.4e} {1:11.4e} {2:11.4e}".format(x[i], y[i], z[i]))
    outputFile.write("{0:14.4e} {1:11.4e} {2:11.4e}\n".format(vx[i], vy[i], vz[i]))

outputFile.close()





