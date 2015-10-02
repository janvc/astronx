#!/usr/bin/env python

# Copyright 2012-2015 Jan von Cosel
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

random.seed()

print("Enter the number of objects:")
Nobj = int(input(">>> "))

print("Enter the maximum allowable mass [kg]:")
maxMass = float(input(">>> "))

print("Enter the maximum allowable velocity [m/s]:")
maxVel = float(input(">>> "))

print("Choose the mass distribution:")
print(" (1) uniform")
print(" (2) exponential\n")
massDist = int(input(">>> "))

print("Choose the type of cloud to be created:")
print(" (1) Rectangular, cartesian")
print(" (2) Rectangular, polar")
print(" (3) Spherical, cartesian")
print(" (4) Spherical, polar\n")
cloudType = int(input(">>> "))

if (cloudType == 3 or cloudType == 4):
    print("Enter the radius of the sphere in meters:\n")
    radius = float(input(">>> "))
else:
    print ("Enter the boundaries of the box [m]:\n")
    dX = float(input("dX > "))
    dY = float(input("dY > "))
    dZ = float(input("dZ > "))

while True:
    print("Enter the filename to write the results to:\n")
    outputFileName = raw_input(">>> ")
    if (not os.path.isfile(outputFileName)):
        break
    print("The file already exists.")

# create the masses:
masses = []
if (massDist == 1):
    for i in range(Nobj):
        masses.append(maxMass * random.random())
elif (massDist == 2):
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
            largestMass = masses[i]     # maybe there's an easier way to do this...
    for i in range(Nobj):
        masses[i] = masses[i] * maxMass / largestMass

x = []
y = []
z = []
vx = []
vy = []
vz = []
if (cloudType == 1):    # rectangular cartesian
    for i in range(Nobj):
        x.append(-dX + 2.0 * random.random() * dX)
        y.append(-dY + 2.0 * random.random() * dY)
        z.append(-dZ + 2.0 * random.random() * dZ)
        vx.append(-maxVel + 2.0 * random.random() * maxVel)
        vy.append(-maxVel + 2.0 * random.random() * maxVel)
        vz.append(-maxVel + 2.0 * random.random() * maxVel)
elif (cloudType == 2):  # rectangular polar
    while (len(x) < Nobj):
        radius = math.sqrt(dX**2 + dY**2 + dZ**2)
        r = radius * random.random()
        theta = math.pi * random.random()
        phi = 2.0 * math.pi * random.random()
        xTmp = r * math.sin(theta) * math.cos(phi)
        yTmp = r * math.sin(theta) * math.sin(phi)
        zTmp = r * math.cos(theta)
        if (abs(xTmp) <= dX and abs(yTmp) <= dY and abs(zTmp) <= dZ):
            x.append(xTmp)
            y.append(yTmp)
            z.append(zTmp)
            vx.append(-maxVel + 2.0 * random.random() * maxVel)
            vy.append(-maxVel + 2.0 * random.random() * maxVel)
            vz.append(-maxVel + 2.0 * random.random() * maxVel)
elif (cloudType == 3):  # spherical cartesian
    while (len(x) < Nobj):
        xTmp = -radius + 2.0 * radius * random.random()
        yTmp = -radius + 2.0 * radius * random.random()
        zTmp = -radius + 2.0 * radius * random.random()
        if (math.sqrt(xTmp**2 + yTmp**2 + zTmp**2) <= radius):
            x.append(xTmp)
            y.append(yTmp)
            z.append(zTmp)
            vx.append(-maxVel + 2.0 * random.random() * maxVel)
            vy.append(-maxVel + 2.0 * random.random() * maxVel)
            vz.append(-maxVel + 2.0 * random.random() * maxVel)
elif (cloudType == 4):  # spherical polar
    for i in range(Nobj):
        r = radius * random.random()
        theta = math.pi * random.random()
        phi = 2.0 * math.pi * random.random()
        x.append(r * math.sin(theta) * math.cos(phi))
        y.append(r * math.sin(theta) * math.sin(phi))
        z.append(r * math.cos(theta))
        vx.append(-maxVel + 2.0 * random.random() * maxVel)
        vy.append(-maxVel + 2.0 * random.random() * maxVel)
        vz.append(-maxVel + 2.0 * random.random() * maxVel)

outputFile = open(outputFileName, "w")
for i in range(Nobj):
    outputFile.write("    obj_{num:05d}".format(num=i+1))
    outputFile.write("{0:13.3e}".format(masses[i]))
    outputFile.write("{0:14.3e} {1:11.3e} {2:11.3e}".format(x[i], y[i], z[i]))
    outputFile.write("{0:14.3e} {1:11.3e} {2:11.3e}\n".format(vx[i], vy[i], vz[i]))

outputFile.close()





