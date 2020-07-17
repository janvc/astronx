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

# call this script like
#
#    $ showtrj <name-directory> <viewing axis>
#

import sys
import math

def gyrate(Nobj, frame):
    gr = 0.0
    for j in range(Nobj):
        vecLen = math.sqrt(float(frame[3 * j + 1 + 0]) * float(frame[3 * j + 1 + 0]) \
                         + float(frame[3 * j + 1 + 1]) * float(frame[3 * j + 1 + 1]) \
                         + float(frame[3 * j + 1 + 2]) * float(frame[3 * j + 1 + 2]))
        gr += vecLen

    gr /= float(Nobj)
    return gr

if __name__ == "__main__":
    trjName = sys.argv[1]
    viewAxis = sys.argv[2].lower()

    # open the trajectory file and find the
    # maximum values to set the plot boundaries
    with open(trjName) as trjFile:
        tmpData = trjFile.readlines()

    Nobj = int((len(tmpData[4].split()) - 1) / 3)
    Nframes = len(tmpData) - 4
    print(Nframes)

    # initialize the boundaries with the position of the 1st object
    xMin = float(tmpData[4].split()[1])
    xMax = float(tmpData[4].split()[1])
    yMin = float(tmpData[4].split()[2])
    yMax = float(tmpData[4].split()[2])
    zMin = float(tmpData[4].split()[3])
    zMax = float(tmpData[4].split()[3])

    for i in range(5, len(tmpData)):
        frame = tmpData[i].split()
        grCurr = gyrate(Nobj, frame)
        if grCurr > abs(xMin):
            xMin = -grCurr
        if grCurr > xMax:
            xMax = grCurr
        if grCurr > abs(yMin):
            yMin = -grCurr
        if grCurr > yMax:
            yMax = grCurr
        if grCurr > abs(zMin):
            zMin = -grCurr
        if grCurr > zMax:
            zMax = grCurr


    # create the plot file
    with open("tmpplt.plt", "w") as plotFile:
        plotFile.write("reset\n")
        plotFile.write("set term x11\n")
        plotFile.write("unset key\n")

        if viewAxis == "x":
            plotFile.write("set xrange [{0:8.1e}:{1:8.1e}]\n".format(yMin, yMax))
            plotFile.write("set yrange [{0:8.1e}:{1:8.1e}]\n".format(zMin, zMax))
        elif viewAxis == "y":
            plotFile.write("set xrange [{0:8.1e}:{1:8.1e}]\n".format(xMin, xMax))
            plotFile.write("set yrange [{0:8.1e}:{1:8.1e}]\n".format(zMin, zMax))
        elif viewAxis == "z":
            plotFile.write("set xrange [{0:8.1e}:{1:8.1e}]\n".format(xMin, xMax))
            plotFile.write("set yrange [{0:8.1e}:{1:8.1e}]\n".format(yMin, yMax))

        plotFile.write("do for [ii=1:{0}] {{\n".format(Nframes - 1))
        plotFile.write('    set title sprintf("frame %i", ii)\n')
        if viewAxis == "x":
            plotFile.write('    p for [jj=1:{0}] "{1}" every ::ii::ii u (column(3*(jj-1)+3)):(column(3*(jj-1)+4)) w p pt 7 ps 2\n'.format(Nobj, trjName))
        elif viewAxis == "y":
            plotFile.write('    p for [jj=1:{0}] "{1}" every ::ii::ii u (column(3*(jj-1)+2)):(column(3*(jj-1)+4)) w p pt 7 ps 2\n'.format(Nobj, trjName))
        elif viewAxis == "z":
            plotFile.write('    p for [jj=1:{0}] "{1}" every ::ii::ii u (column(3*(jj-1)+2)):(column(3*(jj-1)+3)) w p pt 7 ps 2\n'.format(Nobj, trjName))
        #plotFile.write("    pause 0.04\n")
        plotFile.write("}\n")

