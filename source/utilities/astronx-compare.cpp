/*
 * Copyright 2015-2021 Jan von Cosel
 *
 * This file is part of astronx.
 *
 * astronx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * astronx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have recieved a copy of the GNU General Public License
 * along with astronx. If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include <string>
#include <iostream>
#include <boost/program_options.hpp>
#include "trajectory.h"

int main(int argc, char *argv[])
{
    std::string trjAFileName;
    std::string trjBFileName;
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "produce this help message")
            ("first,f", po::value<std::string>(&trjAFileName)->required(), "the first trajectory file")
            ("second,s", po::value<std::string>(&trjBFileName)->required(), "the second trajectory file")
            ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 0;
    }

    po::notify(vm);

    Trajectory trjA(trjAFileName);
    Trajectory trjB(trjBFileName);

    int nObjA = trjA.nObj();
    int nObjB = trjB.nObj();
    if (nObjA != nObjB)
    {
        std::cerr << "Only trajectories with an equal number of objects can be compared.\n";
        return -1;
    }

    int nFramesA = trjA.nFrames();
    int nFramesB = trjB.nFrames();

    // compare trajectories only up to the last frame of the shorter one
    int nComparedFrames = std::min(nFramesA, nFramesB);

    return 0;
}
