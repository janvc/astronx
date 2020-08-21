/*
 * Copyright 2015-2020 Jan von Cosel
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
#include <fstream>
#include <vector>
#include <boost/program_options.hpp>
#include "trajectory.h"

int main(int argc, char *argv[])
{
    std::string trjFileName;
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "produce this help message")
            ("input,i", po::value<std::string>(&trjFileName), "the trajectory file")
            ;

    po::positional_options_description pos_desc;
    pos_desc.add("input", -1);
    po::command_line_parser parser(argc, argv);
    parser.options(desc).positional(pos_desc).allow_unregistered();
    po::parsed_options options = parser.run();
    po::variables_map vm;
    po::store(options, vm);
    po::notify(vm);

    if (vm.count("input"))
    {
        trjFileName = vm["input"].as<std::string>();
    }

    Trajectory trj(trjFileName);
    int Nobj = trj.nObj();
    int nFrames = trj.nFrames();
    std::vector<double> masses = trj.readMasses();

    /*
     * iterate over the frames
     */
    for (int i = 0; i < nFrames; i++)
    {
        TrajectoryFrame frame = trj.readNextFrame();
        double elapsedTime = frame.getElapsedTime();
        std::vector<double> x = frame.getX();
        std::vector<double> v = frame.getV();

        /*
         * calculate the potential energy
         */
        double ePot = 0.0;
        for (int j = 0; j < Nobj; j++)
        {
            for (int k = j + 1; k < Nobj; k++)
            {
                double dx = x[3 * k + 0] - x[3 * j + 0];
                double dy = x[3 * k + 1] - x[3 * j + 1];
                double dz = x[3 * k + 2] - x[3 * j + 2];
                double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

                ePot -= masses[j] * masses[k] / dist;
            }
        }
        ePot *= 6.6726e-11;

        /*
         * calculate the kinetic energy
         */
        double eKin = 0.0;
        for (int j = 0; j < Nobj; j++)
        {
            eKin += 0.5 * masses[j]
                    * (v[3 * j + 0] * v[3 * j + 0]
                     + v[3 * j + 1] * v[3 * j + 1]
                     + v[3 * j + 2] * v[3 * j + 2]);
        }

        std::cout << elapsedTime << " " << ePot << " " << eKin << std::endl;
    }

    return 0;
}
