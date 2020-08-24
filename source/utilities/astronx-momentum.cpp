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
#include <iomanip>
#include <fstream>
#include <vector>
#include <boost/program_options.hpp>
#include "trajectory.h"

int main(int argc, char *argv[])
{
    std::string trjFileName;
    int digits = 8;
    bool printLinearMomentum = false;
    bool printAngularMomentum = false;
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "produce this help message")
            ("linear,l", "produce this help message")
            ("angular,a", "produce this help message")
            ("input,i", po::value<std::string>(&trjFileName), "the trajectory file")
            ("digits,d", po::value<int>(&digits), "number of digits in the output")
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

    printLinearMomentum = vm.count("linear");
    printAngularMomentum = vm.count("angular");

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

        std::cout << std::scientific << std::setprecision(digits) << std::setw(digits + 8) << elapsedTime;

        /*
         * calculate the linear momentum
         */
        if (printLinearMomentum)
        {
            double pxTot = 0.0;
            double pyTot = 0.0;
            double pzTot = 0.0;

            for (int j = 0; j < Nobj; j++)
            {
                double px = v[3 * j + 0] * masses[j];
                double py = v[3 * j + 1] * masses[j];
                double pz = v[3 * j + 2] * masses[j];

                std::cout << std::setprecision(digits) << std::setw(digits + 8) << px
                          << std::setprecision(digits) << std::setw(digits + 8) << py
                          << std::setprecision(digits) << std::setw(digits + 8) << pz;

                pxTot += px;
                pyTot += py;
                pzTot += pz;
            }

            std::cout << std::setprecision(digits) << std::setw(digits + 8) << pxTot
                      << std::setprecision(digits) << std::setw(digits + 8) << pyTot
                      << std::setprecision(digits) << std::setw(digits + 8) << pzTot;
        }



        /*
         * calculate the angular momentum
         */
        if (printAngularMomentum)
        {
            double lxTot = 0.0;
            double lyTot = 0.0;
            double lzTot = 0.0;

            for (int j = 0; j < Nobj; j++)
            {
                double lx = masses[j] * (x[3 * j + 1] * v[3 * j + 2] - x[3 * j + 2] * v[3 * j + 1]);
                double ly = masses[j] * (x[3 * j + 2] * v[3 * j + 0] - x[3 * j + 0] * v[3 * j + 2]);
                double lz = masses[j] * (x[3 * j + 0] * v[3 * j + 1] - x[3 * j + 1] * v[3 * j + 0]);

                std::cout << std::setprecision(digits) << std::setw(digits + 8) << lx
                          << std::setprecision(digits) << std::setw(digits + 8) << ly
                          << std::setprecision(digits) << std::setw(digits + 8) << lz;

                lxTot += lx;
                lyTot += ly;
                lzTot += lz;
            }

            std::cout << std::setprecision(digits) << std::setw(digits + 8) << lxTot
                      << std::setprecision(digits) << std::setw(digits + 8) << lyTot
                      << std::setprecision(digits) << std::setw(digits + 8) << lzTot;
        }

    return 0;
}
