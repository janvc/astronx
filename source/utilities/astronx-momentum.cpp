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
    bool polar = false;
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "produce this help message")
            ("linear,l", "print linear momentum")
            ("angular,a", "print angular momentum")
            ("polar,p", "use polar coordinates")
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
    polar = vm.count("polar");

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

                pxTot += px;
                pyTot += py;
                pzTot += pz;

                double p = std::sqrt(px * px + py * py + pz * pz);
                double theta = std::acos(pz / p);
                double phi = std::atan2(py, px);

                std::cout << std::setprecision(digits) << std::setw(digits + 8) << (!polar ? px : p)
                          << std::setprecision(digits) << std::setw(digits + 8) << (!polar ? py : theta)
                          << std::setprecision(digits) << std::setw(digits + 8) << (!polar ? pz : phi);

            }

            double pTot = std::sqrt(pxTot * pxTot + pyTot * pyTot + pzTot * pzTot);
            double thetaTot = std::acos(pzTot / pTot);
            double phiTot = std::atan2(pyTot, pxTot);

            std::cout << std::setprecision(digits) << std::setw(digits + 8) << (!polar ? pxTot : pTot)
                      << std::setprecision(digits) << std::setw(digits + 8) << (!polar ? pyTot : thetaTot)
                      << std::setprecision(digits) << std::setw(digits + 8) << (!polar ? pzTot : phiTot);
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

                lxTot += lx;
                lyTot += ly;
                lzTot += lz;

                double l = std::sqrt(lx * lx + ly * ly + lz * lz);
                double theta = std::acos(lz / l);
                double phi = std::atan2(ly, lx);

                std::cout << std::setprecision(digits) << std::setw(digits + 8) << (!polar ? lx : l)
                          << std::setprecision(digits) << std::setw(digits + 8) << (!polar ? ly : theta)
                          << std::setprecision(digits) << std::setw(digits + 8) << (!polar ? lz : phi);
            }

            double lTot = std::sqrt(lxTot * lxTot + lyTot * lyTot + lzTot * lzTot);
            double thetaTot = std::acos(lzTot / lTot);
            double phiTot = std::atan2(lyTot, lxTot);

            std::cout << std::setprecision(digits) << std::setw(digits + 8) << (!polar ? lxTot : lTot)
                      << std::setprecision(digits) << std::setw(digits + 8) << (!polar ? lyTot : thetaTot)
                      << std::setprecision(digits) << std::setw(digits + 8) << (!polar ? lzTot : phiTot);
        }

        std::cout << std::endl;
    }

    return 0;
}
