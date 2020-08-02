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
#include <boost/program_options.hpp>

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

    std::cout << "the trajectory file is " << trjFileName << std::endl;

    std::ifstream trjFile = std::ifstream(trjFileName, std::ios::binary);

    double *testData = new double[10];

    trjFile.read(reinterpret_cast<char*>(testData), 10 * sizeof(double));

    for (int i = 0; i < 10; i++)
    {
        std::cout << testData[i] << std::endl;
    }

    return 0;
}
