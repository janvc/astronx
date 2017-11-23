/*
 * Copyright 2015 Jan von Cosel
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


#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include "configuration.h"

namespace Astronx
{

void Configuration::init(int argnum, char *arguments[])
{
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce this help message")
        ("input,i", po::value<std::string>(&m_inputFileName)->required(), "the input file")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argnum, arguments, desc), vm);
    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return;
    }
    po::notify(vm);
}

} // namespace astronx
