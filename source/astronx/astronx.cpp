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


#include<iostream>
#include "configuration.h"

int main(int argc, char *argv[])
{
    if (Astronx::Configuration::get().init(argc, argv) != 0)
        return -1;

    std::ofstream &out = Astronx::Configuration::get().outputFile();

    out << "/-------------------------------------------------------------------------------\\\n";
    out << "|                                  ** AstronX **                                |\n";
    out << "|                                                                               |\n";
    out << "|                A program for the simulation of celestial mechanics            |\n";
    out << "|      /\\                                                              \\\\    // |\n";
    out << "|     //\\\\                Copyright 2012-2018 Jan von Cosel             \\\\  //  |\n";
    out << "|    //  \\\\                                                              \\\\//   |\n";
    out << "|   //====\\\\                  Astronx is free software.                  //\\\\   |\n";
    out << "|  //      \\\\    You can redistribute it and/or modify it under the     //  \\\\  |\n";
    out << "| //        \\\\   terms of the GNU General Public License as published  //    \\\\ |\n";
    out << "|                by the Free Software Foundation, either version 3 of           |\n";
    out << "|                the License, or (at your option) any later version.            |\n";
    out << "|                     Astronx comes with absolutely no warranty.                |\n";
    out << "\\-------------------------------------------------------------------------------/\n";


}

