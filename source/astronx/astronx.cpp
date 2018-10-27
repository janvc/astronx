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
#include <iomanip>
#include "configuration.h"
#include "propagator.h"

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

    if (Astronx::Configuration::get().Verbose())
        std::cout << "AstronX: A program for the simulation of celestial mechanics.\n"
                  << "Copyright 2012-2018 Jan von Cosel. Astronx is free software.\n";

    Astronx::Propagator prop;

    out << "\n\n\n\n\n";
    out << "------------------------------------\n";
    out << "GENERAL INFORMATION ABOUT THE SYSTEM\n";
    out << "------------------------------------\n\n";
    out << "Total mass:\n";
    out << " m = " << std::scientific << std::setprecision(9) << std::setw(15) << Astronx::Configuration::get().TotMass() << " kg\n\n";

    std::array<double,3> cog = prop.com();
    out << "Location of the centre of gravity:\n";
    out << " x = " << std::setprecision(8) << std::setw(15) << cog[0] << " m\n";
    out << " y = " << std::setprecision(8) << std::setw(15) << cog[1] << " m\n";
    out << " z = " << std::setprecision(8) << std::setw(15) << cog[2] << " m\n\n";

    std::array<double,3> lm = prop.linMom();
    out << "Total linear momentum:\n";
    out << " x = " << std::setprecision(8) << std::setw(15) << lm[0] << " kg*m/s\n";
    out << " y = " << std::setprecision(8) << std::setw(15) << lm[1] << " kg*m/s\n";
    out << " z = " << std::setprecision(8) << std::setw(15) << lm[2] << " kg*m/s\n\n";

    std::array<double,3> am = prop.angMom();
    out << "Total angular momentum:\n";
    out << " x = " << std::setprecision(8) << std::setw(15) << am[0] << " kg*m^2/s\n";
    out << " y = " << std::setprecision(8) << std::setw(15) << am[1] << " kg*m^2/s\n";
    out << " z = " << std::setprecision(8) << std::setw(15) << am[2] << " kg*m^2/s\n\n\n\n";

    Astronx::Configuration::get().listParas();

    if (Astronx::Configuration::get().ShiftCOM())
    {
        prop.shiftCom();
        cog = prop.com();

        out << "--------------------------------------------\n";
        out << "SHIFTING THE CENTRE OF GRAVITY TO THE ORIGIN\n";
        out << "--------------------------------------------\n\n";
        out << "Centre of gravity after shifting:\n";
        out << " x = " << std::setprecision(8) << std::setw(15) << cog[0] << " m\n";
        out << " y = " << std::setprecision(8) << std::setw(15) << cog[1] << " m\n";
        out << " z = " << std::setprecision(8) << std::setw(15) << cog[2] << " m\n\n\n\n";
    }

    if (Astronx::Configuration::get().ShiftMom())
    {
        prop.shiftMom();
        lm = prop.linMom();

        out << "--------------------------------------------\n";
        out << "ELIMINATING THE TOTAL MOMENTUM OF THE SYSTEM\n";
        out << "--------------------------------------------\n\n";
        out << "Residual linear momentum:\n";
        out << " x = " << std::setprecision(8) << std::setw(15) << lm[0] << " kg*m/s\n";
        out << " y = " << std::setprecision(8) << std::setw(15) << lm[1] << " kg*m/s\n";
        out << " z = " << std::setprecision(8) << std::setw(15) << lm[2] << " kg*m/s\n\n\n\n";
    }

    prop.propagate();
}

