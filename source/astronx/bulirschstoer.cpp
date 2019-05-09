/*
 * Copyright 2015-2019 Jan von Cosel
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
#include "propagator.h"
#include "configuration.h"
#include "bulirschstoer.h"


namespace Astronx
{

BulirschStoer::BulirschStoer()
    : Propagator()
{
    std::ofstream &out = Configuration::get().outputFile();

    out << "                          ***************************************\n";
    out << "                          * USING THE BULIRSCH-STOER INTEGRATOR *\n";
    out << "                          ***************************************\n\n";

    out << "       elapsed time      large steps      BS steps      small steps      cpu time [ms]\n";
    out << "                         good    bad\n";

    if (Configuration::get().Verbose())
        std::cout << " Starting propagation with the Bulirsch-Stoer integrator\n";
}

BulirschStoer::~BulirschStoer()
{
}

void BulirschStoer::largeStep()
{
    std::cout << "this is BulirschStoer::largeStep()\n";

    if (Configuration::get().Steps())
    {
        std::ofstream &stepsFile = Configuration::get().stepsFile();
        stepsFile << "# trying propagation with step:" << std::setw(20) << m_timeStep << "\n";
    }
}

void BulirschStoer::writeOutputLine()
{
    std::cout << "this is BulirschStoer::writeOutputLine()\n";

    if (Configuration::get().Steps())
    {
        std::ofstream &stepsFile = Configuration::get().stepsFile();
        stepsFile << "# successful / failed steps:" << m_N_ok << m_N_fail << "\n";
        stepsFile << "#\n";
        stepsFile << "######################################################################\n";
        stepsFile << "#\n";
    }
}

}
