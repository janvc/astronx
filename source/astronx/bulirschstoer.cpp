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
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include "propagator.h"
#include "configuration.h"
#include "bulirschstoer.h"


namespace Astronx
{

BulirschStoer::BulirschStoer(const int Npad)
    : Propagator(Npad)
{
    m_N_ok = 0;
    m_N_fail = 0;

    void *dm0, *dm1;
    posix_memalign(&dm0, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm1, 64, 3 * m_Npad * sizeof(double));
    m_x_BSLtmp = (double*) dm0;
    m_v_BSLtmp = (double*) dm1;

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

void BulirschStoer::largeStep(double *x, double *v)
{
    std::cout << "this is BulirschStoer::largeStep()\n";

    if (Configuration::get().Steps())
    {
        std::ofstream &stepsFile = Configuration::get().stepsFile();
        stepsFile << "# trying propagation with step:" << std::setw(20) << m_timeStep << "\n";
    }

    double internalElapsedTime = 0.0;

    while (true)
    {
        if (! Configuration::get().UnRes())
        {
            if (internalElapsedTime + m_timeStep + Configuration::get().minStep() > Configuration::get().tout())
                m_timeStep = Configuration::get().tout() - internalElapsedTime;
        }

        bool success = oneStep();

        if (success)
        {
            if (m_nsteps <= Configuration::get().IncThres())
            {
                double factor = std::pow(Configuration::get().Eps() / m_delta, 1.0 / m_nsteps);
                m_timeStep = m_doneStep * std::min(factor, Configuration::get().MaxInc());
            }
            m_N_ok++;
        }
        else
        {
            m_N_fail++;
            m_timeStep = m_doneStep;
        }

        internalElapsedTime += m_doneStep;

        if (internalElapsedTime >= Configuration::get().tout())
        {
            break;
        }
    }
}

bool BulirschStoer::oneStep()
{
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
