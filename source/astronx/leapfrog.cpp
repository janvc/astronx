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
#include "configuration.h"
#include "leapfrog.h"


namespace Astronx
{

LeapFrog::LeapFrog(System *sys)
    : Propagator(sys)
{
    void *dm1, *dm2, *dm3;
    posix_memalign(&dm1, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm2, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm3, 64, 3 * m_Npad * sizeof(double));
    m_x = (double*) dm1;
    m_v = (double*) dm2;
    m_a = (double*) dm3;

    m_stepSize = Configuration::get().tout() / Configuration::get().nSteps();
    m_nSteps = 0;

    std::ofstream &out = Configuration::get().outputFile();

    out << "                            *********************************\n";
    out << "                            * USING THE LEAPFROG INTEGRATOR *\n";
    out << "                            *********************************\n\n";

    out << "       elapsed time            steps                                     cpu time [ms]\n\n";

    if (Configuration::get().Verbose())
    {
        std::cout << " Starting propagation with the Leapfrog integrator\n";
    }
}

LeapFrog::~LeapFrog()
{
}

/**
 * @brief LeapFrog::largeStep
 * @param x
 * @param v
 * @return
 */
double LeapFrog::largeStep(double *x, double *v)
{
    m_x = x;
    m_v = v;

    m_internalElapsedTime = 0.0;

    // do the initial half-step for the velocity:
    acceleration(m_x, m_a);
    for (int j = 0; j < 3 * m_Npad; j++)
    {
        m_v[j] += 0.5 * m_a[j] * m_stepSize;
    }

    // do the remaining steps except for the last one
    for (int i = 0; i < Configuration::get().nSteps() - 1; i++)
    {
        // position
        for (int j = 0; j < 3 * m_Npad; j++)
        {
            m_x[j] += m_v[j] * m_stepSize;
        }

        // velocity
        acceleration(m_x, m_a);
        for (int j = 0; j < 3 * m_Npad; j++)
        {
            m_v[j] += m_a[j] * m_stepSize;
        }

        m_internalElapsedTime += m_stepSize;
        m_nSteps++;
    }

    // do one more full position step ...
    for (int j = 0; j < 3 * m_Npad; j++)
    {
        m_x[j] += m_v[j] + m_stepSize;
    }

    // ... and one half velocity step
    acceleration(m_x, m_a);
    for (int j = 0; j < 3 * m_Npad; j++)
    {
        m_v[j] += 0.5 * m_a[j] * m_stepSize;
    }

    m_internalElapsedTime += m_stepSize;
    m_nSteps++;

    return m_internalElapsedTime;
}

void LeapFrog::writeOutputLine(const double cpuTimeUsed)
{
    std::ofstream &out = Configuration::get().outputFile();

    out << "       " << std::scientific << std::setprecision(5) << std::setw(11) << m_sys->elapsedTime()
        << "        " << std::setw(10) << Configuration::get().nSteps()
        << "                                      "
        << std::fixed << std::setprecision(1) << std::setw(7) << cpuTimeUsed * 1000.0 << std::endl;
}

void LeapFrog::writeSummary()
{
    std::ofstream &out = Configuration::get().outputFile();

    out << std::endl;
    out << "                     *************************************\n";
    out << "                     *    SUMMARY OF THE CALCULATION     *\n";
    out << "                     *                                   *\n";
    out << "                     * total time [s]     Leapfrog steps *\n";
    out << "                     *                                   *\n";
    out << "                     *" << std::scientific << std::setprecision(3) << std::setw(10) << m_sys->elapsedTime()
        << "         " << std::setw(15) << m_nSteps << " *\n";
    out << "                     *************************************\n";
    out << std::endl;
}

}
