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


#include "propagator.h"
#include "configuration.h"
#include "leapfrog.h"


namespace Astronx
{

LeapFrog::LeapFrog(System *sys)
    : Propagator(sys)
{
    void *dm1, *dm3, *dm5;
    void *dm2, *dm4;
    posix_memalign(&dm1, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm2, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm3, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm4, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm5, 64, 3 * m_Npad * sizeof(double));
    m_xLFtmp = (double*) dm1;
    m_vLFtmp = (double*) dm2;
    m_v12 = (double*) dm3;
    m_a0 = (double*) dm4;
    m_x2 = (double*) dm5;

    m_stepSize = Configuration::get().tout() / Configuration::get().nSteps();

    std::ofstream &out = Configuration::get().outputFile();

    out << "                         *********************************\n";
    out << "                         * USING THE LEAPFROG INTEGRATOR *\n";
    out << "                         *********************************\n\n";

    out << "       elapsed time       steps      cpu time [ms]\n";
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
    m_xLFtmp = x;
    m_vLFtmp = v;

    m_internalElapsedTime = 0.0;

    // do the initial half-step for the velocity:
    acceleration(m_xLFtmp, m_a0);
    for (int j = 0; j < 3 * m_Npad; j++)
    {
        m_v12[j] = m_vLFtmp[j] + 0.5 * m_stepSize * m_a0[j];
    }

    // do the remaining steps except for the last one
    for (int i = 0; i < Configuration::get().nSteps() - 1; i++)
    {
        // position
        for (int j = 0; j < 3 * m_Npad; j++)
        {
            m_xLFtmp[j] += m_v12[j] * m_stepSize;
        }

        // velocity
        for (int j = 0; j < 3 * m_Npad; j++)
        {
        }
    }
}

void LeapFrog::writeOutputLine(const double cpuTimeUsed)
{
}

void LeapFrog::writeSummary()
{
}

}
