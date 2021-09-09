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
#include "rungekutta4.h"


namespace Astronx
{

RungeKutta4::RungeKutta4(System *sys)
    : Propagator(sys)
{
    void *dm1, *dm2, *dm3, *dm4, *dm5, *dm6, *dm7, *dm8;
    void *dm9, *dm10, *dm11, *dm12, *dm13, *dm14, *dm15;
    posix_memalign(&dm1, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm2, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm3, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm4, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm5, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm6, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm7, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm8, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm9, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm10, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm11, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm12, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm13, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm14, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm15, 64, 3 * m_Npad * sizeof(double));
    m_x0 = (double*) dm1;
    m_v0 = (double*) dm2;
    m_a0 = (double*) dm3;
    m_x1 = (double*) dm4;
    m_v1 = (double*) dm5;
    m_a1 = (double*) dm6;
    m_x2 = (double*) dm7;
    m_v2 = (double*) dm8;
    m_a2 = (double*) dm9;
    m_x3 = (double*) dm10;
    m_v3 = (double*) dm11;
    m_a3 = (double*) dm12;
    m_x4 = (double*) dm13;
    m_v4 = (double*) dm14;
    m_a4 = (double*) dm14;

    m_stepSize = Configuration::get().tout() / Configuration::get().nSteps();
    m_nSteps = 0;

    std::ofstream &out = Configuration::get().outputFile();

    out << "                    *************************************************\n";
    out << "                    * USING THE FOURTH ORDER RUNGE-KUTTA INTEGRATOR *\n";
    out << "                    *************************************************\n\n";

    out << "       elapsed time            steps                                     cpu time [ms]\n\n";

    if (Configuration::get().Verbose())
    {
        std::cout << " Starting propagation with the Runge-Kutta integrator of fourth order\n";
    }
}

RungeKutta4::~RungeKutta4()
{
}

double RungeKutta4::largeStep(double *x, double *v)
{
    m_x0 = x;
    m_v0 = v;

    m_internalElapsedTime = 0.0;

    for (int i = 0; i < Configuration::get().nSteps(); i++)
    {
        // first substep:
        acceleration(m_x0, m_a0);
        for (int j = 0; j < 3 * m_Npad; j++)
        {
            m_x1[j] = m_x0[j] + 0.5 * m_stepSize * m_v0[j];
            m_v1[j] = m_v0[j] + 0.5 * m_stepSize * m_a0[j];
        }

        // second substep:
        acceleration(m_x1, m_a1);
        for (int j = 0; j < 3 * m_Npad; j++)
        {
            m_x2[j] = m_x0[j] + 0.5 * m_stepSize * m_v1[j];
            m_v2[j] = m_v0[j] + 0.5 * m_stepSize * m_a1[j];
        }

        // third substep:
        acceleration(m_x2, m_a2);
        for (int j = 0; j < 3 * m_Npad; j++)
        {
            m_x3[j] = m_x0[j] + m_stepSize * m_v2[j];
            m_v3[j] = m_v0[j] + m_stepSize * m_a2[j];
        }

        // fourth substep:
        acceleration(m_x3,  m_a3);

        // calculate new values
        for (int j = 0; j < 3 * m_Npad; j++)
        {
            m_x0[j] += m_stepSize * (m_v0[j] + 2.0 * m_v1[j] + 2.0 * m_v2[j] + m_v3[j]) / 6.0;
            m_v0[j] += m_stepSize * (m_a0[j] + 2.0 * m_a1[j] + 2.0 * m_a2[j] + m_a3[j]) / 6.0;
        }

        m_internalElapsedTime += m_stepSize;
        m_nSteps++;
    }

    return m_internalElapsedTime;
}

void RungeKutta4::writeOutputLine(const double cpuTimeUsed)
{
    std::ofstream &out = Configuration::get().outputFile();

    out << "       " << std::scientific << std::setprecision(5) << std::setw(11) << m_sys->elapsedTime()
        << "        " << std::setw(10) << Configuration::get().nSteps()
        << "                                      "
        << std::fixed << std::setprecision(1) << std::setw(7) << cpuTimeUsed * 1000.0 << std::endl;
}

void RungeKutta4::writeSummary()
{
    std::ofstream &out = Configuration::get().outputFile();

    out << std::endl;
    out << "                    ****************************************\n";
    out << "                    *     SUMMARY OF THE CALCULATION       *\n";
    out << "                    *                                      *\n";
    out << "                    * total time [s]     Runge-Kutta steps *\n";
    out << "                    *                                      *\n";
    out << "                    *" << std::scientific << std::setprecision(3) << std::setw(10) << m_sys->elapsedTime()
        << "            " << std::setw(15) << m_nSteps << " *\n";
    out << "                    ****************************************\n";
    out << std::endl;
}

}
