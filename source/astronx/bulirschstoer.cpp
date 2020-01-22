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
#include <cstring>
#include "propagator.h"
#include "configuration.h"
#include "bulirschstoer.h"
//#include "acceleration.h"


namespace Astronx
{

BulirschStoer::BulirschStoer(const int Npad, System *sys)
    : Propagator(Npad, sys)
{
    m_N_ok = 0;
    m_N_fail = 0;

    void *dm2, *dm3, *dm4, *dm5, *dm6, *dm7, *dm8, *dm9, *dm10, *dm11;
    posix_memalign(&dm2, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm3, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm4, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm5, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm6, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm7, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm8, 64, Configuration::get().MaxSubStep() * 6 * m_Npad * sizeof(double));
    posix_memalign(&dm9, 64, 6 * m_Npad * sizeof(double));
    posix_memalign(&dm10, 64, 6 * m_Npad * sizeof(double));
    posix_memalign(&dm11, 64, 6 * m_Npad * sizeof(double));
    m_x_SubStep = (double*) dm2;
    m_v_SubStep = (double*) dm3;
    m_x_SubFin = (double*) dm4;
    m_v_SubFin = (double*) dm5;
    m_a_BSStart = (double*) dm6;
    m_a_SubInt = (double*) dm7;
    m_extD = (double*) dm8;
    m_extErr = (double*) dm9;
    m_tmpDat = (double*) dm10;
    m_extC = (double*) dm11;

    m_extH.resize(Configuration::get().MaxSubStep());

    std::ofstream &out = Configuration::get().outputFile();

    out << "                          ***************************************\n";
    out << "                          * USING THE BULIRSCH-STOER INTEGRATOR *\n";
    out << "                          ***************************************\n\n";

    out << "       elapsed time      large steps      BS steps      small steps      cpu time [ms]\n";
    out << "                         good    bad\n";

    out.flush();
    if (Configuration::get().Verbose())
    {
        std::cout << " Starting propagation with the Bulirsch-Stoer integrator\n";
    }

    if (Configuration::get().Steps())
    {
        m_stepsFile = std::ofstream(Configuration::get().baseName() + ".stp");
    }
}

BulirschStoer::~BulirschStoer()
{
}

double BulirschStoer::largeStep(double *x, double *v)
{
    m_x_BSLtmp = x;
    m_v_BSLtmp = v;

    if (Configuration::get().Steps())
    {
        m_stepsFile << "# elapsed time:" << std::scientific << std::setprecision(9) << std::setw(20) << std::uppercase << m_sys->elapsedTime() << "\n";
        m_stepsFile << "# trying propagation with step:" << std::setw(20) << m_timeStep << "\n";
    }

    double internalElapsedTime = 0.0;

    while (true)
    {
        if (! Configuration::get().UnRes())
        {
            if (internalElapsedTime + m_timeStep + Configuration::get().minStep()
                    > Configuration::get().tout())
            {
                m_timeStep = Configuration::get().tout() - internalElapsedTime;
            }
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
    return internalElapsedTime;
}

bool BulirschStoer::oneStep()
{
    acceleration(m_x_BSLtmp, m_a_BSStart);
    double gyr = radiusOfGyration(m_x_BSLtmp);
    double v_avg = 0.0;
    for (int i = 0; i < m_Nobj; i++)
    {
        v_avg += std::abs(m_v_BSLtmp[0 * m_Npad + i])
               + std::abs(m_v_BSLtmp[1 * m_Npad + i])
               + std::abs(m_v_BSLtmp[2 * m_Npad + i]);
    }
    v_avg /= (3 * m_Nobj);

    double trialStep = m_timeStep;

    if (Configuration::get().Steps())
    {
        m_stepsFile << "#          Time            Stepsize      Substeps       Error\n";
    }

    bool finished = false;
    bool success = false;
    while (!finished)
    {
        for (int i = 1; i <= Configuration::get().MaxSubStep(); i++)
        {
            subSteps(i, trialStep);

            double extStep = (trialStep / i) * (trialStep / i);

            extrapolate(i, extStep);

            m_delta = 0.0;
            for (int j = 0; j < 3 * m_Nobj; j++)
            {
                m_delta += std::abs(m_extErr[m_Npad * 0 + j] / gyr);
                m_delta += std::abs(m_extErr[m_Npad * 3 + j] / v_avg);
            }
            m_delta /= 6 * m_Nobj;

            m_nsteps = i;

            if (Configuration::get().Steps())
            {
                m_stepsFile << "  " << std::setprecision(8) << std::setw(18) << m_sys->elapsedTime()
                            << std::setprecision(5) << std::setw(17) << trialStep
                            << std::setw(9) << i << std::setw(18) << m_delta << "\n";
            }

            if (m_delta < Configuration::get().Eps())
            {
                finished = true;
                m_doneStep = trialStep;

                memcpy(m_x_BSLtmp, m_x_SubFin, 3 * m_Npad * sizeof(double));
                memcpy(m_v_BSLtmp, m_v_SubFin, 3 * m_Npad * sizeof(double));

                if (trialStep == m_timeStep)
                {
                    success = true;
                }
                break;
            }
        }

        if (!success)
        {
            if (trialStep <= Configuration::get().minStep())
            {
                // TODO: throw a custom exception here
                std::cout << "WARNING: No convergence with minimum stepsize. Aborting!\n";
                return success;
            }

            double factor = std::pow(Configuration::get().Eps() / m_delta, 1.0 / Configuration::get().MaxSubStep());
            factor = std::max(std::min(factor, Configuration::get().RedMin()), Configuration::get().RedMax());

            trialStep *= factor;
            if (trialStep < Configuration::get().minStep())
            {
                trialStep = Configuration::get().minStep();
            }
        }
    }

    return success;
}

void BulirschStoer::subSteps(const int nSteps, const double stepSize)
{
    double smallStep = stepSize / nSteps;

    /*
     * the first step
     */
    for (int i = 0; i < m_Nobj; i++)
    {
        m_x_SubStep[0 * m_Npad + i] = smallStep * (m_v_BSLtmp[0 * m_Npad + i] + 0.5 * smallStep * m_a_BSStart[0 * m_Npad + i]);
        m_x_SubStep[1 * m_Npad + i] = smallStep * (m_v_BSLtmp[1 * m_Npad + i] + 0.5 * smallStep * m_a_BSStart[1 * m_Npad + i]);
        m_x_SubStep[2 * m_Npad + i] = smallStep * (m_v_BSLtmp[2 * m_Npad + i] + 0.5 * smallStep * m_a_BSStart[2 * m_Npad + i]);
        m_x_SubFin[0 * m_Npad + i] = m_x_BSLtmp[0 * m_Npad + i] + m_x_SubStep[0 * m_Npad + i];
        m_x_SubFin[1 * m_Npad + i] = m_x_BSLtmp[1 * m_Npad + i] + m_x_SubStep[1 * m_Npad + i];
        m_x_SubFin[2 * m_Npad + i] = m_x_BSLtmp[2 * m_Npad + i] + m_x_SubStep[2 * m_Npad + i];
        m_a_SubInt[0 * m_Npad + i] = m_a_BSStart[0 * m_Npad + i];
        m_a_SubInt[1 * m_Npad + i] = m_a_BSStart[1 * m_Npad + i];
        m_a_SubInt[2 * m_Npad + i] = m_a_BSStart[2 * m_Npad + i];
    }

    /*
     * the remaining steps
     */
    for (int i = 2; i <= nSteps; i++)
    {
        acceleration(m_x_SubFin, m_a_SubInt);

        /*
         * accumulate the step to reduce roundoff error according to p. 726 Numerical recipes in Fortran 77
         */
        for (int j = 0; j < m_Nobj; j++)
        {
            m_x_SubStep[0 * m_Npad + j] += smallStep * smallStep * m_a_SubInt[0 * m_Npad + j];
            m_x_SubStep[1 * m_Npad + j] += smallStep * smallStep * m_a_SubInt[1 * m_Npad + j];
            m_x_SubStep[2 * m_Npad + j] += smallStep * smallStep * m_a_SubInt[2 * m_Npad + j];
            m_x_SubFin[0 * m_Npad + j] += m_x_SubStep[0 * m_Npad + j];
            m_x_SubFin[1 * m_Npad + j] += m_x_SubStep[1 * m_Npad + j];
            m_x_SubFin[2 * m_Npad + j] += m_x_SubStep[2 * m_Npad + j];
        }
    }

    /*
     * calculate the velocity at the end of the intervall
     */
    acceleration(m_x_SubFin, m_a_SubInt);
    for (int i = 0; i < m_Nobj; i++)
    {
        m_v_SubFin[0 * m_Npad + i] = (m_x_SubStep[0 * m_Npad + i] / smallStep) + 0.5 * smallStep * m_a_SubInt[0 * m_Npad + i];
        m_v_SubFin[1 * m_Npad + i] = (m_x_SubStep[1 * m_Npad + i] / smallStep) + 0.5 * smallStep * m_a_SubInt[1 * m_Npad + i];
        m_v_SubFin[2 * m_Npad + i] = (m_x_SubStep[2 * m_Npad + i] / smallStep) + 0.5 * smallStep * m_a_SubInt[2 * m_Npad + i];
    }
}

void BulirschStoer::extrapolate(const int i_est, const double h_est)
{
    m_extH[i_est - 1] = h_est;

    for (int i = 0; i < 3 * m_Npad; i++)
    {
        m_extErr[0 * m_Npad + i] = m_x_SubFin[i];
        m_extErr[3 * m_Npad + i] = m_v_SubFin[i];
        m_tmpDat[0 * m_Npad + i] = m_x_SubFin[i];
        m_tmpDat[3 * m_Npad + i] = m_v_SubFin[i];
    }

    double extC[6 * m_Npad];

    if (i_est == 1)
    {
        for (int i = 0; i < 3 * m_Npad; i++)
        {
            m_extD[0 * m_Npad + i] = m_x_SubFin[i];
            m_extD[3 * m_Npad + i] = m_v_SubFin[i];
        }
    }
    else
    {
        for (int i = 0; i < 3 * m_Npad; i++)
        {
            extC[0 * m_Npad + i] = m_x_SubFin[i];
            extC[3 * m_Npad + i] = m_v_SubFin[i];
        }

        for (int k = 1; k < i_est; k++)
        {
            double delta = 1.0 / (m_extH[i_est - k - 1] - h_est);
            double tmp1 = h_est * delta;
            double tmp2 = m_extH[i_est - k - 1] * delta;

            for (int j = 0; j < 6 * m_Npad; j++)
            {
                double tmp3 = m_extD[(k - 1) * 6 * m_Npad + j];
                m_extD[(k - 1) * 6 * m_Npad + j] = m_extErr[j];
                delta = extC[j] - tmp3;
                m_extErr[j] = tmp1 * delta;
                extC[j] = tmp2 * delta;
                m_tmpDat[j] += m_extErr[j];
            }
        }

        for (int i = 0; i < 6 * m_Npad; i++)
        {
            m_extD[(i_est - 1) * 6 * m_Npad + i] = m_extErr[i];
        }
    }

    std::memcpy(m_x_SubFin, m_tmpDat + m_Npad * 0, 3 * m_Npad * sizeof(double));
    std::memcpy(m_v_SubFin, m_tmpDat + m_Npad * 3, 3 * m_Npad * sizeof(double));
}

void BulirschStoer::writeOutputLine()
{
    if (Configuration::get().Steps())
    {
        m_stepsFile << "# successful / failed steps:" << m_N_ok << m_N_fail << "\n";
        m_stepsFile << "#\n";
        m_stepsFile << "######################################################################\n";
        m_stepsFile << "#\n";
    }
}

}
