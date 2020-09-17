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
#include "stepsizeunderflow.h"


namespace Astronx
{

BulirschStoer::BulirschStoer(System *sys)
    : Propagator(sys)
{
    m_NlargeOkTotal = 0;
    m_NlargeFailTotal = 0;
    m_NBSStepsTotal = 0;
    m_NsmallStepsTotal = 0;

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
    m_NlargeOk = 0;
    m_NlargeFail = 0;
    m_NBSSteps = 0;
    m_NsmallSteps = 0;

    m_x_BSLtmp = x;
    m_v_BSLtmp = v;

    if (Configuration::get().Steps())
    {
        m_stepsFile << "# elapsed time:" << std::scientific << std::setprecision(9) << std::setw(20) << std::uppercase << m_sys->elapsedTime() << "\n";
        m_stepsFile << "# trying propagation with step:" << std::setw(19) << m_timeStep << "\n";
    }

    m_internalElapsedTime = 0.0;

    while (true)
    {
        if (! Configuration::get().UnRes())
        {
            if (m_internalElapsedTime + m_timeStep + Configuration::get().minStep()
                    > Configuration::get().tout())
            {
                m_timeStep = Configuration::get().tout() - m_internalElapsedTime;
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
            m_NlargeOk++;
        }
        else
        {
            m_NlargeFail++;
            m_timeStep = m_doneStep;
        }

        m_internalElapsedTime += m_doneStep;

        if (m_internalElapsedTime >= Configuration::get().tout())
        {
            break;
        }
    }

    m_NlargeOkTotal += m_NlargeOk;
    m_NlargeFailTotal += m_NlargeFail;
    m_NBSStepsTotal += m_NBSSteps;
    m_NsmallStepsTotal += m_NsmallSteps;

    return m_internalElapsedTime;
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
        m_NBSSteps++;

        for (int i = 1; i <= Configuration::get().MaxSubStep(); i++)
        {
            subSteps(i, trialStep);

            double extStep = (trialStep / i) * (trialStep / i);

            extrapolate(i, extStep);

            m_delta = 0.0;
            for (int j = 0; j < 3 * m_Nobj; j++)
            {
                m_delta += std::abs(m_extErr[m_Npad * 0 + j] / (gyr != 0.0 ? gyr : 1.0));
                m_delta += std::abs(m_extErr[m_Npad * 3 + j] / (v_avg != 0.0 ? v_avg : 1.0));
            }
            m_delta /= 6 * m_Nobj;

            m_nsteps = i;

            if (Configuration::get().Steps())
            {
                m_stepsFile << "  " << std::setprecision(8) << std::setw(18) << m_sys->elapsedTime() + m_internalElapsedTime
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

        if (Configuration::get().Steps())
        {
            m_stepsFile.flush();
        }

        if (!success)
        {
            if (trialStep <= Configuration::get().minStep())
            {
                throw StepsizeUnderflow();
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

    m_NsmallSteps += nSteps;
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
            m_extC[0 * m_Npad + i] = m_x_SubFin[i];
            m_extC[3 * m_Npad + i] = m_v_SubFin[i];
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
                delta = m_extC[j] - tmp3;
                m_extErr[j] = tmp1 * delta;
                m_extC[j] = tmp2 * delta;
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

void BulirschStoer::writeOutputLine(const double cpuTimeUsed)
{
    std::ofstream &out = Configuration::get().outputFile();

    out << "       " << std::scientific << std::setprecision(5) << std::setw(11) << m_sys->elapsedTime()
        << "      " << std::setw(5) << m_NlargeOk << "  " << std::setw(5) << m_NlargeFail
        << "         " << std::setw(5) << m_NBSSteps << "         " << std::setw(7) << m_NsmallSteps
        << "        " << std::fixed << std::setprecision(1) << std::setw(7) << cpuTimeUsed * 1000.0 << std::endl;

    if (Configuration::get().Steps())
    {
        m_stepsFile << "# successful / failed steps:" << std::setw(5) << m_NlargeOk << "," << std::setw(4) << m_NlargeFail << "\n";
        m_stepsFile << "#\n";
        m_stepsFile << "######################################################################\n";
        m_stepsFile << "#\n";
        m_stepsFile.flush();
    }
}

void BulirschStoer::writeSummary()
{
    std::ofstream &out = Configuration::get().outputFile();

    out << std::endl;
    out << "       *****************************************************************\n";
    out << "       *                  SUMMARY OF THE CALCULATION                   *\n";
    out << "       *                                                               *\n";
    out << "       * total time [s]     large steps      BS steps      small steps *\n";
    out << "       *                    good    bad                                *\n";
    out << "       *" << std::scientific << std::setprecision(3) << std::setw(10) << m_sys->elapsedTime()
        << "      " << std::setw(8) << m_NlargeOkTotal << std::setw(7) << m_NlargeFailTotal
        << "      " << std::setw(8) << m_NBSStepsTotal << "      " << std::setw(10) << m_NsmallStepsTotal << "  *\n";
    out << "       *****************************************************************\n";
    out << std::endl;

    if (Configuration::get().Steps())
    {
        m_stepsFile << "# Simulation done!" << std::endl;
    }
}

}
