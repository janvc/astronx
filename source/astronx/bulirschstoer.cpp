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
#include "acceleration.h"


namespace Astronx
{

BulirschStoer::BulirschStoer(const int Npad)
    : Propagator(Npad)
{
    m_N_ok = 0;
    m_N_fail = 0;

    void *dm0, *dm1, *dm2, *dm3, *dm4, *dm5, *dm6, *dm7;
    posix_memalign(&dm0, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm1, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm2, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm3, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm4, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm5, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm6, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm7, 64, 3 * m_Npad * sizeof(double));
    m_x_BSLtmp = (double*) dm0;
    m_v_BSLtmp = (double*) dm1;
    m_x_SubStep = (double*) dm2;
    m_v_SubStep = (double*) dm3;
    m_x_SubFin = (double*) dm4;
    m_v_SubFin = (double*) dm5;
    m_a_BSStart = (double*) dm6;
    m_a_SubInt = (double*) dm7;

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

    m_x_BSLtmp = x;
    m_v_BSLtmp = v;

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
}

bool BulirschStoer::oneStep()
{
    std::cout << "this is BulirschStoer::largeStep()\n";

    acceleration(m_x_BSLtmp, m_a_BSStart);
    double gyr = radiusOfGyration(m_x_BSLtmp);
    double v_avg = 0.0;
    for (int i = 0; i < Configuration::get().Nobj(); i++)
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
    bool success = true;
    while (!finished)
    {
        for (int i = 1; i <= Configuration::get().MaxSubStep(); i++)
        {
            subSteps(i, trialStep);

            double extStep = (trialStep / i) * (trialStep / i);

            extrapolate(i, extStep);
        }
    }
    // possible loop structure:
    //
    // finished = false;
    // success = false;     // did we converge with the trial step?
    // while (!finished)
    // {
    //     for (i from 1 to maxsubstep)
    //     {
    //         substeps;
    //         extrapolate;
    //         if (error < eps)
    //         {
    //             finished = true;
    //             if (step == trialstep)
    //                 success = true;
    //             break;
    //         }
    //     }
    //     if (!success)
    //         reduce stepsize;
    // }
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

void BulirschStoer::extrapolate(const int stepNum, const double squaredStep)
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
