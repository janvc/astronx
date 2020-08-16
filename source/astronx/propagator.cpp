#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <string.h>
#include "propagator.h"
#include "configuration.h"
#include "constants.h"
//#include "acceleration.h"

namespace Astronx
{

Propagator::Propagator()
{
    m_Nobj = Configuration::get().Nobj();
    m_totMass = Configuration::get().TotMass();

    // make sure the array length is divisible by 4 (AVX register length)
    m_Npad = m_Nobj % 4 == 0 ? m_Nobj : ((m_Nobj / 4) + 1) * 4;

    void *dm0, *dm1, *dm2;
    posix_memalign(&dm0, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm1, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm2, 64, m_Npad * sizeof(double));
    m_x_BSlarge = (double*) dm0;
    m_v_BSlarge = (double*) dm1;
    m_masses = (double*) dm2;

    void *dm3, *dm4, *dm5, *dm6, *dm7, *dm8, *dm9;
    posix_memalign(&dm3, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm4, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm5, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm6, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm7, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm8, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm9, 64, 3 * m_Npad * sizeof(double));
    m_x_BSLtmp = (double*) dm3;
    m_v_BSLtmp = (double*) dm4;
    m_a_BSstart = (double*) dm5;
    m_x_SubStep = (double*) dm6;
    m_x_SubFin = (double*) dm7;
    m_v_SubFin = (double*) dm8;
    m_a_SubInt = (double*) dm9;

    for (int i = 0; i < m_Nobj; i++)
    {
        m_x_BSlarge[0 * m_Npad + i] = Configuration::get().XX0()[i];
        m_x_BSlarge[1 * m_Npad + i] = Configuration::get().XY0()[i];
        m_x_BSlarge[2 * m_Npad + i] = Configuration::get().XZ0()[i];
        m_v_BSlarge[0 * m_Npad + i] = Configuration::get().VX0()[i];
        m_v_BSlarge[1 * m_Npad + i] = Configuration::get().VY0()[i];
        m_v_BSlarge[2 * m_Npad + i] = Configuration::get().VZ0()[i];
        m_masses[i] = Configuration::get().masses()[i];
        m_names = Configuration::get().names();
    }

    /*
     * The format of the binary trajectory is as follows:
     * - number of objects (int)
     * - vector with masses of all objects (double)
     * - for each frame:
     *   - elapsed time (double)
     *   - vector with positions of all objects (double)
     *   - vector with velocities of all objects (double)
     */
    m_binTrj = std::ofstream(Configuration::get().baseName() + ".bin.trj", std::ios::binary);
    m_binTrj.write((char*) &m_Nobj, sizeof(int));
    for (int i = 0; i < m_Nobj; i++)
    {
        m_binTrj.write((char*) &m_masses[i], sizeof(double));
    }

    if (Configuration::get().TextTrj())
    {
        m_txtTrj = std::ofstream(Configuration::get().baseName() + ".txt.trj");

        m_txtTrj << "# trajectory in gnuplot-friendly text form\n";
        m_txtTrj << "#\n";
        m_txtTrj << "# time               ";
        for (int i = 0; i < m_Nobj - 1; i++)
        {
            m_txtTrj << m_names[i];
            for (int j = 0; j < 3 * Configuration::get().Ndigit() + 25 - m_names[i].size(); j++)
            {
                m_txtTrj << " ";
            }
        }
        m_txtTrj << std::setw(30) << m_names[m_Nobj - 1] << "\n";
        m_txtTrj << "#                    ";
        for (int i = 0; i < m_Nobj - 1; i++)
        {
            m_txtTrj << "x [m]";
            for (int j = 0; j < Configuration::get().Ndigit() + 3; j++)
            {
                m_txtTrj << " ";
            }
            m_txtTrj << "y [m]";
            for (int j = 0; j < Configuration::get().Ndigit() + 3; j++)
            {
                m_txtTrj << " ";
            }
            m_txtTrj << "z [m]";
            for (int j = 0; j < Configuration::get().Ndigit() + 4; j++)
            {
                m_txtTrj << " ";
            }
        }
        m_txtTrj << "x [m]";
        for (int j = 0; j < Configuration::get().Ndigit() + 3; j++)
        {
            m_txtTrj << " ";
        }
        m_txtTrj << "y [m]";
        for (int j = 0; j < Configuration::get().Ndigit() + 3; j++)
        {
            m_txtTrj << " ";
        }
        m_txtTrj << "z [m]\n";
    }

    if (Configuration::get().Steps())
        m_stepsFile = std::ofstream(Configuration::get().baseName() + ".stp");

    m_timeStep = Configuration::get().initStep();
    m_N_ok_tot = 0;
    m_N_fail_tot = 0;
    m_N_bstotal = 0;
    m_N_rkfixtotal = 0;
    m_N_smalltotal = 0;

    m_elapsedTime = 0.0;

    void *dm10, *dm11, *dm12;
    posix_memalign(&dm10, 64, Configuration::get().MaxSubStep() * m_Npad * sizeof(double) * 6);
    posix_memalign(&dm11, 64, 6 * m_Npad * sizeof(double));
    posix_memalign(&dm12, 64, 6 * m_Npad * sizeof(double));
    m_extD = (double*) dm10;
    m_extErr = (double*) dm11;
    m_tmpDat = (double*) dm12;
    m_extH.resize(Configuration::get().MaxSubStep());
}

std::array<double,3> Propagator::com()
{
    std::array<double,3> result{0.0};

    for (int i = 0; i < m_Nobj; i++)
    {
        result[0] += m_x_BSlarge[0 * m_Npad + i] * m_masses[i] / Configuration::get().TotMass();
        result[1] += m_x_BSlarge[1 * m_Npad + i] * m_masses[i] / Configuration::get().TotMass();
        result[2] += m_x_BSlarge[2 * m_Npad + i] * m_masses[i] / Configuration::get().TotMass();
    }

    return result;
}

std::array<double,3> Propagator::linMom()
{
    std::array<double,3> result{0.0};

    for (int i = 0; i < m_Nobj; i++)
    {
        result[0] += m_v_BSlarge[0 * m_Npad + i] * m_masses[i];
        result[1] += m_v_BSlarge[1 * m_Npad + i] * m_masses[i];
        result[2] += m_v_BSlarge[2 * m_Npad + i] * m_masses[i];
    }

    return result;
}

std::array<double,3> Propagator::angMom()
{
    std::array<double,3> result{0.0};

    for (int i = 0; i < m_Nobj; i++)
    {
        result[0] += m_masses[i] * (m_x_BSlarge[1 * m_Npad + i] * m_v_BSlarge[2 * m_Npad + i] - m_x_BSlarge[2 * m_Npad + i] * m_v_BSlarge[1 * m_Npad + i]);
        result[1] += m_masses[i] * (m_x_BSlarge[2 * m_Npad + i] * m_v_BSlarge[0 * m_Npad + i] - m_x_BSlarge[0 * m_Npad + i] * m_v_BSlarge[2 * m_Npad + i]);
        result[2] += m_masses[i] * (m_x_BSlarge[0 * m_Npad + i] * m_v_BSlarge[1 * m_Npad + i] - m_x_BSlarge[1 * m_Npad + i] * m_v_BSlarge[0 * m_Npad + i]);
    }

    return result;
}

double Propagator::radiusOfGyration(double *__restrict__ x)
{
    // calculate geometric center:
    double cx, cy, cz;
    for (int i = 0; i < m_Nobj; i++)
    {
        cx += x[0 * m_Npad + i];
        cy += x[1 * m_Npad + i];
        cz += x[2 * m_Npad + i];
    }
    cx /= m_Nobj;
    cy /= m_Nobj;
    cz /= m_Nobj;

    double gyr = 0.0;
    for (int i = 0; i < m_Nobj; i++)
    {
        double dist = 0.0;
        dist += (x[0 * m_Npad + i] - cx) * (x[0 * m_Npad + i] - cx);
        dist += (x[1 * m_Npad + i] - cy) * (x[1 * m_Npad + i] - cy);
        dist += (x[2 * m_Npad + i] - cz) * (x[2 * m_Npad + i] - cz);

        gyr += std::sqrt(dist);
    }
    gyr /= m_Nobj;

    return gyr;
}

void Propagator::shiftCom()
{
    std::array<double,3> com = this->com();

    for (int i = 0; i < m_Nobj; i++)
    {
        m_x_BSlarge[0 * m_Npad + i] -= com[0];
        m_x_BSlarge[1 * m_Npad + i] -= com[1];
        m_x_BSlarge[2 * m_Npad + i] -= com[2];
    }
}

void Propagator::shiftMom()
{
    std::array<double,3> mom = this->linMom();

    for (int i = 0; i < m_Nobj; i++)
    {
        m_v_BSlarge[0 * m_Npad + i] -= mom[0] / m_totMass;
        m_v_BSlarge[1 * m_Npad + i] -= mom[1] / m_totMass;
        m_v_BSlarge[2 * m_Npad + i] -= mom[2] / m_totMass;
    }
}

void Propagator::propagate()
{
    std::ofstream &out = Configuration::get().outputFile();

    // write the initial conditions:
    out << "----------------------------------\n";
    out << "INITIAL COORDINATES AND VELOCITIES\n";
    out << "----------------------------------\n\n";

    out << "    name      mass (kg)       X (m)      Y (m)      Z (m)     V_x (m/s)  V_y (m/s)  V_z (m/s)\n\n";

    for (int i = 0; i < m_Nobj; i++)
    {
        out << std::setw(10) << m_names[i] << std::setprecision(3) << std::setw(13) << m_masses[i]
            << std::setw(13) << m_x_BSlarge[0 * m_Npad + i]
            << std::setw(11) << m_x_BSlarge[1 * m_Npad + i]
            << std::setw(11) << m_x_BSlarge[2 * m_Npad + i]
            << std::setw(13) << m_v_BSlarge[0 * m_Npad + i]
            << std::setw(11) << m_v_BSlarge[1 * m_Npad + i]
            << std::setw(11) << m_v_BSlarge[2 * m_Npad + i] << "\n";
    }
    out << "\n\n";

    out << "  ----------------------------------------------------------------------------------\n";
    out << "                                  STARTING THE PROPAGATION\n";
    out << "  ----------------------------------------------------------------------------------\n\n";

    if (Configuration::get().intType() == BS)
    {
        out << "                          ***************************************\n";
        out << "                          * USING THE BULIRSCH-STOER INTEGRATOR *\n";
        out << "                          ***************************************\n\n";

        out << "       elapsed time      large steps      BS steps      small steps      cpu time [ms]\n";
        out << "                         good    bad\n";

        if (Configuration::get().Verbose())
            std::cout << " Starting propagation with the Bulirsch-Stoer integrator\n";
    }
    else if (Configuration::get().intType() == Rk4FixM || Configuration::get().intType() == Rk4FixT)
    {
        out << "                ************************************************************\n";
        out << "                * USING THE FIXED-STEP RUNGE-KUTTA INTEGRATOR OF 4TH ORDER *\n";

        if (Configuration::get().intType() == Rk4FixM)
            out << "                *                      WITH MIDPOINT STEP                  *\n";
        else
            out << "                *                       WITH THIRDS STEP                   *\n";

        out << "                ************************************************************\n\n";

        out << "       elapsed time            RK steps               cpu time [ms]\n\n";

        if (Configuration::get().Verbose())
            std::cout << " Starting propagation with the fixed-step Runge-Kutta integrator of fourth order\n";
    }

    while (true)
    {
        writeToTrj();

        if (Configuration::get().Verbose())
            writeStatus();

        if (Configuration::get().Restart())
        {
            std::ofstream resFile(Configuration::get().baseName() + ".rst");
            resFile << m_elapsedTime << "\n";

            for (int i = 0; i < m_Nobj; i++)
                resFile << m_masses[i]
                        << m_x_BSlarge[0 * m_Npad + i]
                        << m_x_BSlarge[1 * m_Npad + i]
                        << m_x_BSlarge[2 * m_Npad + i]
                        << m_v_BSlarge[0 * m_Npad + i]
                        << m_v_BSlarge[1 * m_Npad + i]
                        << m_v_BSlarge[2 * m_Npad + i] << "\n";

            resFile.close();
        }

        if (m_elapsedTime >= Configuration::get().tfinal())
            break;

        if (Configuration::get().Steps())
        {
            m_stepsFile << "# elapsed time:" << std::scientific << std::setprecision(9) << std::setw(20) << std::uppercase << m_elapsedTime << "\n";
            m_stepsFile << "# trying propagation with step:" << std::setw(20) << m_timeStep << "\n";
        }

        clock_t t0 = clock();
        switch (Configuration::get().intType()) {
        case BS:
            BS_LargeStep();
            break;
        case Rk4FixM:
        case Rk4FixT:
            RK4F_LargeStep();
            break;
        default:
            break;
        }
        clock_t t1 = clock();
        double runTime = ((double) (t1 - t0)) / CLOCKS_PER_SEC;

        if (Configuration::get().Steps())
        {
            m_stepsFile << "# successful / failed steps:" << m_N_ok << m_N_fail << "\n";
            m_stepsFile << "#\n";
            m_stepsFile << "######################################################################\n";
            m_stepsFile << "#\n";
        }

        switch (Configuration::get().intType()) {
        case BS:
            out << "       " << std::setprecision(5) << std::setw(11) << m_elapsedTime
                << "      " << std::setw(5) << m_N_ok << " " << std::setw(5) << m_N_fail << " " << std::setw(5) << m_N_bssteps << " " << std::setw(7) << m_N_smallsteps << " " << runTime * 1000.0 << "\n";
            break;
        case Rk4FixM:
        case Rk4FixT:
            out << m_elapsedTime << runTime * 1000.0 << "\n";
            break;
        default:
            break;
        }
    }

    out << "\n";
    out << "       **************************************************************\n";
    out << "       *                 SUMMARY OF THE CALCULATION                 *\n";
    out << "       *                                                            *\n";

    switch (Configuration::get().intType()) {
    case BS:
        out << "       * total time [s]     large steps      BS steps      small steps *\n";
        out << "       *                    good    bad                                *\n";
        out << m_elapsedTime << m_N_ok_tot << m_N_fail_tot << m_N_bstotal << m_N_smalltotal << "\n";
        break;
    case Rk4FixM:
    case Rk4FixT:
        out << "       * total time [s]            RK4 steps                        *\n";
        out << m_elapsedTime << m_N_rkfixtotal << "\n";
        break;
    default:
        break;
    }

    out << "       **************************************************************\n\n";

    if (Configuration::get().Steps())
        m_stepsFile << "# Simulation done!\n";
}

void Propagator::writeToTrj()
{
    if (Configuration::get().TextTrj())
    {
        m_txtTrj << std::uppercase << std::setprecision(10) << std::scientific << std::setw(18) << m_elapsedTime;

        int prec = Configuration::get().Ndigit();
        int width = prec + 8;
        for (int i = 0; i < m_Nobj; i++)
        {
            m_txtTrj << " "
                     << std::setprecision(prec) << std::setw(width) << m_x_BSlarge[0 * m_Npad + i]
                     << std::setprecision(prec) << std::setw(width) << m_x_BSlarge[1 * m_Npad + i]
                     << std::setprecision(prec) << std::setw(width) << m_x_BSlarge[2 * m_Npad + i];
        }
        m_txtTrj << "\n";
    }

    m_binTrj.write((char*) &m_elapsedTime, sizeof(double));

    for (int i = 0; i < m_Nobj; i++)
    {
        m_binTrj.write((char*) &m_x_BSlarge[0 * m_Npad + i], sizeof(double));
        m_binTrj.write((char*) &m_x_BSlarge[1 * m_Npad + i], sizeof(double));
        m_binTrj.write((char*) &m_x_BSlarge[2 * m_Npad + i], sizeof(double));
    }
    for (int i = 0; i < m_Nobj; i++)
    {
        m_binTrj.write((char*) &m_v_BSlarge[0 * m_Npad + i], sizeof(double));
        m_binTrj.write((char*) &m_v_BSlarge[1 * m_Npad + i], sizeof(double));
        m_binTrj.write((char*) &m_v_BSlarge[2 * m_Npad + i], sizeof(double));
    }
}

void Propagator::writeStatus()
{
    printf("\r Simulated %9.3e out of %9.3e seconds.", m_elapsedTime, Configuration::get().tfinal());
    fflush(stdout);
}

void Propagator::BS_LargeStep()
{
    Configuration::get().inc_BSL();

    double internalElapsedTime = 0.0;

    memcpy(m_x_BSLtmp, m_x_BSlarge, 3 * m_Npad * sizeof(double));
    memcpy(m_v_BSLtmp, m_v_BSlarge, 3 * m_Npad * sizeof(double));

    m_N_ok = 0;
    m_N_fail = 0;
    m_N_bstotal = 0;
    m_N_smalltotal = 0;

    while (true)
    {
        // for a normal propagation, the stepsize must not overshoot, so it
        // will be reduced (adding the minimum step to avoid numerical instabilities):
        if (! Configuration::get().UnRes())
        {
            if (internalElapsedTime + m_timeStep + Configuration::get().minStep() > Configuration::get().tout())
                m_timeStep = Configuration::get().tout() - internalElapsedTime;
        }

        bool success = BS_OneStep();

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
        m_elapsedTime += m_doneStep;

        if (internalElapsedTime >= Configuration::get().tout())
        {
            break;
        }
    }

    memcpy(m_x_BSlarge, m_x_SubFin, 3 * m_Npad * sizeof(double));
    memcpy(m_v_BSlarge, m_v_SubFin, 3 * m_Npad * sizeof(double));
}

int Propagator::BS_OneStep()
{
    Configuration::get().inc_BSO();

    acceleration(m_x_BSLtmp, m_a_BSstart);
    double gyr = radiusOfGyration(m_x_BSLtmp);
    double v_avg = 0.0;
    for (int i = 0; i < m_Nobj; i++)
    {
        v_avg += std::abs(m_v_BSLtmp[0 * m_Npad + i]) + std::abs(m_v_BSLtmp[1 * m_Npad + i]) + std::abs(m_v_BSLtmp[2 * m_Npad + i]);
    }
    v_avg /= (3 * m_Nobj);

    double trialStep = m_timeStep;

    if (Configuration::get().Steps())
    {
        m_stepsFile << "#          Time            Stepsize      Substeps       Error\n";
    }

    bool success = true;
    while (true)
    {
        m_N_bssteps++;

        m_delta = 0.0;
        for (int i = 1; i <= Configuration::get().MaxSubStep(); i++)
        {
            BS_SubSteps(i, trialStep);

            m_N_smallsteps += i;
            double h_est = (trialStep / i) * (trialStep / i);

            BS_Extrapolate(i, h_est);

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
                m_stepsFile << "  " << std::setprecision(8) << std::setw(18) << m_elapsedTime
                            << std::setprecision(5) << std::setw(17) << trialStep
                            << std::setw(9) << i << std::setw(18) << m_delta << "\n";
            }

            if (m_delta < Configuration::get().Eps())
            {
                memcpy(m_x_BSLtmp, m_x_SubFin, 3 * m_Npad * sizeof(double));
                memcpy(m_v_BSLtmp, m_v_SubFin, 3 * m_Npad * sizeof(double));
                m_doneStep = trialStep;
                goto finished;
            }
        }

        success = false;

        if (trialStep <= Configuration::get().minStep())
        {
            std::cout << "WARNING: No convergence with minimum stepsize. Aborting!\n";
            m_underflow = true;
            goto finished;
        }

        double factor = std::pow(Configuration::get().Eps() / m_delta, 1.0 / Configuration::get().MaxSubStep());
        factor = std::max(std::min(factor, Configuration::get().RedMin()), Configuration::get().RedMax());

        if (trialStep * factor < Configuration::get().minStep())
        {
            trialStep = Configuration::get().minStep();
        }
        else
        {
            trialStep *= factor;
        }
    }

    finished:

    return success;
}

void Propagator::BS_SubSteps(const int nsteps, const double stepSize)
{
    Configuration::get().inc_BSS();

    double step = stepSize / nsteps;

    /*
     * the first step
     */
    for (int i = 0; i < m_Nobj; i++)
    {
        m_x_SubStep[0 * m_Npad + i] = step * (m_v_BSLtmp[0 * m_Npad + i] + 0.5 * step * m_a_BSstart[0 * m_Npad + i]);
        m_x_SubStep[1 * m_Npad + i] = step * (m_v_BSLtmp[1 * m_Npad + i] + 0.5 * step * m_a_BSstart[1 * m_Npad + i]);
        m_x_SubStep[2 * m_Npad + i] = step * (m_v_BSLtmp[2 * m_Npad + i] + 0.5 * step * m_a_BSstart[2 * m_Npad + i]);
        m_x_SubFin[0 * m_Npad + i] = m_x_BSLtmp[0 * m_Npad + i] + m_x_SubStep[0 * m_Npad + i];
        m_x_SubFin[1 * m_Npad + i] = m_x_BSLtmp[1 * m_Npad + i] + m_x_SubStep[1 * m_Npad + i];
        m_x_SubFin[2 * m_Npad + i] = m_x_BSLtmp[2 * m_Npad + i] + m_x_SubStep[2 * m_Npad + i];
        m_a_SubInt[0 * m_Npad + i] = m_a_BSstart[0 * m_Npad + i];
        m_a_SubInt[1 * m_Npad + i] = m_a_BSstart[1 * m_Npad + i];
        m_a_SubInt[2 * m_Npad + i] = m_a_BSstart[2 * m_Npad + i];
    }

    /*
     * the remaining steps
     */
    for (int i = 2; i <= nsteps; i++)
    {
        acceleration(m_x_SubFin, m_a_SubInt);

        /*
         * accumulate the step to reduce roundoff error according to p. 726 Numerical recipes in Fortran 77
         */
        for (int j = 0; j < m_Nobj; j++)
        {
            m_x_SubStep[0 * m_Npad + j] += step * step * m_a_SubInt[0 * m_Npad + j];
            m_x_SubStep[1 * m_Npad + j] += step * step * m_a_SubInt[1 * m_Npad + j];
            m_x_SubStep[2 * m_Npad + j] += step * step * m_a_SubInt[2 * m_Npad + j];
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
        m_v_SubFin[0 * m_Npad + i] = (m_x_SubStep[0 * m_Npad + i] / step) + 0.5 * step * m_a_SubInt[0 * m_Npad + i];
        m_v_SubFin[1 * m_Npad + i] = (m_x_SubStep[1 * m_Npad + i] / step) + 0.5 * step * m_a_SubInt[1 * m_Npad + i];
        m_v_SubFin[2 * m_Npad + i] = (m_x_SubStep[2 * m_Npad + i] / step) + 0.5 * step * m_a_SubInt[2 * m_Npad + i];
    }
}

void Propagator::BS_Extrapolate(const int i_est, const double h_est)
{
    Configuration::get().inc_BSE();

    m_extH[i_est - 1] = h_est;

    for (int i = 0; i < 3 * m_Npad; i++)
    {
        m_extErr[0 * m_Npad + i] = m_x_SubFin[i];
        m_extErr[3 * m_Npad + i] = m_v_SubFin[i];
        m_tmpDat[0 * m_Npad + i] = m_x_SubFin[i];
        m_tmpDat[3 * m_Npad + i] = m_v_SubFin[i];
    }

    double c[6 * m_Npad];

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
            c[0 * m_Npad + i] = m_x_SubFin[i];
            c[3 * m_Npad + i] = m_v_SubFin[i];
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
                delta = c[j] - tmp3;
                m_extErr[j] = tmp1 * delta;
                c[j] = tmp2 * delta;
                m_tmpDat[j] += m_extErr[j];
            }
        }

        for (int i = 0; i < 6 * m_Npad; i++)
        {
            m_extD[(i_est - 1) * 6 * m_Npad + i] = m_extErr[i];
        }
    }

    memcpy(m_x_SubFin, m_tmpDat + m_Npad * 0, 3 * m_Npad * sizeof(double));
    memcpy(m_v_SubFin, m_tmpDat + m_Npad * 3, 3 * m_Npad * sizeof(double));
}

void Propagator::RK4F_LargeStep()
{
}

void Propagator::acceleration(double *__restrict__ x, double *__restrict__ a)
{
    int i, j;

    for (j = 0; j < m_Nobj; j++)
    {
        a[0 * m_Npad + j] = 0.0;
        a[1 * m_Npad + j] = 0.0;
        a[2 * m_Npad + j] = 0.0;
    }

    for (i = 0; i < m_Nobj - 1; i++)
    {
        for (j = i + 1; j < m_Nobj; j++)
        {
            double dX = x[0 * m_Npad + j] - x[0 * m_Npad + i];
            double dY = x[1 * m_Npad + j] - x[1 * m_Npad + i];
            double dZ = x[2 * m_Npad + j] - x[2 * m_Npad + i];
            double R2 = dX * dX + dY * dY + dZ * dZ;
            double tmpFac = PhyCon::G / (R2 * sqrt(R2));
            a[0 * m_Npad + i] += m_masses[j] * tmpFac * dX;
            a[1 * m_Npad + i] += m_masses[j] * tmpFac * dY;
            a[2 * m_Npad + i] += m_masses[j] * tmpFac * dZ;
            a[0 * m_Npad + j] -= m_masses[i] * tmpFac * dX;
            a[1 * m_Npad + j] -= m_masses[i] * tmpFac * dY;
            a[2 * m_Npad + j] -= m_masses[i] * tmpFac * dZ;
        }
    }
}

}
