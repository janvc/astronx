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
#include <sys/time.h>
#include "system.h"
#include "configuration.h"
#include "propagator.h"
#include "bulirschstoer.h"
#include "stepsizeunderflow.h"


namespace Astronx
{

System::System()
{
    m_Nobj = Configuration::get().Nobj();
    m_elapsedTime = 0.0;

    // make sure the array length is divisible by 4 (AVX register length)
    m_Npad = m_Nobj % 4 == 0 ? m_Nobj : ((m_Nobj / 4) + 1) * 4;

    void *dm0, *dm1, *dm2;
    posix_memalign(&dm0, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm1, 64, 3 * m_Npad * sizeof(double));
    posix_memalign(&dm2, 64, m_Npad * sizeof(double));
    m_xLarge = (double*) dm0;
    m_vLarge = (double*) dm1;
    m_masses = (double*) dm2;

    for (int i = 0; i < m_Nobj; i++)
    {
        m_xLarge[0 * m_Npad + i] = Configuration::get().XX0()[i];
        m_xLarge[1 * m_Npad + i] = Configuration::get().XY0()[i];
        m_xLarge[2 * m_Npad + i] = Configuration::get().XZ0()[i];
        m_vLarge[0 * m_Npad + i] = Configuration::get().VX0()[i];
        m_vLarge[1 * m_Npad + i] = Configuration::get().VY0()[i];
        m_vLarge[2 * m_Npad + i] = Configuration::get().VZ0()[i];
        m_masses[i] = Configuration::get().masses()[i];
        m_names = Configuration::get().names();
    }

    // calculate the total mass
    m_totMass = 0.0;
    for (int i = 0; i < m_Nobj; i++)
    {
        m_totMass += m_masses[i];
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
        m_txtTrj.flush();
    }
}

System::~System()
{
}

double System::totMass() const
{
    return m_totMass;
}

double System::elapsedTime() const
{
    return m_elapsedTime;
}

std::array<double,3> System::com() const
{
    std::array<double,3> result{0.0};

    for (int i = 0; i < m_Nobj; i++)
    {
        result[0] += m_xLarge[0 * m_Npad + i] * m_masses[i] / m_totMass;
        result[1] += m_xLarge[1 * m_Npad + i] * m_masses[i] / m_totMass;
        result[2] += m_xLarge[2 * m_Npad + i] * m_masses[i] / m_totMass;
    }

    return result;
}

std::array<double,3> System::linMom() const
{
    std::array<double,3> result{0.0};

    for (int i = 0; i < m_Nobj; i++)
    {
        result[0] += m_vLarge[0 * m_Npad + i] * m_masses[i];
        result[1] += m_vLarge[1 * m_Npad + i] * m_masses[i];
        result[2] += m_vLarge[2 * m_Npad + i] * m_masses[i];
    }

    return result;
}

std::array<double,3> System::angMom() const
{
    std::array<double,3> result{0.0};

    for (int i = 0; i < m_Nobj; i++)
    {
        result[0] += m_masses[i] * (m_xLarge[1 * m_Npad + i] * m_vLarge[2 * m_Npad + i] - m_xLarge[2 * m_Npad + i] * m_vLarge[1 * m_Npad + i]);
        result[1] += m_masses[i] * (m_xLarge[2 * m_Npad + i] * m_vLarge[0 * m_Npad + i] - m_xLarge[0 * m_Npad + i] * m_vLarge[2 * m_Npad + i]);
        result[2] += m_masses[i] * (m_xLarge[0 * m_Npad + i] * m_vLarge[1 * m_Npad + i] - m_xLarge[1 * m_Npad + i] * m_vLarge[0 * m_Npad + i]);
    }

    return result;
}

void System::shiftCom()
{
    std::array<double,3> com = this->com();

    for (int i = 0; i < m_Nobj; i++)
    {
        m_xLarge[0 * m_Npad + i] -= com[0];
        m_xLarge[1 * m_Npad + i] -= com[1];
        m_xLarge[2 * m_Npad + i] -= com[2];
    }
}

void System::shiftMom()
{
    std::array<double,3> mom = this->linMom();

    for (int i = 0; i < m_Nobj; i++)
    {
        m_vLarge[0 * m_Npad + i] -= mom[0] / m_totMass;
        m_vLarge[1 * m_Npad + i] -= mom[1] / m_totMass;
        m_vLarge[2 * m_Npad + i] -= mom[2] / m_totMass;
    }
}

void System::writeToTrj()
{
    if (Configuration::get().TextTrj())
    {
        m_txtTrj << std::uppercase << std::setprecision(10) << std::scientific << std::setw(18) << m_elapsedTime;

        int prec = Configuration::get().Ndigit();
        int width = prec + 8;
        for (int i = 0; i < m_Nobj; i++)
        {
            m_txtTrj << " "
                     << std::setprecision(prec) << std::setw(width) << m_xLarge[0 * m_Npad + i]
                     << std::setprecision(prec) << std::setw(width) << m_xLarge[1 * m_Npad + i]
                     << std::setprecision(prec) << std::setw(width) << m_xLarge[2 * m_Npad + i];
        }
        m_txtTrj << "\n";
    }

    std::ofstream &binTrj = Configuration::get().binTrjFile();
    binTrj.write((char*) &m_elapsedTime, sizeof(double));

    for (int i = 0; i < m_Nobj; i++)
    {
        binTrj.write((char*) &m_xLarge[0 * m_Npad + i], sizeof(double));
        binTrj.write((char*) &m_xLarge[1 * m_Npad + i], sizeof(double));
        binTrj.write((char*) &m_xLarge[2 * m_Npad + i], sizeof(double));
    }
    for (int i = 0; i < m_Nobj; i++)
    {
        binTrj.write((char*) &m_vLarge[0 * m_Npad + i], sizeof(double));
        binTrj.write((char*) &m_vLarge[1 * m_Npad + i], sizeof(double));
        binTrj.write((char*) &m_vLarge[2 * m_Npad + i], sizeof(double));
    }
}

void System::writeStatus()
{
    printf("\r Simulated %9.3e out of %9.3e seconds.", m_elapsedTime, Configuration::get().tfinal());
    fflush(stdout);
}

void System::writeRestart()
{
    std::ofstream resFile(Configuration::get().baseName() + ".rst");
    resFile << m_elapsedTime << "\n";

    for (int i = 0; i < m_Nobj; i++)
        resFile << m_masses[i]
                << m_xLarge[0 * m_Npad + i]
                << m_xLarge[1 * m_Npad + i]
                << m_xLarge[2 * m_Npad + i]
                << m_vLarge[0 * m_Npad + i]
                << m_vLarge[1 * m_Npad + i]
                << m_vLarge[2 * m_Npad + i] << "\n";

    resFile.close();
}

void System::getWalltime(double *wcTime)
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    *wcTime = (double) (tp.tv_sec + tp.tv_usec / 1000000.0);
}

void System::propagate()
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
            << std::setw(13) << m_xLarge[0 * m_Npad + i]
            << std::setw(11) << m_xLarge[1 * m_Npad + i]
            << std::setw(11) << m_xLarge[2 * m_Npad + i]
            << std::setw(13) << m_vLarge[0 * m_Npad + i]
            << std::setw(11) << m_vLarge[1 * m_Npad + i]
            << std::setw(11) << m_vLarge[2 * m_Npad + i] << "\n";
    }
    out << "\n\n";

    out << "  ----------------------------------------------------------------------------------\n";
    out << "                                  STARTING THE PROPAGATION\n";
    out << "  ----------------------------------------------------------------------------------\n\n";

    out.flush();

    double totalBegin;
    getWalltime(&totalBegin);
    Propagator *prop;

    switch (Configuration::get().intType()) {
    case BS:
        prop = new BulirschStoer(m_Npad, this);
        break;
    default:
        break;
    }

    try
    {
        while (true)
        {
            writeToTrj();

            if (Configuration::get().Verbose())
            {
                writeStatus();
            }

            if (Configuration::get().Restart())
            {
                writeRestart();
            }

            if (m_elapsedTime >= Configuration::get().tfinal())
            {
                break;
            }

            if (Configuration::get().Steps())
            {
                std::ofstream &stepsFile = Configuration::get().stepsFile();
                stepsFile << "# elapsed time:" << std::scientific << std::setprecision(9) << std::setw(20) << std::uppercase << m_elapsedTime << "\n";
            }

            double t0, t1;
            getWalltime(&t0);
            m_elapsedTime += prop->largeStep(m_xLarge, m_vLarge);
            getWalltime(&t1);
            double cpuTimeUsed = t1 - t0;

            prop->writeOutputLine(cpuTimeUsed);
        }
    }
    catch (StepsizeUnderflow &s)
    {
        std::cout << std::endl << s.what() << std::endl;

        out << std::endl;
        out << "              **************************************************\n";
        out << "              * WARNING: No convergence with minimum stepsize. *\n";
        out << "              *            Aborting the simulation.            *\n";
        out << "              **************************************************\n";
    }

    double totalEnd;
    getWalltime(&totalEnd);

    prop->writeSummary();

    out << "  ----------------------------------------------------------------------------------\n";
    out << "                                 FINISHED THE PROPAGATION\n";
    out << "                           total cpu time:"
        << std::fixed << std::setprecision(3) << std::setw(12) << totalEnd - totalBegin << " seconds\n";
    out << "  ----------------------------------------------------------------------------------\n";

    out << "--------------------------------\n";
    out << "FINAL COORDINATES AND VELOCITIES\n";
    out << "--------------------------------\n\n";

    out << "    name      mass (kg)       X (m)      Y (m)      Z (m)     V_x (m/s)  V_y (m/s)  V_z (m/s)\n\n";

    for (int i = 0; i < m_Nobj; i++)
    {
        out << std::setw(10) << m_names[i] << std::scientific << std::setprecision(3) << std::setw(13) << m_masses[i]
            << std::setw(13) << m_xLarge[0 * m_Npad + i]
            << std::setw(11) << m_xLarge[1 * m_Npad + i]
            << std::setw(11) << m_xLarge[2 * m_Npad + i]
            << std::setw(13) << m_vLarge[0 * m_Npad + i]
            << std::setw(11) << m_vLarge[1 * m_Npad + i]
            << std::setw(11) << m_vLarge[2 * m_Npad + i] << "\n";
    }
    out << "\n\n\n";
}

}
