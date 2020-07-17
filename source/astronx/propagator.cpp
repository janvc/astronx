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


#include <iomanip>
#include <cmath>
#include "configuration.h"
#include "propagator.h"


namespace Astronx
{

Propagator::Propagator(const int Npad, System* sys)
{
    m_sys = sys;
    m_Nobj = Configuration::get().Nobj();
    m_Npad = Npad;
    m_timeStep = Configuration::get().initStep();

    void *dm0;
    posix_memalign(&dm0, 64, m_Npad * sizeof(double));
    m_masses = (double*) dm0;

    for (int i = 0; i < m_Nobj; i++)
    {
        m_masses[i] = Configuration::get().masses()[i];
        m_names = Configuration::get().names();
    }

    m_binTrj = std::ofstream(Configuration::get().baseName() + ".bin.trj", std::ios::binary);
    for (int i = 0;i < m_Nobj; i++)
    {
        m_binTrj.write((char*) &m_masses[i], sizeof (double));
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
}

Propagator::~Propagator()
{
}

double Propagator::largeStep(double *x, double *v)
{
}

void Propagator::writeOutputLine(const double cpuTimeUsed)
{
}

void Propagator::writeSummary()
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
            double tmpFac = PhyCon::G / (R2 * std::sqrt(R2));
            a[0 * m_Npad + i] += m_masses[j] * tmpFac * dX;
            a[1 * m_Npad + i] += m_masses[j] * tmpFac * dY;
            a[2 * m_Npad + i] += m_masses[j] * tmpFac * dZ;
            a[0 * m_Npad + j] -= m_masses[i] * tmpFac * dX;
            a[1 * m_Npad + j] -= m_masses[i] * tmpFac * dY;
            a[2 * m_Npad + j] -= m_masses[i] * tmpFac * dZ;
        }
    }
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

}
