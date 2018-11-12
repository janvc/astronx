#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "propagator.h"
#include "configuration.h"

namespace Astronx
{

Propagator::Propagator()
{
    m_Nobj = Configuration::get().Nobj();
    m_totMass = Configuration::get().TotMass();

    // make sure the array length is divisible by 4 (AVX register length)
    m_Npad = m_Nobj % 4 == 0 ? m_Nobj : ((m_Nobj / 4) + 1) * 4;

    void *dm0, *dm1, *dm2, *dm3, *dm4, *dm5, *dm6;
    posix_memalign(&dm0, 64, m_Npad * sizeof(double));
    posix_memalign(&dm1, 64, m_Npad * sizeof(double));
    posix_memalign(&dm2, 64, m_Npad * sizeof(double));
    posix_memalign(&dm3, 64, m_Npad * sizeof(double));
    posix_memalign(&dm4, 64, m_Npad * sizeof(double));
    posix_memalign(&dm5, 64, m_Npad * sizeof(double));
    posix_memalign(&dm6, 64, m_Npad * sizeof(double));

    m_xx = (double*) dm0;
    m_xy = (double*) dm1;
    m_xz = (double*) dm2;
    m_vx = (double*) dm3;
    m_vy = (double*) dm4;
    m_vz = (double*) dm5;
    m_masses = (double*) dm6;

    for (int i = 0; i < m_Nobj; i++)
    {
        m_xx[i] = Configuration::get().XX0()[i];
        m_xy[i] = Configuration::get().XY0()[i];
        m_xz[i] = Configuration::get().XZ0()[i];
        m_vx[i] = Configuration::get().VX0()[i];
        m_vy[i] = Configuration::get().VY0()[i];
        m_vz[i] = Configuration::get().VZ0()[i];
        m_masses[i] = Configuration::get().masses()[i];
        m_names = Configuration::get().names();
    }

    m_binTrj = std::ofstream(Configuration::get().baseName() + ".bin.trj", std::ios::binary);

    if (Configuration::get().TextTrj())
        m_txtTrj = std::ofstream(Configuration::get().baseName() + ".txt.trj");

    if (Configuration::get().Restart())
        m_restartFile = std::ofstream(Configuration::get().baseName() + ".rst");

    if (Configuration::get().Steps())
        m_stepsFile = std::ofstream(Configuration::get().baseName() + ".stp");
}

std::array<double,3> Propagator::com()
{
    std::array<double,3> result{0.0};

    for (int i = 0; i < m_Nobj; i++)
    {
        result[0] += m_xx[i] * m_masses[i] / Configuration::get().TotMass();
        result[1] += m_xy[i] * m_masses[i] / Configuration::get().TotMass();
        result[2] += m_xz[i] * m_masses[i] / Configuration::get().TotMass();
    }

    return result;
}

std::array<double,3> Propagator::linMom()
{
    std::array<double,3> result{0.0};

    for (int i = 0; i < m_Nobj; i++)
    {
        result[0] += m_vx[i] * m_masses[i];
        result[1] += m_vy[i] * m_masses[i];
        result[2] += m_vz[i] * m_masses[i];
    }

    return result;
}

std::array<double,3> Propagator::angMom()
{
    std::array<double,3> result{0.0};

    for (int i = 0; i < m_Nobj; i++)
    {
        result[0] += m_masses[i] * (m_xy[i] * m_vz[i] - m_xz[i] * m_vy[i]);
        result[1] += m_masses[i] * (m_xz[i] * m_vx[i] - m_xx[i] * m_vz[i]);
        result[2] += m_masses[i] * (m_xx[i] * m_vy[i] - m_xy[i] * m_vx[i]);
    }

    return result;
}

void Propagator::shiftCom()
{
    std::array<double,3> com = this->com();

    for (int i = 0; i < m_Nobj; i++)
    {
        m_xx[i] -= com[0];
        m_xy[i] -= com[1];
        m_xz[i] -= com[2];
    }
}

void Propagator::shiftMom()
{
    std::array<double,3> mom = this->linMom();

    for (int i = 0; i < m_Nobj; i++)
    {
        m_vx[i] -= mom[0] / m_totMass;
        m_vy[i] -= mom[1] / m_totMass;
        m_vz[i] -= mom[2] / m_totMass;
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
            << std::setw(13) << m_xx[i] << std::setw(11) << m_xy[i] << std::setw(11) << m_xz[i]
            << std::setw(13) << m_vx[i] << std::setw(11) << m_vy[i] << std::setw(11) << m_vz[i] << "\n";
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


    switch (Configuration::get().intType()) {
    case BS:
        propBS();
        break;
    default:
        break;
    }
}

void Propagator::propBS()
{
    std::ofstream &out = Configuration::get().outputFile();


}

void Propagator::writeToTrj()
{
    if (Configuration::get().TextTrj())
    {
        m_txtTrj << std::setprecision(10) << std::scientific << std::setw(18) << m_elapsedTime;

        for (int i = 0; i < m_Nobj; i++)
        {
            m_txtTrj << std::setw(18) << m_xx[i]
                     << std::setw(18) << m_xy[i]
                     << std::setw(18) << m_xz[i];
        }
        m_txtTrj << "\n";
    }

    m_binTrj.write((char*) &m_elapsedTime, sizeof(double));

    for (int i = 0; i < m_Nobj; i++)
    {
        m_binTrj.write((char*) &m_xx[i], sizeof(double));
        m_binTrj.write((char*) &m_xy[i], sizeof(double));
        m_binTrj.write((char*) &m_xz[i], sizeof(double));
    }
    for (int i = 0; i < m_Nobj; i++)
    {
        m_binTrj.write((char*) &m_vx[i], sizeof(double));
        m_binTrj.write((char*) &m_vy[i], sizeof(double));
        m_binTrj.write((char*) &m_vz[i], sizeof(double));
    }
}

}
