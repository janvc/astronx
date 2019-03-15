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


#include "system.h"
#include "configuration.h"
#include "propagator.h"
#include "bulirschstoer.h"


namespace Astronx
{

System::System()
{
    m_Nobj = Configuration::get().Nobj();

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
}

System::~System()
{
}

double System::totMass() const
{
    return m_totMass;
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

void System::propagate()
{
    Propagator *prop;

    switch (Configuration::get().intType()) {
    case BS:
        prop = new BulirschStoer();
        break;
    default:
        break;
    }

    prop->largeStep();
}

}
