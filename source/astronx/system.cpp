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
    m_x = (double*) dm0;
    m_v = (double*) dm1;
    m_masses = (double*) dm2;

    for (int i = 0; i < m_Nobj; i++)
    {
        m_x[0 * m_Npad + i] = Configuration::get().XX0()[i];
        m_x[1 * m_Npad + i] = Configuration::get().XY0()[i];
        m_x[2 * m_Npad + i] = Configuration::get().XZ0()[i];
        m_v[0 * m_Npad + i] = Configuration::get().VX0()[i];
        m_v[1 * m_Npad + i] = Configuration::get().VY0()[i];
        m_v[2 * m_Npad + i] = Configuration::get().VZ0()[i];
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

}
