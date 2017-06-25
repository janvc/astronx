/*
 * Copyright 2015 Jan von Cosel
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
 * along with molconv. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <vector>
#include "system.h"

namespace guistronx
{

class SystemPrivate
{
public:
    std::vector<Object*> m_objects;
};

System::System()
    : d(new SystemPrivate)
{
}


size_t System::N_obj() const
{
    return d->m_objects.size();
}

Object *System::get_object(const size_t index) const
{
    return d->m_objects.at(index);
}

double System::mass() const
{
    return calc_mass();
}

Eigen::Vector3d System::centerOfMass() const
{
    return calc_com();
}

Eigen::Vector3d System::centerOfGeometry() const
{
    return calc_cog();
}

Eigen::Vector3d System::linearMomentum() const
{
    return calc_linmom();
}

Eigen::Vector3d System::angularMomentum() const
{
    return calc_angmom();
}

void System::addObject(const Object &newObject)
{
    d->m_objects.push_back(new Object(newObject));
}

Eigen::Vector3d System::calc_com() const
{
}

Eigen::Vector3d System::calc_cog() const
{
}

Eigen::Vector3d System::calc_linmom() const
{
}

Eigen::Vector3d System::calc_angmom() const
{
}

double System::calc_mass() const
{
}

double System::calc_gyrate() const
{
}

}
