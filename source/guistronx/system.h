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

#ifndef SYSTEM_H
#define SYSTEM_H

#include"object.h"

namespace guistronx
{

class SystemPrivate;

class System
{
public:
    System();
    size_t N_obj() const;
    Object *get_object(const size_t index) const;
    double mass() const;
    double radiusOfGyration() const;
    Eigen::Vector3d centerOfMass() const;
    Eigen::Vector3d centerOfGeometry() const;
    Eigen::Vector3d linearMomentum() const;
    Eigen::Vector3d angularMomentum() const;

    void addObject(const Object &newObject);

private:
    Eigen::Vector3d calc_com() const;
    Eigen::Vector3d calc_cog() const;
    Eigen::Vector3d calc_linmom() const;
    Eigen::Vector3d calc_angmom() const;
    double calc_mass() const;
    double calc_gyrate() const;

    SystemPrivate *d;
};

}  // namespace guistronx

#endif // SYSTEM_H
