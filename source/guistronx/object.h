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

#ifndef OBJECT_H
#define OBJECT_H

#include<eigen3/Eigen/Core>

namespace guistronx
{

class ObjectPrivate;

class Object
{
public:
    Object();
    Object(const std::string &init_name, const double &init_mass, const Eigen::Vector3d &init_pos, const Eigen::Vector3d &init_vel);
    Object(const Object &oldObject);

    std::string name() const;
    double mass() const;
    Eigen::Vector3d position() const;
    Eigen::Vector3d velocity() const;
    Eigen::Vector3d acceleration() const;

private:
    ObjectPrivate *d;
};

}  // namespace guistronx

#endif // OBJECT_H
