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


#include"object.h"

namespace guistronx
{

class ObjectPrivate
{
public:
    ObjectPrivate()
    {
    }

    ObjectPrivate(const std::string &init_name, const double &init_mass, const Eigen::Vector3d &init_pos, const Eigen::Vector3d &init_vel)
        : name(init_name)
        , mass(init_mass)
        , position(init_pos)
        , velocity(init_vel)
        , acceleration(Eigen::Vector3d::Zero())
    {
    }

    std::string name;
    double mass;
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    Eigen::Vector3d acceleration;
};


Object::Object()
    : d(new ObjectPrivate())
{
}


Object::Object(const std::string &init_name, const double &init_mass, const Eigen::Vector3d &init_pos, const Eigen::Vector3d &init_vel)
    : d(new ObjectPrivate(init_name, init_mass, init_pos, init_vel))
{
}


Object::Object(const Object &oldObject)
    : d(new ObjectPrivate(oldObject.name(), oldObject.mass(), oldObject.position(), oldObject.velocity()))
{
}


std::string Object::name() const
{
    return d->name;
}

double Object::mass() const
{
    return d->mass;
}


Eigen::Vector3d Object::position() const
{
    return d->position;
}


Eigen::Vector3d Object::velocity() const
{
    return d->velocity;
}


Eigen::Vector3d Object::acceleration() const
{
    return d->acceleration;
}

}  // namespace guistronx
