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


#ifndef SYSTEM_H
#define SYSTEM_H

#include <array>
#include <vector>
#include "propagator.h"


namespace Astronx
{

class System
{
public:
    System();
    ~System();

    double totMass() const;

    std::array<double,3> com() const;
    std::array<double,3> linMom() const;
    std::array<double,3> angMom() const;
    void shiftCom();
    void shiftMom();

    void propagate();

private:
    void writeToTrj();      // write current configuration to file
    void writeStatus();     // write current status to stdout

    int m_Nobj;             // number of objects
    int m_Npad;             // padded array length, divisible by 4

    double m_totMass;       // total mass of the system
    double m_elapsedTime;   // the time

    double *m_masses;   // masses of objects
    double *m_xLarge;   // positions and velocities
    double *m_vLarge;   // between large propagation steps

    std::vector<std::string> m_names;   // object names
};

}

#endif // SYSTEM_H
