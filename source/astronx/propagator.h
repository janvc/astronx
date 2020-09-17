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


#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include "system.h"

namespace Astronx
{

class Propagator
{
public:
    Propagator(System *sys);
    virtual ~Propagator();
    virtual double largeStep(double *x, double *v);
    virtual void writeOutputLine(const double cpuTimeUsed);
    virtual void writeSummary();

protected:
    void acceleration(double *__restrict__ x, double *__restrict__ a);
    double radiusOfGyration(double *__restrict__ x);

    int m_Nobj;
    int m_Npad;
    double m_timeStep;

    double *m_masses;

    System *m_sys;
};

}

#endif // PROPAGATOR_H
