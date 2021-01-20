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


#ifndef LEAPFROG_H
#define LEAPFROG_H

#include "propagator.h"


namespace Astronx
{

/**
 * @brief The LeapFrog class
 */
class LeapFrog : public Propagator
{
public:
    LeapFrog(System *sys);
    ~LeapFrog();

    double largeStep(double *x, double *v);

    void writeOutputLine(const double cpuTimeUsed);
    void writeSummary();

private:
    double m_stepSize;
    double m_internalElapsedTime;

    double *m_xLFtmp;
    double *m_vLFtmp;

    double *m_v12;
    double *m_a0;
    double *m_x2;
};

}

#endif // LEAPFROG_H
