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


#ifndef RUNGEKUTTA4_H
#define RUNGEKUTTA4_H

#include "propagator.h"


namespace Astronx
{

class RungeKutta4 : public Propagator
{
public:
    RungeKutta4(System *sys);
    ~RungeKutta4();

    double largeStep(double *x, double *v);

    void writeOutputLine(const double cpuTimeUsed);
    void writeSummary();

private:
    int m_nSteps;

    double m_stepSize;
    double m_internalElapsedTime;

    double *m_xRKtmp;
    double *m_vRKtmp;
    double *m_a0;

    double *m_x1;
    double *m_v1;
    double *m_a1;

    double *m_x2;
    double *m_v2;
    double *m_a2;

    double *m_x3;
    double *m_v3;
    double *m_a3;

    double *m_x4;
    double *m_v4;
    double *m_a4;
};

}

#endif // RUNGEKUTTA4_H
