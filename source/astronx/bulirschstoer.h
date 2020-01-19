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


#ifndef BULIRSCHSTOER_H
#define BULIRSCHSTOER_H

#include "propagator.h"


namespace Astronx
{

class BulirschStoer : public Propagator
{
public:
    BulirschStoer(const int Npad, System *sys);
    ~BulirschStoer();

    double largeStep(double *x, double *v);
    bool oneStep();
    void subSteps(const int nSteps, const double stepSize);
    void extrapolate(const int i_est, const double h_est);

    void writeOutputLine();

private:
//    int m_Nobj;
//    int m_Npad;
    int m_N_ok;
    int m_N_fail;
    int m_nsteps;

    double m_delta;
    double m_doneStep;

    double *m_x_BSLtmp;
    double *m_v_BSLtmp;
    double *m_x_SubStep;
    double *m_v_SubStep;
    double *m_x_SubFin;
    double *m_v_SubFin;
    double *m_a_BSStart;
    double *m_a_SubInt;

    // auxiliary data for Bulirsch-Stoer extrapolation
    double *m_extD;
    double *m_extErr;
    double *m_tmpDat;
    double *m_extC;
    std::vector<double> m_extH;

    std::ofstream m_stepsFile;
};

}

#endif // BULIRSCHSTOER_H
