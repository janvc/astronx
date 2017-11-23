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
 * along with astronx. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef INPUT_H
#define INPUT_H


namespace Astronx
{

class Configuration
{
public:
    static Configuration &get()
    {
        static Configuration instance;
        return instance;
    }
    void init(int argnum, char *arguments[]);
    std::string get_testString();

private:
    Configuration(){}
    Configuration(const Configuration&);
    Configuration& operator=(const Configuration&);

    std::string m_inputFileName;

    // here comes the data:
    int m_Nobj;         // number of objects
    int m_MaxSubStep;   // max. no. of substeps in one BS step
    int m_IncThres;     // no. of substeps below which the stepsize will be increased
    int m_Ndigit;       // no. of significant digits in text trajectory
    int m_Nstep;        // no. of steps in rk4fix

    double m_tFinal;    // total length of simulation
    double m_tOut;      // writeout interval
    double m_eps;       // error tolerance for the propagation
    double m_MinStep;   // minimum timestep
    double m_MaxInc;    // max. factor by which to increase stepsize
    double m_RedMin;    // minimum factor for the stepsize reduction
    double m_RedMax;    // maximum reduction factor
    double m_InitStep;  // initial timestep
    double m_Mtot;      // total mass of the system

    bool m_restart;     // create a restart file at every write step
    bool m_steps;       // create a steps file
    bool m_textTrj;     // create a text trajectory file
    bool m_UnResProp;   // do an unrestricted propagation
    bool m_overwrite;   // overwrite previous results
    bool m_ShiftCOM;    // shift the center of mass to the origin
    bool m_ShiftMom;    // eliminate the system's total linear momentum
    bool m_verbose;     // write info to the terminal
};

} // namespace Astronx

#endif // INPUT_H

