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


#include <vector>
#include <fstream>
#include "constants.h"

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
    int init(int argnum, char *arguments[]);
    std::string get_testString();

    std::string &inputFile();
    std::string baseName();
    double tfinal();
    double tout();
    double minStep();
    double Eps();
    double RedMin();
    double RedMax();
    double MaxInc();
    std::ofstream &outputFile();
    std::ofstream &binTrjFile();
    std::ofstream &txtTrjFile();
    std::ofstream &stepsFile();
    std::ofstream &restartFile();
    int Nobj();
    int Npad();
    int MaxSubStep();
    int IncThres();
    int Ndigit();
    std::vector<double> XX0();
    std::vector<double> XY0();
    std::vector<double> XZ0();
    std::vector<double> VX0();
    std::vector<double> VY0();
    std::vector<double> VZ0();
    std::vector<std::string> names();
    std::vector<double> masses();
    bool Verbose();
    bool ShiftCOM();
    bool ShiftMom();
    bool TextTrj();
    bool Restart();
    bool Steps();
    bool UnRes();
    double TotMass();
    void listParas();
    IntType intType();
    double initStep();

    void inc_BSL();
    void inc_BSO();
    void inc_BSS();
    void inc_BSE();
private:
    Configuration(){}
    Configuration(const Configuration&);
    Configuration& operator=(const Configuration&);

    void setDefaults();
    int parseInputLine(std::string &inputLine);

    std::string m_inputFileName;
    std::string m_baseName;
    std::string m_outputFileName;
    std::string m_stepsName;

    // here comes the data:
    int m_Nobj;         // number of objects
    int m_Npad;         // Nobj increased to be divisible by 4 (AVX register length)
    int m_MaxSubStep;   // max. no. of substeps in one BS step
    int m_IncThres;     // no. of substeps below which the stepsize will be increased
    int m_Ndigit;       // no. of significant digits in text trajectory
    int m_Nstep;        // no. of steps in rk4fix

    // counters for function calls
    int m_N_BS_LargeStep;
    int m_N_BS_OneStep;
    int m_N_BS_SubSteps;
    int m_N_BS_Extrapolate;

    double m_tFinal;    // total length of simulation
    double m_tOut;      // writeout interval
    double m_eps;       // error tolerance for the propagation
    double m_epsThres;  // threshold to increase stepsize in Runge-Kutta run
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

    std::vector<double> m_masses;       // masses of the objects
    std::vector<std::string> m_names;   // names of the objects
    std::vector<double> m_X0;           // initial x position
    std::vector<double> m_Y0;           // initial y position
    std::vector<double> m_Z0;           // initial z position
    std::vector<double> m_VX0;          // initial velocity along x
    std::vector<double> m_VY0;          // initial velocity along y
    std::vector<double> m_VZ0;          // initial velocity along z

    IntType m_IntType;

    std::ofstream m_outputFile;         // stream corresponding to the output file
    std::ofstream m_binTrjFile;         // stream corresponding to the trajectory file
    std::ofstream m_txtTrjFile;         // stream corresponding to the text trajectory file
    std::ofstream m_stepsFile;          // stream corresponding to the steps file
    std::ofstream m_restartFile;        // stream corresponding to the restart file
};

} // namespace Astronx

#endif // INPUT_H

