#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include <vector>
#include <string>
#include <fstream>
#include <array>


namespace Astronx
{

class Propagator
{
public:
    Propagator();
    std::array<double,3> com();
    std::array<double,3> linMom();
    std::array<double,3> angMom();
    void shiftCom();
    void shiftMom();

    void propagate();

private:
    void propBS();
    void writeToTrj();

    // general properties of the system
    int m_Npad;
    int m_Nobj;
    double m_totMass;

    // the time
    double m_elapsedTime;

    // position and velocity of objects
    double *m_xx;
    double *m_xy;
    double *m_xz;
    double *m_vx;
    double *m_vy;
    double *m_vz;

    double *m_masses;   // masses of objects

    std::vector<std::string> m_names;   // object names

    std::ofstream m_restartFile;
    std::ofstream m_stepsFile;
    std::ofstream m_binTrj;
    std::ofstream m_txtTrj;
};

}

#endif // PROPAGATOR_H
