#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include <vector>
#include <string>
#include <fstream>
#include <array>


namespace Astronx
{

class OldPropagator
{
public:
    OldPropagator();
    std::array<double,3> com();
    std::array<double,3> linMom();
    std::array<double,3> angMom();
    double radiusOfGyration(double *__restrict__ x);
    void shiftCom();
    void shiftMom();

    void propagate();

private:
    void writeToTrj();
    void writeStatus();
    void BS_LargeStep();
    int BS_OneStep();
    void BS_SubSteps(const int nsteps, const double stepSize);
    void BS_Extrapolate(const int i_est, const double h_est);
    void RK4F_LargeStep();
    void acceleration(double *__restrict__ x, double *__restrict__ a);

    // general properties of the system
    int m_Npad;
    int m_Nobj;
    int m_nsteps;
    double m_totMass;

    double m_elapsedTime;   // the time
    double m_timeStep;      // current (trial) integration time step
    double m_doneStep;      // last successful integration step
    double m_delta;         // current integration error

    // integration step counters
    int m_N_ok;
    int m_N_fail;
    int m_N_ok_tot;
    int m_N_fail_tot;
    int m_N_bssteps;
    int m_N_bstotal;
    int m_N_rkfixtotal;
    int m_N_smallsteps;
    int m_N_smalltotal;

    bool m_underflow;

    // position and velocity of objects, including temporary arrays
    double *m_x_BSlarge;
    double *m_v_BSlarge;
    double *m_x_BSLtmp;
    double *m_v_BSLtmp;
    double *m_a_BSstart;
    double *m_x_SubStep;
    double *m_x_SubFin;
    double *m_v_SubFin;
    double *m_a_SubInt;

    // auxiliary data for Bulirsch-Stoer extrapolation
    double *m_extD;
    double *m_extErr;
    double *m_tmpDat;
    std::vector<double> m_extH;

    double *m_masses;   // masses of objects

    std::vector<std::string> m_names;   // object names

    std::ofstream m_restartFile;
    std::ofstream m_stepsFile;
    std::ofstream m_binTrj;
    std::ofstream m_txtTrj;
};

}

#endif // PROPAGATOR_H
