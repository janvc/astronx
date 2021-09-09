#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace PhyCon
{

const double G = 6.6726e-11;

}

namespace Integrators {

enum IntType
{
    BS = 1,
    RK4 = 2,
    LeapFrog = 3
};

}

#endif // CONSTANTS_H
