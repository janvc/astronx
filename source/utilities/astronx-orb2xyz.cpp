/*
 * Copyright 2015-2020 Jan von Cosel
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

#include <string>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>

int main(int argc, char *argv[])
{
    const double G = 6.6726e-11;

    double periapsis;
    double eccentricity;
    double semiMajorAxis;
    double centralMass;
    double trueAnomaly;
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "produce this help message")
            ("periapsis,p", po::value<double>(&periapsis), "periapsis of orbit")
            ("eccentricity,e", po::value<double>(&eccentricity), "eccentricity of orbit")
            ("trueanomaly,t", po::value<double>(&trueAnomaly), "true anomaly of orbit")
            ("semimajoraxis,a", po::value<double>(&semiMajorAxis), "semi-major axis of orbit")
            ("centralmass,m", po::value<double>(&centralMass), "mass of the central body")
            ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // if the periapsis is not given, calculate it
    if (! vm.count("periapsis"))
    {
        if (! (vm.count("semimajoraxis") and vm.count("eccentricity")))
        {
            std::cerr << "Periapsis not given. We need the semi-major axis and the eccentricity to calculate the position\n";
            return 1;
        }

        periapsis = semiMajorAxis * (1.0 - eccentricity);
    }

    /*
     * the distance between the focal point (central body) and the satelite is given by
     *
     *             1 - e**2
     *    r = a --------------
     *           1 + e cos(n)
     *
     * with the semi-major axis a, the eccentricity e and the true anomaly n
     */
    double distance = semiMajorAxis * (1.0 - eccentricity * eccentricity)
                                    / (1.0 + eccentricity * std::cos(trueAnomaly));

    // calculate the speed at the periapsis
    if (! vm.count("centralmass"))
    {
        std::cerr << "Mass of central object not given\n";
        return 1;
    }

    /*
     * the speed along an eliptical depends on the central body's mass M and the distance r
     *           ___________________
     *          /     (  2     1  )
     *    v =  /  G M ( --- - --- )
     *       |/       (  r     a  )
     *
     * where G is the gravitational constant and a is the semi-major-axis
     */
    double speed = std::sqrt(G * centralMass * ((2.0 / distance) - (1.0 / semiMajorAxis)));

    std::cout << "name mass "
              << distance << " " << 0 << " " << 0 << " "
              << 0 << " " << speed << " " << 0
              << std::endl;

    return 0;
}
