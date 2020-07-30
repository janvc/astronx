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
#include <eigen3/Eigen/Core>

Eigen::Matrix3d euler2rot(const double psi, const double theta, const double phi)
{
    Eigen::Matrix3d rot;

    rot(0,0) =  std::cos(psi) * std::cos(phi) - std::cos(theta) * std::sin(phi) * std::sin(psi);
    rot(0,1) =  std::sin(psi) * std::cos(phi) + std::cos(theta) * std::sin(phi) * std::cos(psi);
    rot(0,2) =  std::sin(phi) * std::sin(theta);
    rot(1,0) = -std::cos(psi) * std::sin(phi) - std::cos(theta) * std::cos(phi) * std::sin(psi);
    rot(1,1) = -std::sin(psi) * std::sin(phi) + std::cos(theta) * std::cos(phi) * std::cos(psi);
    rot(1,2) =  std::cos(phi) * std::sin(theta);
    rot(2,0) =  std::sin(theta) * std::sin(psi);
    rot(2,1) = -std::sin(theta) * std::cos(psi);
    rot(2,2) =  std::cos(theta);

    return rot;
}


int main(int argc, char *argv[])
{
    // TODO:
    // 1) set up position vector in local (orbit) coordinate system
    // 2) get transformation matrix to global coordinate system
    // 3) rotate position vector
    // 4) deal with velocity
    //
    // do we actually need the periapsis?

    const double G = 6.6726e-11;
    const double deg2rad = M_PI / 180.0;

    double periapsis;
    double eccentricity = 0.0;
    double semiMajorAxis;
    double centralMass;
    double trueAnomaly = 0.0;
    double inclination = 0.0;
    double longitudeOfAscendingNode = 0.0;
    double argumentOfPeriapsis = 0.0;

    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "produce this help message")
            ("centralmass,m", po::value<double>(&centralMass)->required(), "mass of the central body")
            ("semimajoraxis,a", po::value<double>(&semiMajorAxis)->required(), "semi-major axis of orbit")
            ("periapsis,p", po::value<double>(&periapsis), "periapsis of orbit")
            ("eccentricity,e", po::value<double>(&eccentricity), "eccentricity of orbit")
            ("trueanomaly,t", po::value<double>(&trueAnomaly), "true anomaly of orbit")
            ("inclination,i", po::value<double>(&inclination), "inclination of the orbit")
            ("lonasc,l", po::value<double>(&longitudeOfAscendingNode), "longitude of the orbit's ascending node")
            ("argper,b", po::value<double>(&argumentOfPeriapsis), "argument (angle) of orbit's periapsis")
            ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);


    /*
     * transform the angles to radians
     */
    trueAnomaly *= deg2rad;
    inclination *= deg2rad;
    longitudeOfAscendingNode *= deg2rad;
    argumentOfPeriapsis *= deg2rad;


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


    /*
     * calculate the position vector in the orbital coordinage system
     */
    Eigen::Vector3d localPosition(
                distance * std::cos(trueAnomaly),   // x coordinate
                distance * std::sin(trueAnomaly),   // y coordinate
                0.0);                               // z coordinate is always zero


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


    /*
     * the flight path angle gamma is the angle between the velocity vector and the vector
     * that is perpendicular to the position vector
     */
    double tanGamma = eccentricity * std::sin(trueAnomaly) / (1.0 + eccentricity * std::cos(trueAnomaly));
    double gamma = std::atan(tanGamma);


    /*
     * the angle between the velocity vector and the x axis is given by
     *
     * pi/4 + theta - gamma
     */
    double velAngle = M_PI / 2.0 + trueAnomaly - gamma;

    Eigen::Vector3d localVelocity(
                speed * std::cos(velAngle),
                speed * std::sin(velAngle),
                0.0);


    /*
     * create the rotation matrix and transform the vectors
     */
    Eigen::Matrix3d rotMat = euler2rot(longitudeOfAscendingNode, inclination, argumentOfPeriapsis);

    Eigen::Vector3d globalPosition = rotMat.transpose() * localPosition;

    Eigen::Vector3d globalVelocity = rotMat.transpose() * localVelocity;

    std::cout << "name mass "
              << globalPosition(0) << " " << globalPosition(1) << " " << globalPosition(2) << " "
              << globalVelocity(0) << " " << globalVelocity(1) << " " << globalVelocity(2)
              << std::endl;



    return 0;
}
