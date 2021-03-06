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


#ifndef TRAJECTORYFRAME_H
#define TRAJECTORYFRAME_H

#include <vector>


class TrajectoryFrame
{
public:
    TrajectoryFrame();

    double getElapsedTime();
    std::vector<double> getX();
    std::vector<double> getV();

    void setElapsedTime(const double time);
    void setX(std::vector<double> newX);
    void setV(std::vector<double> newV);

private:
    double elapsedTime;
    std::vector<double> x;
    std::vector<double> v;
};

#endif // TRAJECTORYFRAME_H
