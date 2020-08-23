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


#include "trajectoryframe.h"

TrajectoryFrame::TrajectoryFrame()
{
}

double TrajectoryFrame::getElapsedTime()
{
    return elapsedTime;
}

std::vector<double> TrajectoryFrame::getX()
{
    return x;
}

std::vector<double> TrajectoryFrame::getV()
{
    return v;
}

void TrajectoryFrame::setElapsedTime(const double time)
{
    elapsedTime = time;
}

void TrajectoryFrame::setX(std::vector<double> newX)
{
    x = newX;
}

void TrajectoryFrame::setV(std::vector<double> newV)
{
    v = newV;
}
