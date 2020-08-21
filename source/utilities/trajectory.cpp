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


#include "trajectory.h"

Trajectory::Trajectory(std::string fileName)
{
    m_trjFile = std::ifstream(fileName, std::ios::binary);

    // determine the length of the file in Bytes
    m_trjFile.seekg(0, m_trjFile.end);
    int length = m_trjFile.tellg();
    m_trjFile.seekg(0, m_trjFile.beg);

    m_trjFile.read(reinterpret_cast<char*>(&m_nObj), sizeof(int));

    // determine the number of frames
    //
    //        subtract header from size
    //                                    a frame contains 6 doubles
    //                                    per object + time
    //        |-----------------------|   |---------------|
    m_nFrames = (length - 4 - (8 * m_nObj)) / ((48 * m_nObj) + 8);
}

int Trajectory::nObj()
{
    return m_nObj;
}

int Trajectory::nFrames()
{
    return m_nFrames;
}

std::vector<double> Trajectory::readMasses()
{
    double *buffer = new double[m_nObj];
    m_trjFile.read(reinterpret_cast<char*>(buffer), m_nObj * sizeof(double));
    std::vector<double> masses;

    for (int i = 0; i < m_nObj; i++)
    {
        masses.push_back(buffer[i]);
    }
    delete[] buffer;

    return masses;
}

TrajectoryFrame Trajectory::readNextFrame()
{
    double elapsedTime;
    m_trjFile.read(reinterpret_cast<char*>(&elapsedTime), sizeof(double));

    double *buffer = new double[3 * m_nObj];
    m_trjFile.read(reinterpret_cast<char*>(buffer), 3 * m_nObj * sizeof(double));
    std::vector<double> newX;
    for (int i = 0; i < m_nObj; i++)
    {
        newX.push_back(buffer[i]);
    }

    m_trjFile.read(reinterpret_cast<char*>(buffer), 3 * m_nObj * sizeof(double));
    std::vector<double> newV;
    for (int i = 0; i < m_nObj; i++)
    {
        newV.push_back(buffer[i]);
    }
    delete[] buffer;

    TrajectoryFrame frame;
    frame.setElapsedTime(elapsedTime);
    frame.setX(newX);
    frame.setV(newV);

    return frame;
}
