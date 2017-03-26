/*
 * Copyright 2012-2017 Jan von Cosel
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

#include<math.h>

void acceleration_c_(int *Nobj, double *x, double *a, double *G, double *mass)
{
#pragma omp parallel
{
    double aPriv[3 * *Nobj];

    int i, j;

    #pragma omp for
    for (j = 0; j < 3 * *Nobj; j++)
    {
        a[j] = 0.0;
    }

    for (j = 0; j < 3 * *Nobj; j++)
    {
        aPriv[j] = 0.0;
    }


    #pragma omp for schedule(guided)
    for (i = 0; i < *Nobj - 1; i++)
    {
        for (j = i + 1; j < *Nobj; j++)
        {
            double dX = x[3 * j + 0] - x[3 * i + 0];
            double dY = x[3 * j + 1] - x[3 * i + 1];
            double dZ = x[3 * j + 2] - x[3 * i + 2];
            double R2 = dX * dX + dY * dY + dZ * dZ;
            double mFac = mass[i] * mass[j] / (R2 * sqrt(R2));
            aPriv[3 * i + 0] += mFac * dX;
            aPriv[3 * i + 1] += mFac * dY;
            aPriv[3 * i + 2] += mFac * dZ;
            aPriv[3 * j + 0] -= mFac * dX;
            aPriv[3 * j + 1] -= mFac * dY;
            aPriv[3 * j + 2] -= mFac * dZ;
        }
    }

    #pragma omp critical
    {
        for (i = 0; i < 3 * *Nobj; i++)
            a[i] += aPriv[i];
    }

    #pragma omp barrier


    #pragma omp for
    for (i = 0; i < *Nobj; i++)
    {
        double fFac = *G / mass[i];
        a[3 * i + 0] *= fFac;
        a[3 * i + 1] *= fFac;
        a[3 * i + 2] *= fFac;
    }
}
}


void acceleration_c2_(int *Nobj, double *x, double *a, double *G, double *mass)
{
    #pragma omp parallel
    {
        double aPriv[3 * *Nobj];

        int i, j;

        #pragma omp for
        for (j = 0; j < 3 * *Nobj; j++)
            a[j] = 0.0;

        for (j = 0; j < 3 * *Nobj; j++)
            aPriv[j] = 0.0;

        #pragma omp for schedule(static)
        for (i = 0; i < *Nobj - 1; i++)
            for (j = i + 1; j < *Nobj; j++)
            {
                double dX = x[3 * j + 0] - x[3 * i + 0];
                double dY = x[3 * j + 1] - x[3 * i + 1];
                double dZ = x[3 * j + 2] - x[3 * i + 2];
                double R2 = dX * dX + dY * dY + dZ * dZ;
                double tmpFac = *G / (R2 * sqrt(R2));
                aPriv[3 * i + 0] += mass[j] * tmpFac * dX;
                aPriv[3 * i + 1] += mass[j] * tmpFac * dY;
                aPriv[3 * i + 2] += mass[j] * tmpFac * dZ;
                aPriv[3 * j + 0] -= mass[i] * tmpFac * dX;
                aPriv[3 * j + 1] -= mass[i] * tmpFac * dY;
                aPriv[3 * j + 2] -= mass[i] * tmpFac * dZ;
            }

        #pragma omp critical
            for (i = 0; i < 3 * *Nobj; i++)
                a[i] += aPriv[i];
    }
}

