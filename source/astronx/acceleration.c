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
#include<x86intrin.h>

/*
 * basic version with two loops
 */
void acc_1c_(int *restrict Nobj, double *restrict x, double *restrict a, double *restrict G, double *restrict mass)
{
    int i, j;

    for (j = 0; j < 3 * *Nobj; j++)
        a[j] = 0.0;

    for (i = 0; i < *Nobj - 1; i++)
    {
        for (j = i + 1; j < *Nobj; j++)
        {
            double dX = x[3 * j + 0] - x[3 * i + 0];
            double dY = x[3 * j + 1] - x[3 * i + 1];
            double dZ = x[3 * j + 2] - x[3 * i + 2];
            double R2 = dX * dX + dY * dY + dZ * dZ;
            double mFac = mass[i] * mass[j] / (R2 * sqrt(R2));
            a[3 * i + 0] += mFac * dX;
            a[3 * i + 1] += mFac * dY;
            a[3 * i + 2] += mFac * dZ;
            a[3 * j + 0] -= mFac * dX;
            a[3 * j + 1] -= mFac * dY;
            a[3 * j + 2] -= mFac * dZ;
        }
    }

    for (i = 0; i < *Nobj; i++)
    {
        double fFac = *G / mass[i];
        a[3 * i + 0] *= fFac;
        a[3 * i + 1] *= fFac;
        a[3 * i + 2] *= fFac;
    }
}

/*
 * basic version with one loop
 */
void acc_2c_(int *restrict Nobj, double *restrict x, double *restrict a, double *restrict G, double *restrict mass)
{
    int i, j;

    for (j = 0; j < 3 * *Nobj; j++)
        a[j] = 0.0;

    for (i = 0; i < *Nobj - 1; i++)
        for (j = i + 1; j < *Nobj; j++)
        {
            double dX = x[3 * j + 0] - x[3 * i + 0];
            double dY = x[3 * j + 1] - x[3 * i + 1];
            double dZ = x[3 * j + 2] - x[3 * i + 2];
            double R2 = dX * dX + dY * dY + dZ * dZ;
            double tmpFac = *G / (R2 * sqrt(R2));
            a[3 * i + 0] += mass[j] * tmpFac * dX;
            a[3 * i + 1] += mass[j] * tmpFac * dY;
            a[3 * i + 2] += mass[j] * tmpFac * dZ;
            a[3 * j + 0] -= mass[i] * tmpFac * dX;
            a[3 * j + 1] -= mass[i] * tmpFac * dY;
            a[3 * j + 2] -= mass[i] * tmpFac * dZ;
        }
}

/*
 * three arrays for position and acceleration, one loop
 */
void acc_2c2_(int *restrict Nobj,
              double *restrict xx, double *restrict xy, double *restrict xz,
              double *restrict ax, double *restrict ay, double *restrict az,
              double *restrict G, double *restrict mass)
{
    int i, j;

    for (j = 0; j < *Nobj; j++)
    {
        ax[j] = 0.0;
        ay[j] = 0.0;
        az[j] = 0.0;
    }

    for (i = 0; i < *Nobj - 1; i++)
        for (j = i + 1; j < *Nobj; j++)
        {
            double dX = xx[j] - xx[i];
            double dY = xy[j] - xy[i];
            double dZ = xz[j] - xz[i];
            double R2 = dX * dX + dY * dY + dZ * dZ;
            double tmpFac = *G / (R2 * sqrt(R2));
            ax[i] += mass[j] * tmpFac * dX;
            ay[i] += mass[j] * tmpFac * dY;
            az[i] += mass[j] * tmpFac * dZ;
            ax[j] -= mass[i] * tmpFac * dX;
            ay[j] -= mass[i] * tmpFac * dY;
            az[j] -= mass[i] * tmpFac * dZ;
        }
}

/*
 * three arrays, avx intrinsics
 */
void acc_2c2avx_(int *restrict Nobj,
                 double *restrict xx, double *restrict xy, double *restrict xz,
                 double *restrict ax, double *restrict ay, double *restrict az,
                 double *restrict G, double *restrict mass)
{
    int i, j, k;

    for (j = 0; j < *Nobj; j++)
    {
        ax[j] = 0.0;
        ay[j] = 0.0;
        az[j] = 0.0;
    }

    for (i = 0; i < *Nobj; i += 4)
    {
        for (k = 0; k < 4; k++)
        {
            /*
             * First, do the diagonal part, that cannot easily be treated
             * with vector instructions:
             *
             *   0_1_2_3__4________8________ ...
             * 0|0 x x x |y y y y |y y y y | ...
             * 1|z 0 x x |y y y y |y y y y | ...
             * 2|z z 0 x |y y y y |y y y y | ...
             * 3|z z z 0 |y y y y |y y y y | ...
             *  |--------|--------|--------| ...
             * 4|        |0 x x x |y y y y | ...
             *  |        |z 0 x x |y y y y | ...
             *  |        |z z 0 x |y y y y | ...
             *  |        |z z z 0 |y y y y | ...
             *  |--------|--------|--------| ...
             * 8|        |        |0 x x x |
             *  .        .        .        .
             *  .        .        .        .
             *  .        .        .        .
             *
             */
//            for (j = i + k + 1; j < i + 4; j++)
//            {
//                double dX = xx[j] - xx[i+k];
//                double dY = xy[j] - xy[i+k];
//                double dZ = xz[j] - xz[i+k];
//                double R2 = dX * dX + dY * dY + dZ * dZ;
//                double tmpFac = *G / (R2 * sqrt(R2));
//                ax[i+k] += mass[j] * tmpFac * dX;
//                ay[i+k] += mass[j] * tmpFac * dY;
//                az[i+k] += mass[j] * tmpFac * dZ;
//                ax[j] -= mass[i+k] * tmpFac * dX;
//                ay[j] -= mass[i+k] * tmpFac * dY;
//                az[j] -= mass[i+k] * tmpFac * dZ;
//            }
            if (k == 0)
            {
                // interaction between objects 0 and 1
                double dX = xx[i+1] - xx[i];
                double dY = xy[i+1] - xy[i];
                double dZ = xz[i+1] - xz[i];
                double R2 = dX * dX + dY * dY + dZ * dZ;
                double tmpFac = *G / (R2 * sqrt(R2));
                ax[i] += mass[i+1] * tmpFac * dX;
                ay[i] += mass[i+1] * tmpFac * dY;
                az[i] += mass[i+1] * tmpFac * dZ;
                ax[i+1] -= mass[i] * tmpFac * dX;
                ay[i+1] -= mass[i] * tmpFac * dY;
                az[i+1] -= mass[i] * tmpFac * dZ;

                // interaction between objects 0 and 2/3
                __m128d x0, y0, z0, x1, y1, z1, x2, y2, z2;
                __m128d tmp0, tmp1;

                x0 = _mm_load_pd1(xx + i);  // load the position of object 'i' into '0'
                y0 = _mm_load_pd1(xy + i);
                z0 = _mm_load_pd1(xz + i);
                x1 = _mm_load_pd(xx + i + 2);   // load objects 'i+2' and 'i+3' into '1'
                y1 = _mm_load_pd(xy + i + 2);
                z1 = _mm_load_pd(xz + i + 2);

                x1 = _mm_sub_pd(x1, x0);    // put the distances into '1'
                y1 = _mm_sub_pd(y1, y0);
                z1 = _mm_sub_pd(z1, z0);
                x0 = _mm_mul_pd(x1, x1);    // and their squares into '0'
                y0 = _mm_mul_pd(y1, y1);
                z0 = _mm_mul_pd(z1, z1);

                tmp0 = _mm_add_pd(x0, y0);  // r2 = dx^2 + dy^2 + dz^2
                tmp0 = _mm_add_pd(tmp0, z0);

                tmp1 = _mm_sqrt_pd(tmp0);   // fac = G / (r2 * sqrt(r2))
                tmp1 = _mm_mul_pd(tmp0, tmp1);
                tmp0 = _mm_load_pd1(G);
                tmp0 = _mm_div_pd(tmp0, tmp1);  // 'tmp0' contains the conversion factor

                x0 = _mm_mul_pd(x1, tmp0);  // put 'factor * dx' into '0'
                y0 = _mm_mul_pd(y1, tmp0);
                z0 = _mm_mul_pd(z1, tmp0);

                tmp0 = _mm_load_pd1(mass + i);      // mass of 'i' into 'tmp0'
                tmp1 = _mm_load_pd(mass + i + 2);   // massES of 'i+2' and 'i+3' into 'tmp1'

                x1 = _mm_mul_pd(x0, tmp1);  // multiply with massES of 'i+2' and 'i+3'
                y1 = _mm_mul_pd(y0, tmp1);
                z1 = _mm_mul_pd(z0, tmp1);
                x1 = _mm_hadd_pd(x1, x1);   // add the two values together
                y1 = _mm_hadd_pd(y1, y1);
                z1 = _mm_hadd_pd(z1, z1);
                x2 = _mm_load_pd1(ax + i);  // load acceleration of reference object
                y2 = _mm_load_pd1(ay + i);
                z2 = _mm_load_pd1(az + i);
                x2 = _mm_add_pd(x2, x1);    // add to it
                y2 = _mm_add_pd(y2, y1);
                z2 = _mm_add_pd(z2, z1);
                _mm_store_sd(ax + i, x2);   // and store it again
                _mm_store_sd(ay + i, y2);
                _mm_store_sd(az + i, z2);

                x1 = _mm_mul_pd(x0, tmp0);  // now multiply with mass of reference
                y1 = _mm_mul_pd(y0, tmp0);
                z1 = _mm_mul_pd(z0, tmp0);
                x2 = _mm_load_pd(ax + i + 2);   // load acceleration of other objects
                y2 = _mm_load_pd(ay + i + 2);
                z2 = _mm_load_pd(az + i + 2);
                x2 = _mm_sub_pd(x2, x1);        // subtract from them
                y2 = _mm_sub_pd(y2, y1);
                z2 = _mm_sub_pd(z2, z1);
                _mm_store_pd(ax + i + 2, x2);   // and store again
                _mm_store_pd(ay + i + 2, y2);
                _mm_store_pd(az + i + 2, z2);
            }
            else if (k == 1)
            {
                // interaction between objects 1 and 2/3
                __m128d x0, y0, z0, x1, y1, z1, x2, y2, z2;
                __m128d tmp0, tmp1;

                x0 = _mm_load_pd1(xx + i + 1);  // load the position of object 'i+1' into '0'
                y0 = _mm_load_pd1(xy + i + 1);
                z0 = _mm_load_pd1(xz + i + 1);
                x1 = _mm_load_pd(xx + i + 2);   // load objects 'i+2' and 'i+3' into '1'
                y1 = _mm_load_pd(xy + i + 2);
                z1 = _mm_load_pd(xz + i + 2);

                x1 = _mm_sub_pd(x1, x0);    // put the distances into '1'
                y1 = _mm_sub_pd(y1, y0);
                z1 = _mm_sub_pd(z1, z0);
                x0 = _mm_mul_pd(x1, x1);    // and their squares into '0'
                y0 = _mm_mul_pd(y1, y1);
                z0 = _mm_mul_pd(z1, z1);

                tmp0 = _mm_add_pd(x0, y0);  // r2 = dx^2 + dy^2 + dz^2
                tmp0 = _mm_add_pd(tmp0, z0);

                tmp1 = _mm_sqrt_pd(tmp0);   // fac = G / (r2 * sqrt(r2))
                tmp1 = _mm_mul_pd(tmp0, tmp1);
                tmp0 = _mm_load_pd1(G);
                tmp0 = _mm_div_pd(tmp0, tmp1);  // 'tmp0' contains the conversion factor

                x0 = _mm_mul_pd(x1, tmp0);  // put 'factor * dx' into '0'
                y0 = _mm_mul_pd(y1, tmp0);
                z0 = _mm_mul_pd(z1, tmp0);

                tmp0 = _mm_load_pd1(mass + i + 1);      // mass of 'i+1' into 'tmp0'
                tmp1 = _mm_load_pd(mass + i + 2);   // massES of 'i+2' and 'i+3' into 'tmp1'

                x1 = _mm_mul_pd(x0, tmp1);  // multiply with massES of 'i+2' and 'i+3'
                y1 = _mm_mul_pd(y0, tmp1);
                z1 = _mm_mul_pd(z0, tmp1);
                x1 = _mm_hadd_pd(x1, x1);   // add the two values together
                y1 = _mm_hadd_pd(y1, y1);
                z1 = _mm_hadd_pd(z1, z1);
                x2 = _mm_load_pd1(ax + i + 1);  // load acceleration of reference object
                y2 = _mm_load_pd1(ay + i + 1);
                z2 = _mm_load_pd1(az + i + 1);
                x2 = _mm_add_pd(x2, x1);    // add to it
                y2 = _mm_add_pd(y2, y1);
                z2 = _mm_add_pd(z2, z1);
                _mm_store_sd(ax + i + 1, x2);   // and store it again
                _mm_store_sd(ay + i + 1, y2);
                _mm_store_sd(az + i + 1, z2);

                x1 = _mm_mul_pd(x0, tmp0);  // now multiply with mass of reference
                y1 = _mm_mul_pd(y0, tmp0);
                z1 = _mm_mul_pd(z0, tmp0);
                x2 = _mm_load_pd(ax + i + 2);   // load acceleration of other objects
                y2 = _mm_load_pd(ay + i + 2);
                z2 = _mm_load_pd(az + i + 2);
                x2 = _mm_sub_pd(x2, x1);        // subtract from them
                y2 = _mm_sub_pd(y2, y1);
                z2 = _mm_sub_pd(z2, z1);
                _mm_store_pd(ax + i + 2, x2);   // and store again
                _mm_store_pd(ay + i + 2, y2);
                _mm_store_pd(az + i + 2, z2);
            }
            else if (k == 2)
            {
                // interaction between objects 2 and 3
                double dX = xx[i+3] - xx[i+2];
                double dY = xy[i+3] - xy[i+2];
                double dZ = xz[i+3] - xz[i+2];
                double R2 = dX * dX + dY * dY + dZ * dZ;
                double tmpFac = *G / (R2 * sqrt(R2));
                ax[i+2] += mass[i+3] * tmpFac * dX;
                ay[i+2] += mass[i+3] * tmpFac * dY;
                az[i+2] += mass[i+3] * tmpFac * dZ;
                ax[i+3] -= mass[i+2] * tmpFac * dX;
                ay[i+3] -= mass[i+2] * tmpFac * dY;
                az[i+3] -= mass[i+2] * tmpFac * dZ;
            }

            /*
             * Now do the remaining part of the row with AVX functions
             */
            __m256d X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2;
            __m256d Tmp0, Tmp1;
            __m256d AccIX, AccIY, AccIZ;

            // put the position of object "i + k" into all
            // four slots of the AVX registers
            X0 = _mm256_broadcast_sd(xx + i + k);
            Y0 = _mm256_broadcast_sd(xy + i + k);
            Z0 = _mm256_broadcast_sd(xz + i + k);
            AccIX = _mm256_setzero_pd();
            AccIY = _mm256_setzero_pd();
            AccIZ = _mm256_setzero_pd();

            for (j = i + 4; j < *Nobj; j += 4)
            {
                // dX = xx[j] - xx[i]
                X1 = _mm256_load_pd(xx + j);
                Y1 = _mm256_load_pd(xy + j);
                Z1 = _mm256_load_pd(xz + j);
                X1 = _mm256_sub_pd(X1, X0);
                Y1 = _mm256_sub_pd(Y1, Y0);
                Z1 = _mm256_sub_pd(Z1, Z0);
                X0 = _mm256_mul_pd(X1, X1);
                Y0 = _mm256_mul_pd(Y1, Y1);
                Z0 = _mm256_mul_pd(Z1, Z1); // distances in '1', squares in '0'

                // dR2 = dX**2 + dY**2 + dZ**2
                Tmp0 = _mm256_add_pd(X0, Y0);
                Tmp0 = _mm256_add_pd(Tmp0, Z0);

                // tmpFac = G / (dR2 * sqrt(dR2))
                Tmp1 = _mm256_sqrt_pd(Tmp0);
                Tmp1 = _mm256_mul_pd(Tmp0, Tmp1);
                Tmp0 = _mm256_broadcast_sd(G);
                Tmp0 = _mm256_div_pd(Tmp0, Tmp1);

                X0 = _mm256_mul_pd(X1, Tmp0);
                Y0 = _mm256_mul_pd(Y1, Tmp0);
                Z0 = _mm256_mul_pd(Z1, Tmp0);

                /*
                 * First, compute the acceleration of object(s) 'j'
                 */
                Tmp0 = _mm256_broadcast_sd(mass + i + k);
                X1 = _mm256_mul_pd(X0, Tmp0);
                Y1 = _mm256_mul_pd(Y0, Tmp0);
                Z1 = _mm256_mul_pd(Z0, Tmp0);
                X2 = _mm256_load_pd(ax + j);
                Y2 = _mm256_load_pd(ay + j);
                Z2 = _mm256_load_pd(az + j);
                X2 = _mm256_sub_pd(X2, X1);
                Y2 = _mm256_sub_pd(Y2, Y1);
                Z2 = _mm256_sub_pd(Z2, Z1);
                _mm256_store_pd(ax + j, X2);
                _mm256_store_pd(ay + j, Y2);
                _mm256_store_pd(az + j, Z2);

                /*
                 * Now, accumulate the contributions to object 'i+k'
                 * in a dedicated register, that can be reduced
                 * at the end
                 */
                Tmp0 = _mm256_load_pd(mass + j);
                X1 = _mm256_mul_pd(X0, Tmp0);
                Y1 = _mm256_mul_pd(Y0, Tmp0);
                Z1 = _mm256_mul_pd(Z0, Tmp0);
                AccIX = _mm256_add_pd(AccIX, X1);
                AccIY = _mm256_add_pd(AccIY, Y1);
                AccIZ = _mm256_add_pd(AccIZ, Z1);
            }

            /*
             * Reduce the register(s) containing the contributions
             * to object 'i+k'
             */
            AccIX = _mm256_hadd_pd(AccIX, AccIX);
            AccIY = _mm256_hadd_pd(AccIY, AccIY);
            AccIZ = _mm256_hadd_pd(AccIZ, AccIZ);
            ax[i+k] += ((double*) &AccIX)[0] + ((double*) &AccIX)[2];
            ay[i+k] += ((double*) &AccIY)[0] + ((double*) &AccIY)[2];
            az[i+k] += ((double*) &AccIZ)[0] + ((double*) &AccIZ)[2];
        }
    }
}

void acc_co1_(int *restrict Nobj, double *restrict x, double *restrict a, double *restrict G, double *restrict mass)
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


void acc_co2_(int *restrict Nobj, double *restrict x, double *restrict a, double *restrict G, double *restrict mass)
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

