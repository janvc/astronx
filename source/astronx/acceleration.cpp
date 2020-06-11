#include "acceleration.h"


namespace Astronx {

void Acceleration::acceleration(double *__restrict__ x, double *__restrict__ a)
{

}

void Acceleration::acc_ser(double *__restrict__ x, double *__restrict__ a)
{
    int i, j;

    for (j = 0; j < m_Nobj; j++)
    {
        a[0 * m_Npad + j] = 0.0;
        a[1 * m_Npad + j] = 0.0;
        a[2 * m_Npad + j] = 0.0;
    }

    for (i = 0; i < m_Nobj - 1; i++)
    {
        for (j = i + 1; j < m_Nobj; j++)
        {
            double dX = x[0 * m_Npad + j] - x[0 * m_Npad + i];
            double dY = x[1 * m_Npad + j] - x[1 * m_Npad + i];
            double dZ = x[2 * m_Npad + j] - x[2 * m_Npad + i];
            double R2 = dX * dX + dY * dY + dZ * dZ;
            double tmpFac = PhyCon::G / (R2 * std::sqrt(R2));
            a[0 * m_Npad + i] += m_masses[j] * tmpFac * dX;
            a[1 * m_Npad + i] += m_masses[j] * tmpFac * dY;
            a[2 * m_Npad + i] += m_masses[j] * tmpFac * dZ;
            a[0 * m_Npad + j] -= m_masses[i] * tmpFac * dX;
            a[1 * m_Npad + j] -= m_masses[i] * tmpFac * dY;
            a[2 * m_Npad + j] -= m_masses[i] * tmpFac * dZ;
        }
    }
}

/*
 * three arrays, avx intrinsics, blocked loop
 */
void Acceleration::acc_avx_block_nok_(const int Nobj, const int Npad,
                        const double *restrict xx, const double *restrict xy, const double *restrict xz,
                        double *restrict ax, double *restrict ay, double *restrict az,
                        const double *restrict mass)
{
    int i, j, l;

    for (j = 0; j < Npad; j++)
    {
        ax[j] = 0.0;
        ay[j] = 0.0;
        az[j] = 0.0;
    }

#define BLOCK_SIZE 10000

    for (l = 0; l < Nobj; l += BLOCK_SIZE)
    {
        int ijmax = l + BLOCK_SIZE > Nobj ? Nobj : l + BLOCK_SIZE;

        for (i = 0; i < ijmax; i++)
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
            if (i >= l)
            {
                double dX, dY, dZ, R2, tmpFac;
                __m128d x0, y0, z0, x1, y1, z1, x2, y2, z2;
                __m128d tmp0, tmp1;

                if (i % 4 == 0)
                {
                    // interaction between objects 0 and 1
                    dX = xx[i+1] - xx[i];
                    dY = xy[i+1] - xy[i];
                    dZ = xz[i+1] - xz[i];
                    R2 = dX * dX + dY * dY + dZ * dZ;
                    tmpFac = G / (R2 * sqrt(R2));
                    ax[i] += mass[i+1] * tmpFac * dX;
                    ay[i] += mass[i+1] * tmpFac * dY;
                    az[i] += mass[i+1] * tmpFac * dZ;
                    ax[i+1] -= mass[i] * tmpFac * dX;
                    ay[i+1] -= mass[i] * tmpFac * dY;
                    az[i+1] -= mass[i] * tmpFac * dZ;

                    // interaction between objects 0 and 2/3
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
                    tmp0 = _mm_set_pd1(G);
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
                else if (i % 4 == 1)
                {
                    // interaction between objects 1 and 2/3
                    x0 = _mm_load_pd1(xx + i);  // load the position of object 'i+1' into '0'
                    y0 = _mm_load_pd1(xy + i);
                    z0 = _mm_load_pd1(xz + i);
                    x1 = _mm_load_pd(xx + i + 1);   // load objects 'i+2' and 'i+3' into '1'
                    y1 = _mm_load_pd(xy + i + 1);
                    z1 = _mm_load_pd(xz + i + 1);

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
                    tmp0 = _mm_set_pd1(G);
                    tmp0 = _mm_div_pd(tmp0, tmp1);  // 'tmp0' contains the conversion factor

                    x0 = _mm_mul_pd(x1, tmp0);  // put 'factor * dx' into '0'
                    y0 = _mm_mul_pd(y1, tmp0);
                    z0 = _mm_mul_pd(z1, tmp0);

                    tmp0 = _mm_load_pd1(mass + i);      // mass of 'i+1' into 'tmp0'
                    tmp1 = _mm_load_pd(mass + i + 1);   // massES of 'i+2' and 'i+3' into 'tmp1'

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
                    x2 = _mm_load_pd(ax + i + 1);   // load acceleration of other objects
                    y2 = _mm_load_pd(ay + i + 1);
                    z2 = _mm_load_pd(az + i + 1);
                    x2 = _mm_sub_pd(x2, x1);        // subtract from them
                    y2 = _mm_sub_pd(y2, y1);
                    z2 = _mm_sub_pd(z2, z1);
                    _mm_store_pd(ax + i + 1, x2);   // and store again
                    _mm_store_pd(ay + i + 1, y2);
                    _mm_store_pd(az + i + 1, z2);
                }
                else if (i % 4 == 2)
                {
                    // interaction between objects 2 and 3
                    dX = xx[i+1] - xx[i];
                    dY = xy[i+1] - xy[i];
                    dZ = xz[i+1] - xz[i];
                    R2 = dX * dX + dY * dY + dZ * dZ;
                    tmpFac = G / (R2 * sqrt(R2));
                    ax[i] += mass[i+1] * tmpFac * dX;
                    ay[i] += mass[i+1] * tmpFac * dY;
                    az[i] += mass[i+1] * tmpFac * dZ;
                    ax[i+1] -= mass[i] * tmpFac * dX;
                    ay[i+1] -= mass[i] * tmpFac * dY;
                    az[i+1] -= mass[i] * tmpFac * dZ;
                }
            }

            /*
             * Now do the remaining part of the row with AVX functions
             */
            __m256d X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2;
            __m256d Tmp0, Tmp1;
            __m256d AccIX, AccIY, AccIZ;

            AccIX = _mm256_setzero_pd();
            AccIY = _mm256_setzero_pd();
            AccIZ = _mm256_setzero_pd();

            int jstart = i >= l ? i - (i % 4) + 4 : l;

            for (j = jstart; j < ijmax; j += 4)
            {
                // put the position of object "i + k" into all
                // four slots of the AVX registers
                X0 = _mm256_broadcast_sd(xx + i);
                Y0 = _mm256_broadcast_sd(xy + i);
                Z0 = _mm256_broadcast_sd(xz + i);

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
                Tmp0 = _mm256_set1_pd(G);
                Tmp0 = _mm256_div_pd(Tmp0, Tmp1);

                X0 = _mm256_mul_pd(X1, Tmp0);
                Y0 = _mm256_mul_pd(Y1, Tmp0);
                Z0 = _mm256_mul_pd(Z1, Tmp0);

                /*
                 * First, compute the acceleration of object(s) 'j'
                 */
                Tmp0 = _mm256_broadcast_sd(mass + i);
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
            } // j

            /*
             * Reduce the register(s) containing the contributions
             * to object 'i+k'
             */
            AccIX = _mm256_hadd_pd(AccIX, AccIX);
            AccIY = _mm256_hadd_pd(AccIY, AccIY);
            AccIZ = _mm256_hadd_pd(AccIZ, AccIZ);
            ax[i] += ((double*) &AccIX)[0] + ((double*) &AccIX)[2];
            ay[i] += ((double*) &AccIY)[0] + ((double*) &AccIY)[2];
            az[i] += ((double*) &AccIZ)[0] + ((double*) &AccIZ)[2];
        } // i
    } // l
}

}
