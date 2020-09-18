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
 * along with molconv. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef ACCEL_H
#define ACCEL_H


void acc_1c_(int *restrict Nobj, double *restrict x, double *restrict a, const double *restrict G, double *restrict mass);
void acc_2c_(int *restrict Nobj, double *restrict x, double *restrict a, const double *restrict G, double *restrict mass);
void acc_2c2_(int *restrict Nobj,
              double *restrict xx, double *restrict xy, double *restrict xz,
              double *restrict ax, double *restrict ay, double *restrict az,
              const double *restrict G, double *restrict mass);
void acc_2c2avx_(int *restrict Nobj,
                 double *restrict xx, double *restrict xy, double *restrict xz,
                 double *restrict ax, double *restrict ay, double *restrict az,
                 const double *restrict G, double *restrict mass);
void acc_co1_(int *restrict Nobj, double *restrict x, double *restrict a, const double *restrict G, double *restrict mass);
void acc_co2_(int *restrict Nobj, double *restrict x, double *restrict a, const double *restrict G, double *restrict mass);

#endif // ACCEL_H
