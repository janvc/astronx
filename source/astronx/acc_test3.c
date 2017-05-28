#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char *argv[])
{
    int Nobj = atoi(argv[1]);
    int Niter = atoi(argv[2]);
    const double G = 6.6726e-11;

    double *xx;
    double *xy;
    double *xz;
    double *ax;
    double *ay;
    double *az;
    double *mass;

    void *dummy1;
    void *dummy2;
    void *dummy3;
    void *dummy4;
    void *dummy5;
    void *dummy6;
    void *dummy7;
    posix_memalign(&dummy1, 32, Nobj * sizeof(double));
    posix_memalign(&dummy2, 32, Nobj * sizeof(double));
    posix_memalign(&dummy3, 32, Nobj * sizeof(double));
    posix_memalign(&dummy4, 32, Nobj * sizeof(double));
    posix_memalign(&dummy5, 32, Nobj * sizeof(double));
    posix_memalign(&dummy6, 32, Nobj * sizeof(double));
    posix_memalign(&dummy7, 32, Nobj * sizeof(double));
    xx = (double*) dummy1;
    xy = (double*) dummy2;
    xz = (double*) dummy3;
    ax = (double*) dummy4;
    ay = (double*) dummy5;
    az = (double*) dummy6;
    mass = (double*) dummy7;

    srand(Nobj);
    for (int i = 0; i < Nobj; i++)
    {
        mass[i] = (rand() / (RAND_MAX / 1.0e25));
        xx[i] = -1.0e10 + (rand() / (RAND_MAX / 2.0e10));
        xy[i] = -1.0e10 + (rand() / (RAND_MAX / 2.0e10));
        xz[i] = -1.0e10 + (rand() / (RAND_MAX / 2.0e10));
    }

    clock_t t0, t1;
    t0 = clock();
    for (int i = 0; i < Niter; i++)
        acc_2c2avx_(&Nobj, xx, xy, xz, ax, ay, az, &G, mass);
    t1 = clock();
    double cpu_time_used = ((double) (t1 - t0)) / CLOCKS_PER_SEC;

    printf("%d %d %f %f %f %f\n", Nobj, Niter, ax[0], ay[0], az[0], cpu_time_used);

    if (argc == 4)
    {
        for (int i = 0; i < Nobj; i++)
            printf("%7d %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e \n", i+1, mass[i], xx[i], xy[i], xz[i], ax[i], ay[i], az[i]);
        for (int i = 0; i < Nobj; i++)
            printf("%22.15e\n", mass[i]);
        for (int i = 0; i < Nobj; i++)
            printf("%22.15e\n", xx[i]);
        for (int i = 0; i < Nobj; i++)
            printf("%22.15e\n", xy[i]);
        for (int i = 0; i < Nobj; i++)
            printf("%22.15e\n", xz[i]);
        for (int i = 0; i < Nobj; i++)
            printf("%22.15e\n", ax[i]);
        for (int i = 0; i < Nobj; i++)
            printf("%22.15e\n", ay[i]);
        for (int i = 0; i < Nobj; i++)
            printf("%22.15e\n", az[i]);
    }

    return 0;
}
