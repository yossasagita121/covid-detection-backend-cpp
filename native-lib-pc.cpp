// #include <jni.h>
#include <string>
#include "detection.h"
#include "functions.h"
// extern "C" JNIEXPORT jstring JNICALL

using namespace std;

std::string native_lib_main(double data[DIMENSION][SAMPLE_MAX])
{
    float Threshold = 1.15, s = 0;
    std::string report;
    singularity sing;

    int i = 0, mode = 1;
    std::string hello = "Hello from C++";
    std::string positif = "Anda Terdeteksi Positif";
    std::string negatif = "Anda Terdeteksi Negatif";
    std::string ragu = "Tidak Konklusif";

    report = hello + "\n";

    get_data(data);
    proses();

    // std::string report;
    for (i = 0; i < EMBED; i++)
    {

        report = report + "\t" + std::to_string(i + 1) + "\t" + std::to_string(localdimension[i]) + "\t" + std::to_string(localsize[i]) + "\n";
    }

    s = measure(sing);

    report = report + "Dim: " + std::to_string(sing.dim) + ", Size: " + std::to_string(sing.size) + "\t" + " Dispersi: " + std::to_string(sing.dispersi) + "\n";

    if (sing.dispersi > 0.1)
    {
        report = report + ragu;
    }
    else
    {
        if (s > Threshold)
        {
            report = report + positif;
        }
        else
        {
            report = report + negatif;
        }
    }
    return report;
}

float measure(singularity &s)
{
    float dim = 0, size = 0, sum = 0, variansi = 0;
    int i;

    for (i = 0; i < EMBED; i++)
    {
        if (i > 2)
        {
            dim = dim + localdimension[i];
            size = size + localsize[i];
            sum = sum + 1.0;
        }
    }

    if (sum != 0)
    {
        s.dim = dim / sum;
        s.size = size / sum;
    }
    else
    {
        s.dim = 0;
        s.size = 0;
    }

    sum = 0.0;

    for (i = 0; i < EMBED; i++)
    {
        if (i > 2)
        {
            variansi = variansi + (localdimension[i] - s.dim) * (localdimension[i] - s.dim);

            sum = sum + 1.0;
        }
    }

    if (s.dim != 0)
    {
        s.dispersi = variansi / s.dim;
    }
    else
    {
        s.dispersi = 0;
    }

    return s.dim;
};

int proses(void)
{
    char dimset = 0, rescale_set = 0, eps_min_set = 0, eps_max_set = 0;
    char smaller;
    long lnorm;
    int maxembed;
    long i1, j1, x, y, sn, n, i, j, n1, n2;
    double eps, EPSMAX1, maxinterval;
    time_t mytime, lasttime;

    double min, interval;
    if (rescale_set)
    {
        for (i = 0; i < DIM; i++)
            rescale_data(series[i], SAMPLE_MAX, &min, &interval);
        maxinterval = 1.0;
    }
    else
    {
        maxinterval = 0.0;
        for (i = 0; i < DIM; i++)
        {
            min = interval = series[i][0];
            for (j = 1; j < SAMPLE_MAX; j++)
            {
                if (min > series[i][j])
                    min = series[i][j];
                if (interval < series[i][j])
                    interval = series[i][j];
            }
            interval -= min;
            if (interval > maxinterval)
                maxinterval = interval;
        }
    }
    if (!eps_max_set)
        EPSMAX *= maxinterval;
    if (!eps_min_set)
        EPSMIN *= maxinterval;
    EPSMAX = (fabs(EPSMAX) < maxinterval) ? fabs(EPSMAX) : maxinterval;
    EPSMIN = (fabs(EPSMIN) < EPSMAX) ? fabs(EPSMIN) : EPSMAX / 2.;
    EPSMAX1 = EPSMAX;

    howoften1 = HOWOFTEN - 1;
    maxembed = DIM * EMBED - 1;

    epsinv = 1.0 / EPSMAX;
    epsfactor = pow(EPSMAX / EPSMIN, 1.0 / (double)howoften1);
    lneps = log(EPSMAX);
    lnfac = log(epsfactor);

    epsm[0] = EPSMAX;
    norm[0] = 0.0;
    for (i = 1; i < HOWOFTEN; i++)
    {
        norm[i] = 0.0;
        epsm[i] = epsm[i - 1] / epsfactor;
    }
    imin = 0;

    scramble();
    for (i = 0; i < (SAMPLE_MAX - (EMBED - 1) * DELAY); i++)
        oscr[scr[i]] = i;

    for (i = 0; i < DIM * EMBED; i++)
        for (j = 0; j < HOWOFTEN; j++)
            found[i][j] = 0.0;

    nmax = SAMPLE_MAX - DELAY * (EMBED - 1);

    for (i1 = 0; i1 < NMAX; i1++)
    {
        boxc1[i1] = -1;
        for (j1 = 0; j1 < NMAX; j1++)
            box[i1][j1] = -1;
    }
    time(&lasttime);
    lnorm = 0;
    // return 1;
    for (n = 1; n < nmax; n++)
    {
        smaller = 0;
        sn = scr[n - 1];
        if (DIM > 1)
        {
            x = (long)(series[0][sn] * epsinv) & imax;
            y = (long)(series[1][sn] * epsinv) & imax;
        }
        else
        {
            x = (long)(series[0][sn] * epsinv) & imax;
            y = (long)(series[0][sn + DELAY] * epsinv) & imax;
        }
        list[sn] = box[x][y];
        box[x][y] = sn;
        listc1[sn] = boxc1[x];
        boxc1[x] = sn;

        i = imin;
        while (found[maxembed][i] >= MAXFOUND)
        {
            smaller = 1;
            if (++i > howoften1)
                break;
        }
        if (smaller)
        {
            imin = i;
            if (imin <= howoften1)
            {
                EPSMAX = epsm[imin];
                epsinv = 1.0 / EPSMAX;
                for (i1 = 0; i1 < NMAX; i1++)
                {
                    boxc1[i1] = -1;
                    for (j1 = 0; j1 < NMAX; j1++)
                        box[i1][j1] = -1;
                }
                for (i1 = 0; i1 < n; i1++)
                {
                    sn = scr[i1];
                    if (DIM > 1)
                    {
                        x = (long)(series[0][sn] * epsinv) & imax;
                        y = (long)(series[1][sn] * epsinv) & imax;
                    }
                    else
                    {
                        x = (long)(series[0][sn] * epsinv) & imax;
                        y = (long)(series[0][sn + DELAY] * epsinv) & imax;
                    }
                    list[sn] = box[x][y];
                    box[x][y] = sn;
                    listc1[sn] = boxc1[x];
                    boxc1[x] = sn;
                }
            }
        }

        if (imin <= howoften1)
        {
            lnorm = n;
            if (MINDIST > 0)
            {
                sn = scr[n];
                n1 = (sn - (long)MINDIST >= 0) ? sn - (long)MINDIST : 0;
                n2 = (sn + MINDIST < SAMPLE_MAX - (EMBED - 1) * DELAY) ? sn + MINDIST : SAMPLE_MAX - (EMBED - 1) * DELAY - 1;
                for (i1 = n1; i1 <= n2; i1++)
                    if ((oscr[i1] < n))
                        lnorm--;
            }

            if (EMBED * DIM > 1)
                make_c2_dim(n);
            make_c2_1(n);
            for (i = imin; i < HOWOFTEN; i++)
                norm[i] += (double)(lnorm);
        }
        // fprintf(stderr,"outer loop\n");
        if (
            ((time(&mytime) - lasttime) > WHEN) ||
            n == (nmax - 1) || (imin > howoften1))
        {
            time(&lasttime);
            for (i = 0; i < DIM * EMBED; i++)
            {

                eps = EPSMAX1;
                localdimension[i] = 0;
                validscale = 0;
                for (j = 1; j < HOWOFTEN; j++)
                {
                    eps /= epsfactor;
                    if ((found[i][j] > 0.0) && (found[i][j - 1] > 0.0))
                    {
                        if (eps > 0.01 && eps < 0.1)
                        {
                            localdimension[i] += log(found[i][j - 1] / found[i][j] / norm[j - 1] * norm[j]) / lnfac;
                            validscale++;
                        }
                    }
                }
                localdimension[i] = localdimension[i] / validscale;
            }

            for (i = 0; i < EMBED * DIM; i++)
            {
                eps = EPSMAX1 * epsfactor;
                localsize[i] = 0;
                validscale = 0;
                for (j = 0; j < HOWOFTEN; j++)
                {
                    eps /= epsfactor;
                    if (norm[j] > 0.0)
                    {
                        if ((found[i][j] > 0.0) && (found[i][j - 1] > 0.0))
                        {

                            if (eps > 0.01 && eps < 0.1)
                            {
                                localsize[i] += found[i][j] / norm[j] * pow(eps, -localdimension[i]);
                                validscale++;
                            }
                        }
                    }
                }
                localsize[i] = localsize[i] / validscale;
                // fprintf(stderr,"\t%ld\t%f\t%f\n",i+1, localdimension[i], localsize[i]);
            }
            // fprintf(stderr,"\n");

            if (imin > howoften1)
                return (0);
        }
    }
    return 0;
}

void scramble(void)
{
    long i, j, k, m;
    unsigned long rnd, rndf, hlength, allscr = 0;
    long scfound[SAMPLE_MAX - (EMBED - 1) * DELAY], scnhelp[SAMPLE_MAX - (EMBED - 1) * DELAY], scnfound;
    long scbox[SCBOX], lswap, element, scbox1 = SCBOX - 1;
    double rz[SAMPLE_MAX - (EMBED - 1) * DELAY];
    double schelp[SAMPLE_MAX - (EMBED - 1) * DELAY];
    double sceps = (double)SCBOX - 0.001, swap;

    hlength = SAMPLE_MAX - (EMBED - 1) * DELAY;

    if (sizeof(long) == 8)
    {
        rndf = 13 * 13 * 13 * 13;
        rndf = rndf * rndf * rndf * 13;
        rnd = 0x849178L;
    }
    else
    {
        rndf = 69069;
        rnd = 0x234571L;
    }
    for (i = 0; i < 1000; i++)
        rnd = rnd * rndf + 1;

    for (i = 0; i < hlength; i++)
        rz[i] = (double)(rnd = rnd * rndf + 1) / ULONG_MAX;

    for (i = 0; i < SCBOX; i++)
        scbox[i] = -1;
    for (i = 0; i < hlength; i++)
    {
        m = (int)(rz[i] * sceps) & scbox1;
        scfound[i] = scbox[m];
        scbox[m] = i;
    }
    for (i = 0; i < SCBOX; i++)
    {
        scnfound = 0;
        element = scbox[i];
        while (element != -1)
        {
            scnhelp[scnfound] = element;
            schelp[scnfound++] = rz[element];
            element = scfound[element];
        }

        for (j = 0; j < scnfound - 1; j++)
            for (k = j + 1; k < scnfound; k++)
                if (schelp[k] < schelp[j])
                {
                    swap = schelp[k];
                    schelp[k] = schelp[j];
                    schelp[j] = swap;
                    lswap = scnhelp[k];
                    scnhelp[k] = scnhelp[j];
                    scnhelp[j] = lswap;
                }
        for (j = 0; j < scnfound; j++)
            scr[allscr + j] = scnhelp[j];
        allscr += scnfound;
    }
}

void make_c2_dim(int n)
{
    char small;
    long i, j, k, x, y, i1, i2, j1, element, n1, maxi, count, hi;
    double hs[EMBED * DIM], max, dx;

    n1 = scr[n];

    count = 0;
    for (i1 = 0; i1 < EMBED; i1++)
    {
        i2 = i1 * DELAY;
        for (j = 0; j < DIM; j++)
            hs[count++] = series[j][n1 + i2];
    }

    x = (int)(hs[0] * epsinv) & imax;
    y = (int)(hs[1] * epsinv) & imax;

    for (i1 = x - 1; i1 <= x + 1; i1++)
    {
        i2 = i1 & imax;
        for (j1 = y - 1; j1 <= y + 1; j1++)
        {
            element = box[i2][j1 & imax];
            while (element != -1)
            {
                if (labs((long)(element - n1)) > MINDIST)
                {
                    count = 0;
                    max = 0.0;
                    maxi = howoften1;
                    small = 0;
                    for (i = 0; i < EMBED; i++)
                    {
                        hi = i * DELAY;
                        for (j = 0; j < DIM; j++)
                        {
                            dx = fabs(hs[count] - series[j][element + hi]);
                            if (dx <= EPSMAX)
                            {
                                if (dx > max)
                                {
                                    max = dx;
                                    if (max < EPSMIN)
                                    {
                                        maxi = howoften1;
                                    }
                                    else
                                    {
                                        maxi = (lneps - log(max)) / lnfac;
                                    }
                                }
                                if (count > 0)
                                    for (k = imin; k <= maxi; k++)
                                        found[count][k] += 1.0;
                            }
                            else
                            {
                                small = 1;
                                break;
                            }
                            count++;
                        }
                        if (small)
                            break;
                    }
                }
                element = list[element];
            }
        }
    }
}

void make_c2_1(int n)
{
    int i, x, i1, maxi;
    long element, n1;
    double hs, max;
    // unsigned long MINDIST=0,MAXFOUND=1000;
    n1 = scr[n];
    hs = series[0][n1];

    x = (int)(hs * epsinv) & imax;

    for (i1 = x - 1; i1 <= x + 1; i1++)
    {
        element = boxc1[i1 & imax];
        while (element != -1)
        {
            if (labs(element - n1) > MINDIST)
            {
                max = fabs(hs - series[0][element]);
                if (max <= EPSMAX)
                {
                    if (max < EPSMIN)
                        maxi = howoften1;
                    else
                        maxi = (lneps - log(max)) / lnfac;
                    for (i = imin; i <= maxi; i++)
                        found[0][i] += 1.0;
                }
            }
            element = listc1[element];
        }
    }
}

void rescale_data(double *x, unsigned long l, double *min, double *interval)
{
    int i;

    *min = *interval = x[0];

    for (i = 1; i < l; i++)
    {
        if (x[i] < *min)
            *min = x[i];
        if (x[i] > *interval)
            *interval = x[i];
    }
    *interval -= *min;

    if (*interval != 0.0)
    {
        for (i = 0; i < l; i++)
            x[i] = (x[i] - *min) / *interval;
    }
    else
    {
        fprintf(stderr, "rescale_data: data ranges from %e to %e. It makes\n"
                        "\t\tno sense to continue. Exiting!\n\n",
                *min, *min + (*interval));
        exit(RESCALE_DATA_ZERO_INTERVAL);
    }
}
