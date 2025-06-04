#include "corr_f.h"

double g2L(const double s)
{
    return 0.25 * (-5 * s * s * s + 3 * s * s + 3 * s - 1);
}

double g2R(const double s)
{
    return 0.25 * (5 * s * s * s + 3 * s * s - 3 * s - 1);
}

double dg2L(const double s)
{
    return 0.25 * (-15 * s * s + 6 * s + 3);
}

double dg2R(const double s)
{
    return 0.25 * (15 * s * s + 6 * s - 3);
}
