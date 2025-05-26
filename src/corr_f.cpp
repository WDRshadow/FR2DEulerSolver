#include "corr_f.h"

double g2(const double s)
{
    return -1.0 / 4.0 * (1 - s * s);
}

double dg2(const double s)
{
    return 0.5 * s;
}

double g2L(const double s)
{
    return 1.0 / 2.0 * (1.0 - s) * (1.0 + g2(s));
}

double g2R(const double s)
{
    return 1.0 / 2.0 * (1.0 + s) * (1.0 - g2(s));
}

double dg2L(const double s)
{
    return 1.0 / 2.0 * (-1 - g2(s) + dg2(s) - s * dg2(s));
}

double dg2R(const double s)
{
    return 1.0 / 2.0 * (1 - g2(s) - dg2(s) - s * dg2(s));
}
