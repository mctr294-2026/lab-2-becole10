#include <functional>
#include <cmath>

#include "roots.hpp"

bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root)

{
    double tol = 1e-6; // tolerance
    int itns = std::ceil(std::log2((std::abs(a - b)) / tol)); // number of iterations
    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0)
    {
        return false; // no crossing between a and b
    }

    for (int i = 0; i <= itns + 1; i++) // loops itns+1 times
    {

        double mid = (a + b) / 2; // midpoint between a and b
        double fmid = f(mid);

        if (std::abs(b-a) < tol || std::abs(fmid) < tol)
        {
            *root = mid;
            return true;
        }

        if (fa * fmid < 0)
        {
            b = mid;
            fb = f(b);
        }
        else
        {
            a = mid;
            fa = f(a);
        }
    }
    return false;
}

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root)

{
    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0)
    {
        return false; // no crossing between a and b
    }

    double tol = 1e-6; // tolerance
    int itns = 1000000;

    for (int i = 0; i <= itns; i++) // loops itns+1 times
    {

        double mid = (a*fb - b*fa) / (fb - fa); // "midpoint" between a and b
        double fmid = f(mid);

        if (std::abs(b-a) < tol || std::abs(fmid) < tol)
        {
            *root = mid;
            return true;
        }

        if (fa * fmid < 0)
        {
            b = mid;
            fb = f(b);
        }
        else
        {
            a = mid;
            fa = f(a);
        }
    }
    return false;
}

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root)

{
    double tol = 1e-6; // tolerance
    double c2 = c+1; // placeholder
    int itns = 1000; // iterations of for loop

    for (int i = 0; i <= itns; i++)
    {
        if (g(c) == 0)
        {
            return false; // Divide by zero
        }

        c2 = c - f(c)/g(c);

        if (c2 < a || c2 > b)
        {
            return false; // outside interval
        }

        if (std::abs(c2-c) < tol)
        {
            *root = c2;
            return true;
        }
        else
        {
            c = c2;
        }
    }
    return false;
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root)

{
    double tol = 1e-6; // tolerance
    double c2 = c+1; // placeholder
    int itns = 1000; // iterations of for loop

    for (int i = 0; i <= itns; i++)
    {
        if (f(c2)-f(c) == 0)
        {
            return false; // divide by zero
        }

        double c3 = c2 - f(c2)*((c2-c)/(f(c2)-f(c)));

        if (c3 < a || c3 > b)
        {
            return false; // outside interval
        }

        if (std::abs(c3-c2) < tol)
        {
            *root = c3;
            return true;
        }
        else
        {
            c = c2;
            c2 = c3;
        }
    }
    return false;
}