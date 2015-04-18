#include "method.h"

double Solver::f(double x, double y, double z, const Param &param) {
    return -x * param["sigma"] + y * param["sigma"];
}

double Solver::g(double x, double y, double z, const Param &param) {
    return -x * z + x * param["r"] - y;
}

double Solver::h(double x, double y, double z, const Param &param) {
    return x * y - z * param["b"];
}

std::vector<std::vector<double>> ExplicitEuler::solve(const Param &param)
{
    double x = param["x0"];
    double y = param["y0"];
    double z = param["z0"];
    const double step = param["step"];

    std::vector<std::vector<double>> ans(3);

    ans[0].push_back(x);
    ans[1].push_back(y);
    ans[2].push_back(z);

    for(int i = 0; i < (int)param["steps"]; ++i)
    {
        double x1 = x + step * f(x, y ,z, param);
        double y1 = y + step * g(x, y, z, param);
        double z1 = z + step * h(x, y, z, param);

        z = z1;
        y = y1;
        x = x1;

        ans[0].push_back(x);
        ans[1].push_back(y);
        ans[2].push_back(z);
    }

    return ans;
}

std::vector<std::vector<double>> ImplicitEuler::solve(const Param &param) {
    double x = param["x0"];
    double y = param["y0"];
    double z = param["z0"];
    const double step = param["step"];

    std::vector<std::vector<double>> ans(3);

    ans[0].push_back(x);
    ans[1].push_back(y);
    ans[2].push_back(z);

    for(int i = 0; i < (int)param["steps"]; ++i)
    {
        double x1 = x + step * f(x, y ,z, param);
        double y1 = y + step * g(x, y, z, param);
        double z1 = z + step * h(x, y, z, param);

        double x2 = x + step / 2 * (f(x, y, z, param) + f(x1, y1, z1, param));
        double y2 = y + step / 2 * (g(x, y, z, param) + g(x1, y1, z1, param));
        double z2 = z + step / 2 * (h(x, y, z, param) + h(x1, y1, z1, param));

        z = z2;
        y = y2;
        x = x2;

        ans[0].push_back(x);
        ans[1].push_back(y);
        ans[2].push_back(z);
    }

    return ans;
}

std::vector<std::vector<double>> RungeKutta::solve(const Param &param)
{
    double x = param["x0"];
    double y = param["y0"];
    double z = param["z0"];
    const double step = param["step"];

    std::vector<std::vector<double>> ans(3, std::vector<double>((size_t)param["steps"]));

    for(int i = 0; i < (int)param["steps"]; ++i)
    {
        double k1 = f(x, y, z, param);
        double q1 = g(x, y, z, param);
        double p1 = h(x, y, z, param);

        double k2 = f(x + step / 2 * k1, y + step / 2 * q1, z + step / 2 * p1, param);
        double q2 = g(x + step / 2 * k1, y + step / 2 * q1, z + step / 2 * p1, param);
        double p2 = h(x + step / 2 * k1, y + step / 2 * q1, z + step / 2 * p1, param);

        double k3 = f(x + step / 2 * k2, y + step / 2 * q2, z + step / 2 * p2, param);
        double q3 = g(x + step / 2 * k2, y + step / 2 * q2, z + step / 2 * p2, param);
        double p3 = h(x + step / 2 * k2, y + step / 2 * q2, z + step / 2 * p2, param);

        double k4 = f(x + step * k3, y + step * q3, z + step * p3, param);
        double q4 = g(x + step * k3, y + step * q3, z + step * p3, param);
        double p4 = h(x + step * k3, y + step * q3, z + step * p3, param);

        double x1 = x + step / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
        double y1 = y + step / 6 * (q1 + 2 * q2 + 2 * q3 + q4);
        double z1 = z + step / 6 * (p1 + 2 * p2 + 2 * p3 + p4);

        z = z1;
        y = y1;
        x = x1;

        ans[0][i] = x;
        ans[1][i] = y;
        ans[2][i] = z;
    }

    return ans;
}

std::vector<std::vector<double>> Adams::solve(const Param &param) {
    const double step = param["step"];

    RungeKutta runge_kutta;

    auto ans = runge_kutta.solve(param);
    double x = ans[0][3];
    double y = ans[1][3];
    double z = ans[2][3];

    for(int i = 3; i < (int)param["steps"] - 1; ++i)
    {
        double x1 = x + step / 24 * (55 * f(ans[0][i], ans[1][i], ans[2][i], param) - 59 * f(ans[0][i - 1], ans[1][i - 1], ans[2][i - 1], param) + 37 * f(ans[0][i - 2], ans[1][i - 2], ans[2][i - 2], param) - 9 * f(ans[0][i - 3], ans[1][i - 3], ans[2][i - 3], param));
        double y1 = y + step / 24 * (55 * g(ans[0][i], ans[1][i], ans[2][i], param) - 59 * g(ans[0][i - 1], ans[1][i - 1], ans[2][i - 1], param) + 37 * g(ans[0][i - 2], ans[1][i - 2], ans[2][i - 2], param) - 9 * g(ans[0][i - 3], ans[1][i - 3], ans[2][i - 3], param));
        double z1 = z + step / 24 * (55 * h(ans[0][i], ans[1][i], ans[2][i], param) - 59 * h(ans[0][i - 1], ans[1][i - 1], ans[2][i - 1], param) + 37 * h(ans[0][i - 2], ans[1][i - 2], ans[2][i - 2], param) - 9 * h(ans[0][i - 3], ans[1][i - 3], ans[2][i - 3], param));

        z = z1;
        y = y1;
        x = x1;

        ans[0][i + 1] = x;
        ans[1][i + 1] = y;
        ans[2][i + 1] = z;
    }

    return ans;
}




