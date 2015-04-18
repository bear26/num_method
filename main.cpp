#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>

#include "method.h"

#include <boost/program_options.hpp>

std::string path;

void write(const std::vector<std::vector<double>> &ans, const Param &param)
{
    const double step = param["step"];

    FILE *f_x = fopen((path + ".x(t)").data(), "w");
    FILE *f_y = fopen((path + ".y(t)").data(), "w");
    FILE *f_z = fopen((path + ".z(t)").data(), "w");
    FILE *f_xyz = fopen((path + ".xyz").data(), "w");

    for(size_t i = 0; i < ans[0].size(); ++i)
    {
        fprintf(f_x, "%.9lf %.9lf\n", i * step, ans[0][i]);
        fprintf(f_y, "%.9lf %.9lf\n", i * step, ans[1][i]);
        fprintf(f_z, "%.9lf %.9lf\n", i * step, ans[2][i]);
        fprintf(f_xyz, "%.9lf %.9lf %.9lf\n", ans[0][i], ans[1][i], ans[2][i]);
    }

    fclose(f_x);
    fclose(f_y);
    fclose(f_z);
}

//! возвращает номер метода
//! 0 - явный эйлера, 1 - неявный эйлера, 2 - рунге - кутты, 3 - адамса
int parse(int argc, char *argv[], Param &param)
{
    double sigma;
    double b;
    double r;
    double step;
    int n;
    double x0, y0, z0;
    int num_m;
    boost::program_options::options_description options;
    options.add_options()("help", "produce help message")
            ("sigma,S",boost::program_options::value<double>(&sigma)->default_value(10), "sigma value")
            ("b,B", boost::program_options::value<double>(&b)->default_value(8.0 / 3), "b")
            ("r,R", boost::program_options::value<double>(&r)->required(), "r")
            ("step,s", boost::program_options::value<double>(&step)->required(), "step length")
            ("count_steps,n", boost::program_options::value<int>(&n)->required(), "count steps")
            ("x0,x", boost::program_options::value<double>(&x0)->required(), "x0")
            ("y0,y", boost::program_options::value<double>(&y0)->required(), "y0")
            ("z0,z", boost::program_options::value<double>(&z0)->required(), "z0")
            ("path,p", boost::program_options::value<std::string>(&path)->required(), "path to output file")
            ("num_method,m", boost::program_options::value<int>(&num_m)->required(), "0 - explicit Euler\n1 - implicit Euler\n2 - Runge Kutta\n3 - explicit Adams");

    try {
        boost::program_options::variables_map vm;
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), vm);
        boost::program_options::notify(vm);
    }
    catch(boost::program_options::error)
    {
        std::cout << options << std::endl;
        exit(1);
    }

    param.set("sigma", sigma);
    param.set("b", b);
    param.set("r", r);
    param.set("step", step);
    param.set("steps", n);
    param.set("x0", x0);
    param.set("y0", y0);
    param.set("z0", z0);

    return num_m;
}

int main(int argc, char *argv[]) {

    std::vector<Solver*> methods;
    ExplicitEuler m1;
    ImplicitEuler m2;
    RungeKutta m3;
    Adams m4;

    methods.push_back(&m1);
    methods.push_back(&m2);
    methods.push_back(&m3);
    methods.push_back(&m4);

    Param param;
    write(methods[parse(argc, argv, param)]->solve(param), param);
    return 0;
}