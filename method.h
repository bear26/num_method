#ifndef NUM_METHOD_METHOD_H
#define NUM_METHOD_METHOD_H

#include <vector>
#include <map>
#include <string>

struct Param
{
public:
    template <class T>
    inline void set(const std::string &key, T const &val) { _v_map[key] = std::to_string(val);}

    inline std::string get(const std::string &key) const { return _v_map.find(key)->second; }

    inline double operator[](const std::string &key) const { return std::stod(_v_map.find(key)->second); }

private:
    std::map<std::string, std::string> _v_map;
};

class Solver
{
public:
    virtual std::vector<std::vector<double>> solve(const Param &param) = 0;

protected:
    double f(double x, double y, double z, const Param &param);
    double g(double x, double y, double z, const Param &param);
    double h(double x, double y, double z, const Param &param);
};

class ExplicitEuler : public Solver
{
public:
    virtual std::vector<std::vector<double>> solve(const Param &param);
};

class ImplicitEuler : public Solver
{
public:
    virtual std::vector<std::vector<double>> solve(const Param &param);
};

class RungeKutta : public Solver
{
public:
    virtual std::vector<std::vector<double>> solve(const Param &param);
};

class Adams : public Solver
{
public:
    virtual std::vector<std::vector<double>> solve(const Param &param);
};

#endif //NUM_METHOD_METHOD_H
