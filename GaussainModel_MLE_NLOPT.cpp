#include <iostream>
#include <vector>
#include <nlopt.hpp>
#include <cmath>
#include <fstream>

double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
    std::vector<double> *d = reinterpret_cast<std::vector<double>*>(my_func_data);

    double mu = x[0];
    double var = x[1];
    int data_size = d->size();

    if(grad) {
        grad[0] = 0.0;
        grad[1] = - data_size/(2.0*var);
        for(int i = 0; i < data_size; i++) {
            grad[0] += ((*d)[i] - mu)/var;
            grad[1] += ((*d)[i] - mu)*((*d)[i] - mu)/(2.0*var*var);
        }
    }
    //double func = data_size/sqrt(2*3.1415926*var);
    double func = -data_size*log(2*3.1415926*var)/2.0;
    for(int i = 0; i < data_size; i++) {
        func -= (((*d)[i] - mu)*((*d)[i] - mu))/(2*var);
    }
    printf("size = %d\n", data_size);
    // printf("mu = %f, var = %f\n", mu, var);
    // printf("grad[0] = %f, grad[1] = %f\n", grad[0], grad[1]);
    // printf("\n");
    return func;
}

 int main() { 

    std::vector<double> data;
    // Read data from csv file....
   // std::ifstream stockfile("hq_pure.csv");
    std::ifstream stockfile("sam.csv");
    std::string value;
    if(stockfile.good()) {
        while(getline(stockfile, value, '\n')) {
            data.push_back(std::stod(value));
        }
    }
    // use nlopt package to optimize the data...
    nlopt_opt opt;
    // NLOPT_LD_LBFGS Algorithms supprt unconstrained problems 
    opt = nlopt_create(NLOPT_LD_LBFGS, 2); /* algorithm and dimensionality */
    nlopt_set_max_objective(opt, myfunc, &data);
    nlopt_set_xtol_rel (opt, 1.0e-15);

    double x[2] = {0.5,2};  /* some initial guess */
    double minf; /* the minimum objective value, upon return */
    nlopt_optimize(opt, x, &minf);
    std::cout << "found maximum value at f(" << x[0] << ", " << x[1] << ") = " << minf << "\n" << std::endl;

    nlopt_destroy(opt);
    return 0;
 }
