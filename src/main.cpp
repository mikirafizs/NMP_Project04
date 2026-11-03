#include <iostream>
#include <vector>
#include <cmath>
#include "sparse.hpp"


const size_t n = 500;
const double delta = 1;
const int itmax = 10000;
const  double epsilon = 1e-15;
double omega = 1;


int main(){
    std::cout<<"Enter omega value: ";
    std::cin>>omega;
    std::vector <double> rho = rho_population(n);
    std::vector <double> v (n, 1);

    size_t nnz_max = 3*n;

    std::vector <double> a (nnz_max);
    std::vector <size_t> col (nnz_max);
    std::vector <size_t> row (n+1);
    std::vector <double> d (n, 0.0);
    std::vector <double> r (n);

    for(std::size_t i = 0; i < n; i++){
        d[i] = -2/(delta*delta);
    }

    CSR_matrix(a, col, row, delta, n);


    double norm = 0.0;
    std::vector <double> y (n);
    int k = 0;
    do{
        k++;
        SOR(d, rho, v, omega, a, col, row);
        matvec(v, y, a, col, row, n);
        for(std::size_t i = 0; i < n; i++){
            r[i] = rho[i]-y[i];
        }
        norm = sqrt(vdot(r,r));
        if(k%1000==0)
        std::cout<<"iteration: "<<k<<" norm v: "<<sqrt(vdot(v,v))<<" norm r: "<<norm<<"\n";
    } while (k<itmax && norm > epsilon);
    
    return 0;
}