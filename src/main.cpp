#include <iostream>
#include <vector>
#include <cmath>
#include "sparse.hpp"
#include <fstream>
#include <sstream>
#include <string>


const size_t n = 500;
const double delta = 1;
const int itmax = 10000;
const  double epsilon = 1e-15;
double omega;


int main(){
    std::cout<<"Enter omega value: ";
    std::cin>>omega;
    std::vector <double> rho = rho_population(n);
    std::vector <double> v (n, 0.0);
    /*for (size_t i = 0; i < n; i++)
    {
        v[i] = 100*(double)rand()/RAND_MAX;
    }*/
    

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
    std::ostringstream name;
    name << "../data/sor_omega_"<<omega<<".csv";
    std::ofstream solution(name.str());
    if(!solution.is_open()){
        std::cerr << "Cannot open " << name.str() << "\n";
        return 1;
    }
    
    /*std::ofstream solution;
    solution.open("../data/random_vec_sol.csv");
    if(!solution.is_open()){
        std::cerr << "Cannot open \n";
        return 1;
    }*/
    do{
        k++;
        SOR(d, rho, v, omega, a, col, row);
        matvec(v, y, a, col, row, n);
        for(std::size_t i = 0; i < n; i++){
            r[i] = rho[i]-y[i];
        }
        norm = sqrt(vdot(r,r));
        
        solution << k << "," << sqrt(vdot(v,v)) << "," << norm << "\n";
    } while (k<itmax && norm > epsilon);
    
    std::ostringstream name2;
    name2 << "../data/vector_omega_" << omega << ".csv";
    std::ofstream vec(name2.str());
    if(!vec.is_open()){
        std::cerr << "Cannot open " << name2.str() << "\n";
        return 1;
    }
    for(size_t i = 0; i < n; i++)
    vec << i*delta << "," << rho[i] << "," << v[i] << "\n";

    return 0;
}