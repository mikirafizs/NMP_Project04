#ifndef SPARSE_HPP
#define SPARSE_HPP
#include <cmath>
#include <vector>
#include <iostream>

inline void SOR(std::vector <double> &d, std::vector <double> &rho, std::vector <double> &v, double omega, const std::vector <double> &a, const std::vector <std::size_t> &col, const std::vector <std::size_t> &row){
    for(size_t i = 0; i < rho.size(); i++){
        double c = 0;
        for(size_t l = row[i]; l < row[i+1]; l++){
            size_t j = col[l];
            c+=a[l]*v[j];
        }
        v[i] = v[i] + (omega/d[i])*(rho[i] - c);
    }
}

inline std::vector <double> rho_population(size_t n){
    std::vector <double> rho (n, 0);
    const double PI = 3.14159265358979323846;
    for(std::size_t i = 0; i < n; i++){
        rho[i] = (1/10)*pow(sin((2*PI*i)/(n-1)), 10);
    }
    return rho;
}

inline double vdot(std::vector <double> a, std::vector <double> b){
    if(a.size() != b.size()){
        std::cerr << "Vector size does not match!";
        return 1;
    }
    double sum = 0;
    for(size_t i = 0; i < a.size(); i++){
        sum += a[i]*b[i];
    }
    return sum;
}

inline void CSR_matrix(std::vector <double> &a, std::vector <size_t> &col, std::vector <size_t> &row, double delta, std::size_t n){
    std::size_t nnz_max = 3*n;
    size_t nnz = 0;

    for(std::size_t i = 0; i < n; i++){
        row[i] = nnz; //because of the size_t type
        if(i > 0){
            a[nnz] = 1/(delta*delta);
            col[nnz] = i-1;
            row[i] = nnz;
            nnz++;
        }
        a[nnz] = -2/(delta*delta);
        col[nnz] = i;
        if(row[i] < 0){
            row[i] = nnz;
        }
        nnz++;
        if(i < n-1){
            a[nnz] = 1/(delta*delta);
            col[nnz] = i+1;
            nnz++;
        }
    }
    row[n] = nnz;
}

inline void matvec(std::vector <double> x, std::vector <double> &y, std::vector <double> a, std::vector <size_t> col, std::vector <size_t> row, size_t n){
    for(std::size_t i = 0; i < n; i++){
        y[i] = 0;
        for(size_t l = row[i]; l < row[i+1]; l++){
            size_t j = col[l];
            y[i] = y[i]+a[l]*x[j];
        }
    }
}


#endif