#include <Rcpp.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

double phi2(NumericVector Om, double tau_sq, double sigma_sq, double rho){
    double eta = 0;
    for(int i = 1; i < Om.length(); i++){
        double nSamp = i+1;
        double sigma_n_sq = 2*sigma_sq/(pow(nSamp,2)-nSamp) + (4*nSamp-8)*rho/(pow(nSamp,2)-nSamp);
        eta = eta + (Om[i]/(tau_sq + sigma_n_sq));
    }
    return eta;
}

// [[Rcpp::export]]
List recursive_search(int A, int M, int R, NumericVector Om, const NumericVector& N, 
             const double& tau_sq, const double& sigma_sq, const double& rho){

    NumericVector search_grid (N[A-1]+R+1);
    for(int i=0; i < search_grid.length(); ++i){
        search_grid[i] = i;
    }
    
    if (A == 2){
        double a = std::floor(M/(double)2);
        double b = max(search_grid);
        Om[A-1] = std::min(a,b);
        return List::create(Om,phi2(Om,tau_sq,sigma_sq,rho));
    }
    
    else {
        if (M <= 0){
            return List::create(Om,phi2(Om,tau_sq,sigma_sq,rho));
        }
        else {
            double eta_max = 0;
            NumericVector Om_max;
            search_grid = search_grid[A*search_grid <= M];
            for(int i = 0; i < search_grid.length(); ++i){
                int Om_A = search_grid[i];
                Om[A-1] = Om_A;
                int new_M = M-(A*Om_A);
                int new_R = N[A-1] + R - Om_A;
                NumericVector Om_c = clone(Om);
                List temp_res = recursive_search(A-1,new_M,new_R,Om_c,N,tau_sq,sigma_sq,rho);
                double eta = temp_res[1];
                if(eta > eta_max){
                    eta_max = eta;
                    Om_max = temp_res[0];
                }
            }
            return List::create(Om_max,eta_max);
        }
    }
}