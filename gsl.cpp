// [[Rcpp::depends(RcppGSL)]]
#include <iostream>
#include <RcppGSL.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <unistd.h>            // getpid

Rcpp::IntegerVector rmn( Rcpp::NumericVector p, gsl_rng* r){//rmn(unsigned int N, Rcpp::NumericVector p, gsl_rng* r)

    size_t K = p.size();

    Rcpp::IntegerVector x(K);
    //gsl_ran_multinomial(r, K, N, p.begin(), (unsigned int *) x.begin());
    gsl_ran_multinomial(r, K, 1, p.begin(), (unsigned int *) x.begin());
    //std::cout << x;
    //std::cout << " ";
    return x;             // return results vector
}
/*
Rcpp::IntegerVector gsl_mmm_1(Rcpp::NumericMatrix P, gsl_rng* r){
    //size_t K = N.size();

    int i;
    Rcpp::IntegerVector x(K);
    for(i=0; i<K; i++){
        //x += rmn(P(Rcpp::_, i), r);
        x = rmn(P(Rcpp::_, i), r);
        std::cout << x;
    }
    //std::cout << x;
    return x;
}
*/

// [[Rcpp::export]]
Rcpp::IntegerMatrix gsl_mmm(Rcpp::NumericMatrix P){
    int j;
    
    Rcpp::IntegerMatrix X(P.nrow(), P.ncol());
    for(j=0; j<P.nrow(); j++){
        gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
        long seed = rand()/(((double)RAND_MAX + 1)/10000000) * getpid();
        gsl_rng_set (r, seed);
        //X(Rcpp::_, j) = gsl_mmm_1(X_(Rcpp::_,j), P, r);
        X(j,Rcpp::_) = rmn(P(j,Rcpp::_), r);
        gsl_rng_free (r);
    }
    
    return X;
}