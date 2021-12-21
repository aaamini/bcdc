#include <iostream>
#include <vector>
#include <armadillo>

// #include "utils.h"
// #include "sampling.h"
// #include "beta_calcs.h"
// #include "dpsbm.h"


arma::uvec sample_int_vec(const int K, const int size) {
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(0, K-1);

    std::vector<int> out(size);
    std::generate(std::begin(out), std::end(out), [&](){
        return distrib(gen);
    });

    return arma::conv_to<arma::uvec>::from(out);    
}

int main() {

    std::vector<arma::sp_mat> A(3);
    std::generate(std::begin(A), std::end(A), [](){ return arma::sprandu(4,4,0.5); } );

    arma::mat B(5, 5, arma::fill::randn);
    // for (auto &e: A) {
    //     std::cout << e;
    // }
    auto idx1 = arma::span(0,2);
    auto idx2 = arma::span(2,4);
    auto idx3 = arma::regspace<arma::uvec>(0,2);

    // auto idx4 = arma::join_vert(arma::regspace(1,2), 
    //     arma::regspace(3,4));
    

    std::cout << B << "\n" // << B(idx3, idx4) 
        << "\n" << B(idx1,3) << "\n";

    std::cout << B.row(1) * B.col(2) << "\n";

    auto temp = B.row(1) * B.col(2);    

    std::cout << temp << "\n" << typeid(temp).name();

    // std::cout << models[1].z.t();

}
