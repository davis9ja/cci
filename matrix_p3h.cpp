#include "dmcpp.hpp"
#include <Eigen/Eigenvalues>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

matrix<double> matrix_p3h(int nholes, int nparticles,
                          double zBc, double oBc, double tBc, double phBc) {
    int num_sp, num_states;
    matrix<int> basis;
    matrix<double> H;
    vector<int> bra, ket_delta_up, ket_delta_dn, ket_g, ket_pb, ket_pb_conj;
    
    num_sp = nholes+nparticles;
    basis = gen_basis(nholes, nparticles);
    num_states = basis.size1();
    H = matrix<double>(num_states,num_states);
    for (int i = 0; i < num_states; i++)
        for(int j = 0; j < num_states; j++)
            H(i,j) = 0.0;

    for (int i = 0; i < num_states; i++) {
        bra = matrix_row<matrix<int>>(basis, i);
        std::cout << bra << std::endl;
        for (int j = 0; j < num_states; j++) {
            ket_delta_up = matrix_row<matrix<int>>(basis, j);
            ket_delta_dn = matrix_row<matrix<int>>(basis, j);
            ket_g = matrix_row<matrix<int>>(basis, j);
            ket_pb = matrix_row<matrix<int>>(basis, j);
            ket_pb_conj = matrix_row<matrix<int>>(basis, j);
            
            for (int p = 0; p < num_sp; p+=2) {
                annihilator(ket_delta_up, p);
                creator(ket_delta_up, p);

                annihilator(ket_delta_dn, p+1);
                creator(ket_delta_dn, p+1);

                H(i,j) += oBc*(int)(p/2)*(inner_product(bra,ket_delta_up) + inner_product(bra,ket_delta_dn));

                for (int q = 0; q < num_sp; q+=2) {
                    annihilator(ket_g, q);
                    annihilator(ket_g, q+1);
                    creator(ket_g, p+1);
                    creator(ket_g, p);

                    H(i,j) += -tBc/2.*inner_product(bra,ket_g);

                    for (int pp = 0; pp < num_sp; pp+=2) {

                        // enforce conditional sum in pair-breaking term
                        if (pp != p) {
                            annihilator(ket_pb, q);
                            annihilator(ket_pb, q+1);
                            creator(ket_pb, pp+1);
                            creator(ket_pb, p);

                            annihilator(ket_pb_conj, pp);
                            annihilator(ket_pb_conj, p+1);
                            creator(ket_pb_conj, q+1);
                            creator(ket_pb_conj, q);
                            
                            H(i,j) += -phBc/2.*(inner_product(bra,ket_pb) + inner_product(bra,ket_pb_conj));
                        }
                    }
                }
            }
        }
    }
    
    return H;
    
}

int main() {
    int nholes, nparticles;
    double d, g, pb;
    matrix<double> H;

    nholes = 4;
    nparticles = 4;
    d = 1.0;
    g = 0.5;
    pb = 0.1;


    H = matrix_p3h(nholes,nparticles, 0.0, d, g, pb);

    const int dim = 36;
    Eigen::Matrix<double,dim,dim> mat;
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            mat(i,j) = H(i,j);

    std::cout << H << std::endl;

    // for (int i = 0; i < dim; i++) {
    //     for (int j = 0; j < dim; j++) {
    //         //std::cout << H(i,j) << " ";
    //         printf('%-0.4f', (double)H(i,j));
    //     }
    //     std::cout << "\n";
    // }
    // std::cout << mat << std::endl;
    // Eigen::EigenSolver<Eigen::Matrix<double,dim,dim>> eig(mat);
    // //eig.compute(mat);

    // std::cout << eig.eigenvalues() << std::endl;
}
