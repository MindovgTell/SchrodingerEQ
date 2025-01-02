#ifndef CRANKNICOLSON_HPP   
#define CRANKNICOLSON_HPP   

#include <iostream>
#include <complex>
#include <fstream>
#include <cmath>
#include "../eigen/Eigen/Dense"
#include "../eigen/Eigen/Sparse"
//#include "../eigen/Eigen/IterativeLinearSolvers"
//#include "../eigen/Eigen/SparseCholesky"


class CrankNicolson
{
private:
    Eigen::SparseMatrix<std::complex<double> > m_A;
    Eigen::SparseMatrix<std::complex<double> > m_B;

    Eigen::VectorXcd m_Psi;
    Eigen::MatrixXd m_V;
    int m_size, m_T, t_step;
    double m_delta_t, m_h_step,V_0;
    std::complex<double> m_r;

public:
    //Constructor
    CrankNicolson(double h, double deltat, double T, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, double v_0, int slits=0);
    //Wave funciton
    std::complex<double> gauss_wave_packet(double sigma_x, double sigma_y, double x, double y, double x_c, double y_c, double p_x, double p_y = 0);

    //Methods for creating matrixes
    void init_start_state(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y);
    int get_m_index(int i,int j, int M);

    void init_time_evolution_matrices();

    void init_Mat_A(std::complex<double> r,Eigen::VectorXcd& d);
    void init_Mat_B(std::complex<double> r,Eigen::VectorXcd& d);

    Eigen::MatrixXd vec_to_mat(const Eigen::VectorXd& vec);

    Eigen::VectorXd prob(Eigen::VectorXcd&  vec);

    void save_to_csv(std::string filename, Eigen::MatrixXd mat);
    void simulation();
    
    // Functions for create potential
    Eigen::MatrixXd create_potential_box();
    Eigen::MatrixXd create_one_slit();
    Eigen::MatrixXd create_double_slit();
    Eigen::MatrixXd create_triple_slit();
};


#endif 