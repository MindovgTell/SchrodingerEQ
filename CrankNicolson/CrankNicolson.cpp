#include "CrankNicolson.hpp"

using namespace std::complex_literals;

//Constructor
CrankNicolson::CrankNicolson(double h, double deltat, double T, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, double v_0, int slits){
    int M = 1/h;
    this->m_h_step = h;
    this->m_delta_t = deltat;
    this->m_T = T;
    this->V_0 = v_0;
    this->m_r = -1i*m_delta_t/(2*std::pow(m_h_step,2));
    this->t_step = std::round(T/deltat) + 1;
    this->m_size = M;

    if(slits == 0){
        this->m_V = this->create_potential_box();
    }
    else if(slits == 1){
        this->m_V = this->create_one_slit();
    }
    else if(slits == 2){
        this->m_V = this->create_double_slit();
    }
    save_to_csv("./Matrice/potential.csv", this->m_V);

    this->init_time_evolution_matrices();

    this->init_start_state(x_c,y_c,sigma_x,sigma_y,p_x,p_y);
}

//Indexes in vector representation of wave function
int CrankNicolson::get_m_index(int i, int j, int M){
    return (i-1)*(M-2) + j-1;
}

//The value of the wave funciton in the specific point on grid
std::complex<double> CrankNicolson::gauss_wave_packet(double sigma_x, double sigma_y, double x, double y, double x_c, double y_c, double p_x, double p_y){
    std::complex<double> i(0, 1); // Define the imaginary unit
    double exponent = -(pow(x - x_c ,2) / (2 * pow(sigma_x,2))) - (pow(y - y_c, 2) / (2*pow(sigma_y, 2)));
    std::complex<double> phase = i * (p_x * (x - x_c) + p_y * (y - y_c));
    return std::exp(exponent + phase); 
}

//Function for initializing wave function
void CrankNicolson::init_start_state(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y)
{
    int size = std::pow(this->m_size-2,2);
    Eigen::VectorXcd U(size);
    std::complex<double> psum = 0;
    for(int i = 1; i != m_size-1; ++i){
        double x = i*m_h_step;
        for(int j = 1; j != m_size-1; ++j){
            double y = j*m_h_step;
            std::complex<double> c = gauss_wave_packet(sigma_x,sigma_y,x,y,x_c,y_c,p_x,p_y);
            int index = get_m_index(i,j,this->m_size);
            U(index) = c;
            psum += std::real(std::conj(c)*c);
        }
    }

    std::complex<double> normalization_factor = 1.0 / std::sqrt(psum);
    U = normalization_factor * U;
    this->m_Psi = U;
}

//Function for constructing the right-hand and left-hand sides matrices from Crank Nicolson algorithm
void CrankNicolson::init_time_evolution_matrices(){
    int mat_size = pow(this->m_size-2,2);
    Eigen::VectorXcd a(mat_size);
    Eigen::VectorXcd b(mat_size);

    for(int k = 1; k < this->m_size-1; ++k){
        for(int l = 1; l < this->m_size-1; ++l){
            int index = get_m_index(k,l,m_size);
            a(index) = (1.0 - 4.0*this->m_r + 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));
            b(index) = (1.0 + 4.0*this->m_r - 1.0i*(m_delta_t/2)*std::complex<double>(m_V(l,k)));
        }
    }
    this->init_Mat_A(m_r,a);
    this->init_Mat_B(-m_r,b); 
}

//Function which return Eigen::MatrixXd based on Eigen::VectorXd
Eigen::MatrixXd CrankNicolson::vec_to_mat(const Eigen::VectorXd& vec) {   
    int size = std::sqrt(vec.size());
    Eigen::MatrixXd mat(size,size);
    Eigen::Map<const Eigen::MatrixXd> mat_map(vec.data(), size, size);
    mat = mat_map;
    
    return mat;
}

//Function for calculate the density of probability in each point of space
Eigen::VectorXd CrankNicolson::prob(Eigen::VectorXcd &vec)
{
    int size = std::pow(this->m_size-2,2);
    Eigen::VectorXd pr(size);
    for(int i = 0; i != size; ++i){
        pr(i) = std::pow(vec(i).real(),2) + std::pow(vec(i).imag(),2);
    }
    return pr;
}

//Function for initialization left-hand side matrix according to Crank Nicolson algorithm
void CrankNicolson::init_Mat_A(std::complex<double> r, Eigen::VectorXcd& d){
    int S = d.size();
    int s = std::sqrt(S);

    typedef Eigen::Triplet<std::complex<double> > T;
    std::vector<T> tripletList;
    tripletList.reserve(5*s);

    tripletList.push_back(T(0, 0, d(0)));
    tripletList.push_back(T(0, 1, r));
    tripletList.push_back(T(S - 1, S - 2, r));
    tripletList.push_back(T(S - 1, S - 1, d(S-1)));

    for (int i = 1; i < S-1; ++i) {
        if(i + s  < S){
            tripletList.push_back(T(i-1, i+s-1, r));
            tripletList.push_back(T(i+s-1, i-1, r));
        } 
        tripletList.push_back(T(i, i,d(i)));

        if(i%s == 0){
            std::complex<double> z(0.,0.);
            tripletList.push_back(T(i, i - 1, z));
            tripletList.push_back(T(i-1, i, z));
        }
        else {
            tripletList.push_back(T(i, i - 1, r));
            tripletList.push_back(T(i - 1, i, r));
        }
    }

    Eigen::SparseMatrix<std::complex<double> > A(S,S);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    this->m_A = A;
}

//Function for initialization right-hand side matrix according to Crank Nicolson algorithm
void CrankNicolson::init_Mat_B(std::complex<double> r, Eigen::VectorXcd& d) {
    int S = d.size();
    int s = std::sqrt(S);

    typedef Eigen::Triplet<std::complex<double> > T;
    std::vector<T> tripletList;
    tripletList.reserve(5*s);

    tripletList.push_back(T(0, 0, d(0)));
    tripletList.push_back(T(0, 1, r));
    tripletList.push_back(T(S - 1, S - 2, r));
    tripletList.push_back(T(S - 1, S - 1, d(S-1)));

    for (int i = 1; i < S-1; ++i) {
        if(i + s  < S){
            tripletList.push_back(T(i-1, i+s-1, r));
            tripletList.push_back(T(i+s-1, i-1, r));
        } 
        tripletList.push_back(T(i, i,d(i)));

        if(i%s == 0){
            std::complex<double> z(0.,0.);
            tripletList.push_back(T(i, i - 1, z));
            tripletList.push_back(T(i-1, i, z));
        }
        else {
            tripletList.push_back(T(i, i - 1, r));
            tripletList.push_back(T(i - 1, i, r));
        }
    }

    Eigen::SparseMatrix<std::complex<double> > B(S,S);
    B.setFromTriplets(tripletList.begin(), tripletList.end());
    this->m_B = B;
}

//Function for saving Eigen Matrices as csv tabels 
void CrankNicolson::save_to_csv(std::string filename, Eigen::MatrixXd mat){

    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ",","\n");
    int size = std::sqrt(mat.size());
    Eigen::MatrixXd vec(1,size);

    std::ofstream file(filename);
    if(file.is_open()){
        //file << vec.format(CSVFormat) << '\n';
        file << mat.format(CSVFormat);
        file.close();
    }
}

// Function for solving systems of equations for each time step dependently on start conditions
void CrankNicolson::simulation()
{
    int size = this->m_Psi.size();

    Eigen::VectorXcd x(size);
    Eigen::VectorXcd b(size);

    x.setZero();
    b.setZero();

    // Save initial data before the loop
    Eigen::VectorXd p_init = prob(m_Psi);
    save_to_csv("./Matrice/matrix0.csv", vec_to_mat(p_init));

    // Prepare the right-hand side for the time-stepping
    b = (this->m_Psi);

    // Set up the sparse LU solver
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> lg;
    lg.compute(m_A);
    

    for (int i = 1; i < this->t_step; ++i) {
        // Update the right-hand side vector b
        b = (this->m_B) * b;

        // Solve the system A * x = b
        x = lg.solve(b);

        // Update b for the next iteration
        b = x;

        // Calculate probability and save to CSV
        Eigen::VectorXd p = prob(x);
        Eigen::MatrixXd solution = vec_to_mat(p);
        save_to_csv("./Matrice/matrix" + std::to_string(i) + ".csv", solution);   
        
        if(i % 4 == 0)
            std::cout << "Simulation step: #" << i << '\n';
    }
}

//Function for creating potential box without the wall
Eigen::MatrixXd CrankNicolson::create_potential_box(){
    Eigen::MatrixXd V(this->m_size,this->m_size);
    V.setZero();
    Eigen::VectorXd S(this->m_size);
    S.fill(V_0);
    V.col(0) = S;
    V.col(this->m_size-1) = S;
    V.row(0) = S;
    V.row(this->m_size-1) = S;

    return V;
}

//Function for creating potential box with one slit wall in the middle
Eigen::MatrixXd CrankNicolson::create_one_slit(){
    Eigen::MatrixXd V = create_potential_box();
    int center_index = (m_size)*0.5; //200 * 0.5 = 100
    int x_thickness = 0.02/m_h_step;// indices i in x direction: (0.02/0.005) + 1 = 5
    int x_start = center_index - x_thickness/2;
    int x_end = center_index + x_thickness/2;
    int aperture = (0.05/m_h_step) ;// (0.05/0.005) + 1 = 11
    int start = center_index - aperture/2;
    int end = start + aperture + 1;
    Eigen::VectorXd S(this->m_size);
    S.fill(this->V_0);

    for (int i = x_start; i < x_end+1; ++i)
    {
        V.col(i) = S;
        for(int j = start; j < end; ++j)
            V(j,i) = 0;
    }
    return V;
}

//Function for creating potential box with double slit wall in the middle
Eigen::MatrixXd CrankNicolson::create_double_slit(){

    Eigen::MatrixXd V = create_potential_box();

    int center_index = (this->m_size)*0.5; //200 * 0.5 = 100
    int x_thickness = 0.02/m_h_step;// indices i in x direction: (0.02/0.005) + 1 = 5
    int x_start = center_index - x_thickness/2;
    int x_end = center_index + x_thickness/2;
    int aperture = (0.05/m_h_step) + 1 ;// (0.05/0.005) + 1 = 11
    int center_wall_length = (0.05/m_h_step);// (0.05/0.005) = 10. From j=95 to j=105
    int start = center_index - (center_wall_length)/2 - aperture; // j=84 for lower aperture and j=106 for upper aperture
    int end = start + aperture;

    Eigen::VectorXd S(this->m_size);
    S.fill(this->V_0);

    for (int i = x_start; i < x_end+1; i++) {
        V.col(i) = S;
        for(int j = start; j < end+1; j++) {
            V(j,i) = 0;
            V(j+center_wall_length+1+aperture,i) = 0;
        }
    }
    return V;
}


Eigen::MatrixXd CrankNicolson::create_triple_slit(){
    Eigen::MatrixXd V = create_potential_box();
    
    int center_index = (this->m_size)*0.5; //200 * 0.5 = 100
    int x_thickness = 0.02/m_h_step;// indices i in x direction: (0.02/0.005) + 1 = 5
    int x_start = center_index - x_thickness/2;
    int x_end = center_index + x_thickness/2;
    int aperture = (0.05/m_h_step) + 1 ;// (0.05/0.005) + 1 = 11
    int center_wall_length = (0.05/m_h_step);// (0.05/0.005) = 10. From j=95 to j=105
    int start = center_index - (center_wall_length)/2 - aperture; // j=84 for lower aperture and j=106 for upper aperture
    int end = start + aperture;

    Eigen::VectorXd S(this->m_size);
    S.fill(this->V_0);

    for (int i = x_start; i < x_end+1; i++) {
        V.col(i) = S;
        for(int j = start; j < end+1; j++) {
            V(j,i) = 0;
            V(j+center_wall_length+1+aperture,i) = 0;

        }
    }
    
    return V;    
}

