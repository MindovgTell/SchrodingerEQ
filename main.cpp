#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <complex>
#include "./CrankNicolson/CrankNicolson.hpp"

using namespace std::complex_literals;

void simulation(std::string inputfile);

int main(int argc, char const *argv[]){

    if(argc != 2){
        std::string executable = argv[0];

        std::cerr << "Error: Wrong number of input parameters" << '\n';
        std::cerr << "Usage:" << executable << " input_file.txt" << '\n';
        return 1;
    }
    simulation(argv[1]);
  
    return 0;
}


void simulation(std::string inputfile){
    std::ifstream input_data(inputfile);

    std::string line;
    std::getline(input_data,line);//Skip first line in file

    double Problem,h,deltat,T,x_c,sigma_x,p_x,y_c,sigma_y,p_y,v_0, slits;

    std::getline(input_data,line);
    std::stringstream str_stream(line);
    str_stream >> Problem >> h >>deltat >>T >> x_c >> sigma_x >> p_x >> y_c >> sigma_y >> p_y >> v_0 >> slits;

    int width = 10;
    std::cout << std::setw(width) << "Problem" << std::setw(width) << "h" << std::setw(width) << "deltat" << std::setw(width) << "T" << std::setw(width) << "x_c" << std::setw(width)
    << "sigma_x" << std::setw(width) << "p_x" << std::setw(width) << "y_c" << std::setw(width) << "sigma_y" << std::setw(width) << "p_y" << std::setw(width) << "v_0"
    << std::setw(width) << "slits" << std::endl;

    std::cout << std::setw(width) << Problem << std::setw(width) << h << std::setw(width) << deltat << std::setw(width) << T << std::setw(width) << x_c << std::setw(width)
    << sigma_x << std::setw(width) << p_x << std::setw(width) << y_c << std::setw(width) << sigma_y << std::setw(width) << p_y << std::setw(width) << v_0 << std::setw(width) << slits
    << std::endl;

    CrankNicolson Crank(h, deltat, T, x_c, y_c, sigma_x, sigma_y,p_x, p_y,v_0, slits);

    Crank.simulation();

}

