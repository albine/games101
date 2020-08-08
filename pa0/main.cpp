#include<cmath>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<iostream>

#define PI 3.14159265

int main(){

    // Basic Example of cpp
    std::cout << "Example of cpp \n";
    float a = 1.0, b = 2.0;
    std::cout << a << std::endl;
    std::cout << a/b << std::endl;
    std::cout << std::sqrt(b) << std::endl;
    std::cout << std::acos(-1) << std::endl;
    std::cout << std::sin(30.0/180.0*acos(-1)) << std::endl;


    // vector definition
    // P = (2, 1) in homogenous coordinate 
    Eigen::Vector3f p(2.0f, 1.0f, 1.0f);

    // matrix definition
    Eigen::Matrix3f m;
    m << std::cos(PI / 4), std::sin(- PI / 4), 1.0, std::sin(PI / 4), std::cos(PI / 4), 2.0, 0.0, 0.0, 1.0 ;

    // matrix multiply vector i * v
    std::cout << "Example of output \n";
    std::cout << m * p << std::endl;

    return 0;
}