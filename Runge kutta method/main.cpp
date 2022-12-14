//
//  main.cpp
//  Runge kutta method
//
//  Created by WC on 2022/12/1.
//
#include "Runge kutta method.cpp"
#include <fstream>

typedef std::vector <double > state_type;
typedef runge_kutta3 < state_type > rk3_type;
typedef runge_kutta4 < state_type > rk4_type;

struct lorenz {
    const double sigma, R, b;
    lorenz(const double sigma, const double R, const double b): sigma(sigma), R(R), b(b)   { }
    void operator ()( const state_type &x, state_type &dxdt, double t) {
        dxdt [0] = sigma * ( x[1] - x[0] );
        dxdt [1] = R * x[0] - x[1] - x[0] * x[2];
        dxdt [2] = -b * x[2] + x[0] * x[1];
    }
};

int main () {
    std::ofstream lor_out;
    rk3_type stepper3;
    rk4_type stepper4;
    
    const int steps = 4500;
    const double dt = 0.01;
    lorenz system (10.0 , 28.0, 10.0/3.0);
    state_type x = {10.0, 1.0, 1.0};//initial condition
    
    //using 3nd order runge kutta
    lor_out.open("/Users/WC/Desktop/Runge kutta method/Runge kutta method/lorenz_3.txt");
    lor_out << std::fixed;
    for( size_t n=0 ; n<steps ; ++n ) {
    stepper3.do_step(system , x, n*dt , dt);
        lor_out<< x[0] << ' ' << x[1] << ' ' << x[2] << '\n';
    }
        lor_out.close();

    x = {10.0, 1.0, 1.0};
    
    //using 4nd order runge kutta
    lor_out.open("/Users/WC/Desktop/Runge kutta method/Runge kutta method/lorenz_4.txt");
    lor_out << std::fixed;
    for( size_t n=0 ; n<steps ; ++n ) {
    stepper4.do_step(system , x, n*dt , dt);
        lor_out<< x[0] << ' ' << x[1] << ' ' << x[2] << '\n';
    }
        lor_out.close();
}
