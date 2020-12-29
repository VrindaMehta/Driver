#include "basic.h"
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <vector>
using namespace std;
typedef complex<double> dcomp;

#define lapack_complex_float complex<float>
#define lapack_complex_double complex<double>

#include <mkl_lapacke.h>


int main(int argc, char *argv[])    //argv[]=[readfile,t,dt]
{

    ::iota = -1;
    ::iota = sqrt(::iota);
    readfunc(argv[1]);
    int num_ev = 4;
    int i, j;
    double t, dt, step;
    //vector<complex<double> > initialstate(p);

    t = stod(argv[2]);
    dt = stod(argv[3]);
    cout<<endl;
    step = t/dt;
    cout<<"#t =  "<<t<<endl;
    cout<<"#dt = "<<dt<<endl;
    cout<<"#Number of time steps = "<<step<<endl;
    cout<<endl;


    //results_spectrum(dt, num_ev, step);
    results_productevolution(dt, step);

    delete [] h_x;
    delete [] h_z;
    delete [] J_x;
    delete [] J_z;
    delete [] K_x;
    delete [] flag;
    return 0;

}
