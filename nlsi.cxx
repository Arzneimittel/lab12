#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
// #include <fftw3.h>
// #include <cmath>

//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double eta, const double sigma, const double dx,
          const int Nx);

void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin);

void lin_step(complex <double>* const f1, complex <double>* const f0,
          const double dt, const double dx, const int Nx);
void nonlin_step(complex <double>* const f1, complex <double>* const f0,
          const double dt, const int Nx);
//-----------------------------------
int main(){

	const int Nx = 4000;
	const double L = 800;
	const double xmin = 0;
	const double Tend = 50;
	const double dx = L / (Nx - 1);
	const double dt = dx  / 10;
	const int Na = 10;
	double t=0.0;
	int Nk = int(Tend / Na / dt + 0.5);

	const double eta = 0.2;

	stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h = new cmplx[Nx];
	

	init(psi0, eta, dx, dt,Nx);

	writeToFile(psi0,"psi_0", dx,Nx,xmin);


	for (int i = 1; i <= Na; i++) {

		for (int j = 1; j <= Nk-2; j++) {
		  
		  lin_step(psi1,psi0,0.5*dt,dx,Nx);
                  nonlin_step(psi0,psi1,dt,Nx);
		  lin_step(psi1,psi0,0.5*dt,dx,Nx);
	      h = psi0;
 	      psi0 = psi1;
 	      psi1 = h;
	      t +=dt;
		  
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin);
	}
	delete[] psi0;
	delete[] psi1;
        
	return 0;
}
//-----------------------------------
void nonlin_step(complex <double>* const f1, complex <double>* const f0,
          const double dt, const int Nx){
	 for(int i=0;i<Nx;i++){
	   
	   f1[i] = f0[i] * exp(cmplx(0.0, - pow(norm(f0[i]),1) * dt));
	   
	}
	  }


void lin_step(complex <double>* const f1, complex <double>* const f0,
          const double dt, const double dx, const int Nx)
{

  complex<double>* d=new complex <double>[Nx];
  complex<double>* l=new complex <double>[Nx];
  complex<double>* u=new complex <double>[Nx];
  

  for(int i=0;i<Nx;i++) d[i] = cmplx(1.0, - 2.0*dt/(dx*dx));
  for(int i=0;i<Nx;i++) u[i] = cmplx(0.0,+ dt/(dx*dx));
  for(int i=0;i<Nx;i++) l[i] = cmplx(0.0,+ dt/(dx*dx));
  
  // Forward substitution
   complex <double> tri, zw1, zw2;
  
  for(int i=0; i<Nx-1; i++){
    tri = l[0] / d[i] ;
    zw1 = tri * u[i];
    d[i+1] -= zw1;
    
    zw2 = tri * f0[i];
    f0[i+1] -=  zw2;
   
    
  }
  // Backward substitution
  f1[Nx-1] = f0[Nx-1]/d[Nx-1];
  
  for(int i=Nx-2; i>=0; i--){
    f1[i] = (f0[i] - u[0]*f1[i+1]) / d[i];   
  }
  

  
  delete[] d;
  delete[] u;
  delete[] l;
}
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin)
{
	ofstream out(s.c_str());
	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx;
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag() << endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double eta,  const double dx, const double dt,
          const int Nx)
{
	const double x0 = dx*Nx * 0.5;
	const double f = sqrt(2) * eta;
	for(int i=0;i<Nx; i++){
		double x = i*dx - x0;
		psi0[i] = 2*f/cosh(eta * x);
	}
}
