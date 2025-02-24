#ifndef PORT_CONTROLLER_H
#define PORT_CONTROLLER_H

#include <Eigen/Dense>
#include "readFile.h"
#include "LieAlgebra.h"
#include "Spline.h"
#include "ConfigFile.h"
#include "Chameleon.h"
#include "SmoothStep.h"
#include "readFile.h"
#include "DEBUGGER.h"

using namespace std;
using namespace Eigen;

class PortController 
{
	public:

	int Na, Nm;

	typedef Matrix<double, Dynamic, Dynamic> Mnd;
	typedef Matrix<double, 6, Dynamic> Mad;
	typedef Matrix<double, Dynamic, 1> Vnd;
	typedef Matrix<double, 6, 1> V6d;
	typedef Matrix<double, 2, 1> V2d;

	double t;	// time instance
	double h;   // time-stepping
	double T;	// simulation time
	Vnd q; 		// generalized coordinates
	Vnd p; 		// generalized momenta
	Vnd u; 		// generalized control output
	Vnd uk;		// kalman observer output

	V2d H;	    // Hamiltonian
	Vnd dHdq;   // Hamiltonian-Grad q
	Vnd dHdp;   // Hamiltonian-Grad p

	Mnd Dee;	    // dampings matrix
	Mnd Sa;		    // actuation matrix
	double Kp, Kd;  // controller gains
	double K1, K2;  // virtual stiffness
	double Lambda;  // artif. damping gain
	double L;       // observer gain
	double LambdaK; // artif. observer damping gain
	V7d gd;		    // desired pos 

	S1d spline_X;

	PortController();

	void setup(
		const char* str);

	void EnergyShaping(
		Mad J,
		V7d g);

	void KalmanCorrection(
		Mad J,
		V7d g,
		V7d z,
		V6d ed);

};

#endif