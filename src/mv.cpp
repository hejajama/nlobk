#include "mv.hpp"
#include <sstream>
#include <string>
#include <cmath>
#include "nlobk_config.hpp"
#include <iostream>
#include "solver.hpp"   
#include <tools/tools.hpp>

using std::cerr; using std::endl;

double MV::DipoleAmplitude(double r, double b)
{
/*
 * cerr << "Note: MVgamma model is modified!!!" << endl;
	
	// rcmv from 1507.03651
	BKSolver solver;	// Becaus Alpha_s is defined there
*/
	// GBW
    /*
    double q0sqr=std::pow(2.0*0.241,2); // Qs_0^2 = 2lambdaQCD
	double n = 1.0 - std::exp(-r*r*q0sqr);
	if (n<0) return 0;
	if (n>1) return 1;
    return n;
     */
	

	if (b>1e-10)
		cerr << "Impact parameter is not supported!" << LINEINFO << endl;
	const double e = 2.7182818;
	///TODO: some algorithm to determina small r, e.g. when one has to linearize
    if (r < 2e-6)   ///NOTE: factor 1/4 "correctly", not as in AAMS paper
            return std::pow(SQR(r)*qs0sqr, anomalous_dimension)/4.0
            * std::log( 1.0/(r*lambdaqcd) + ec*e) ;
    return 1.0 - std::exp(-std::pow(SQR(r)*qs0sqr, anomalous_dimension)/4.0
            * std::log( 1.0/(r*lambdaqcd) + ec*e) );
}


void MV::SetQsqr(double qsqr)
{
	qs0sqr=qsqr;
}

void MV::SetAnomalousDimension(double gamma_)
{
	anomalous_dimension=gamma_;
}

void MV::SetLambdaQcd(double lambda)
{
	lambdaqcd=lambda;
}

void MV::SetE(double ec_)
{
	ec=ec_;
}

double MV::GetE()
{
    return ec;
}

std::string MV::GetString()
{
	std::stringstream ss;
	ss << "MV model, Q_s0^2 = " << qs0sqr << " GeV^2, \\gamma = " << anomalous_dimension
		<< ", coefficient of E inside Log is " << ec 
		<< ", x0=" << x0 << ", \\Lambda_QCD = " << lambdaqcd << " GeV";
	return ss.str();
}

/*
 * Set some reasonable parameters
 * That is, MV1 model for nucleus
 */
MV::MV()
{
	qs0sqr = 0.10;
	x0=0.01;
	ec=1.0;
	lambdaqcd=config::LAMBDAQCD;
	anomalous_dimension=1;
}
