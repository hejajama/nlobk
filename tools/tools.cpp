/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2015
 */

#include "tools.hpp"
#include "config.hpp"
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <gsl/gsl_integration.h>

using namespace Amplitude;

/*
 * Str to double/int
 */
double StrToReal(std::string str)
{
    std::stringstream buff(str);
    double tmp;
    buff >> tmp;
    return tmp;
}

int StrToInt(std::string str)
{
    std::stringstream buff(str);
    int tmp;
    buff >> tmp;
    return tmp;
}

// GSL Error handler
int errors;
void ErrHandler(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno)
{
    
    // Errors related to convergence of integrals are handled when
    // gsl_integration functions are called, don't do anything with them here
    // 14 = failed to reach tolerance
    // 18 = roundoff error prevents tolerance from being achieved
    // 11 = maximum number of subdivisions reached
    if (gsl_errno == 14 or gsl_errno == 18 or gsl_errno == 11)
        return;

    // 15: underflows
    if (gsl_errno == 15 ) return;
    // 16: overflows
    if (gsl_errno == 16 ) return;


    errors++;
    std::cerr << file << ":"<< line <<": Error " << errors << ": " <<reason
            << " (code " << gsl_errno << ")." << std::endl;
}


/* Returns index i for which
 * vec[i]<=val
 * Assumes that vec[i]<vec[i+1]
 * If such index can't be found, returns -1
 */

int FindIndex(double val, std::vector<double> &vec)
{
    if (val < vec[0]) return -1;
    
    int ind=-1;
    
    uint start=0; uint end=vec.size()-1;
    while(end-start>5)
    {
        int tmp = static_cast<int>((start+end)/2.0);
        
        if (vec[tmp]>=val)
            end=tmp;
        else
            start=tmp;
    }
    
    
    for (uint i=start; i<=end; i++)
    {
        if (vec[i]<=val and vec[i+1]>val)
        {
            ind=i;
            break;
        }
    }
    if (ind == -1) return vec.size()-1;
    return ind;
}

/*
 * Nuclear density profile
 * Normalization: \int d^2 b T_A(b)=1
 * All lenght units: 1/GeV!
 */

// Woodls-Saxon distribution
double w_s_normalization=1;
double w_s_normalization_A=-1;
double W_S(double r, int A)
{
	if (A != w_s_normalization_A)
	{
		cerr << "Woods-Saxon distribution is not initialized for A="<< A << " " << LINEINFO << endl;
		return 0;
	}
	// R = 1.12 fm A^(1/3) - 0.86 fm * A^(-1/3)
	double ra = 1.12 * std::pow(A, 1.0/3.0) - 0.86 * std::pow(A, -1.0/3.0);
	double delta = 0.54 * FMGEV;
	ra *= FMGEV;	// fm => 1/GeV

	return w_s_normalization / (std::exp((r - ra)/delta)+1);
}

struct inthelper_ta
{
	double b;
	int A;
};

// Return W_S(\sqrt{ b^2 + z^2 } )
double inthelperf_ta(double z, void* p)
{
	inthelper_ta* par = (inthelper_ta*)p;
	return W_S(std::sqrt( SQR(par->b) + SQR(z) ), par->A); 
}

double T_A(double b, int A)
{
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(10);
	double res, abserr;
	inthelper_ta par; par.A=A; par.b=b;
	gsl_function f; f.params=&par; f.function=inthelperf_ta;
	int status = gsl_integration_qag(&f, 0, 99, 0, 0.001, 10, GSL_INTEG_GAUSS61,
		w, &res, &abserr);
    gsl_integration_workspace_free(w);
	if(status)
		cerr << "T_A integration failed at " << LINEINFO <<", result " << res
			<< ", relerr " << std::abs(res-abserr)/res <<", A=" << A << endl;
	return 2.0*res;	// 2.0 as we integrate z in [0,\infty]
}

// WS initialization, calculates normalization s.t. 
// \int d^3 b WS(b)=1
double inthelperf_ws(double r, void* p)
{
	return r*r*W_S(r, *((int*)p));
}
	
void InitializeWSDistribution(int A)
{
	w_s_normalization = 1;
	w_s_normalization_A = A;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(10);
	double res, abserr;
	gsl_function f; f.params=&A; f.function=inthelperf_ws;
	int status = gsl_integration_qag(&f, 0, 99, 0, 0.001, 10, GSL_INTEG_GAUSS61,
		w, &res, &abserr);
	if (status)
		cerr << "WS normalization integral failed at " << LINEINFO <<", result "
			<< res << " relerr " << std::abs(res-abserr)/res << endl;
	
	res *= 4.0*M_PI;
	
	w_s_normalization = 1.0/res;
    gsl_integration_workspace_free(w);
}

/*
 * Subtract the smallest element of the vector from the array
 * Forcest that vec[min]=0
 */
void SubtractMinimum(std::vector<double> &array)
{
	if (array.size()==0) return;
	
	double min = array[0];
	for (uint i=1; i<array.size(); i++)
	{
		if (array[i]<min)
			min=array[i];
	}
	
	for (uint i=0; i<array.size(); i++)
		array[i]=array[i]-min;
}

std::string PartonToString(Amplitude::Parton p)
{
    switch(p)
    {
        case U:
            return "u";
		case UBAR:
			return "ubar";
        case D:
            return "d";
		case DBAR:
			return "dbar";
        case S:
            return "s";
		case SBAR:
			return "sbar";
        case LIGHT:
            return "Light (u+d+s)";
        case C:
            return "c";
		case CBAR:
			return "cbar";
        case B:
            return "b";
		case BBAR:
			return "bpar";
        case G:
            return "g";        
        default:
            return "Unknown";
    }
}
