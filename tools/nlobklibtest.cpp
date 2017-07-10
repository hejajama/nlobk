
#include <iostream>
#include <amplitudelib/amplitudelib.hpp>
#include "solver.hpp"
#include "ic.hpp"
#include "mv.hpp"
#include "ic_datafile.hpp"
#include "dipole.hpp"
#include "solver.hpp"
#include <csignal>
#include <iostream>
#include <iomanip>
#include <gsl/gsl_errno.h>
#include <ctime>
#include <unistd.h>
#include <tools/tools.hpp>
#include <sstream>

void ErrHandler(const char * reason,
                const char * file,
                int line,
                int gsl_errno);


int main()
{
    gsl_set_error_handler(&ErrHandler);
    
    // ***Fit parameters***
    // Initial condition
    double qs0sqr = 1.0;    // Q_s,0^2 at x=0.01 (GeV^2)
    double e_c = 1.0;
    double anomalous_dimension = 1.0;  // probably want to keep constant
    double alphas_scaling = 1.0;     // C^2 in the expression for alpha_s
    
    
    
    MV ic;
    ic.SetQsqr(qs0sqr);
    ic.SetAnomalousDimension(anomalous_dimension);
    ic.SetE(e_c);       // e_c of MVe parametrization
    ic.SetLambdaQcd(0.241);
    Dipole dipole(&ic);
    
    // Set other configurations
    double maxy = 1.2;      // Solve BK up to this rapidity
    config::RC_LO = config::BALITSKY_LO; // Balitsky running coupling for LO kernel
    config::RESUM_RC = config::RESUM_RC_PARENT; // Parent dipole in the resummation
    config::RESUM_DLOG = true; // Resum doulbe logs
    config::RESUM_SINGLE_LOG = true; // Resum single logs
    config::NF=3;   // Only light quarks
    config::KSUB = 0.6;  // Optimal value for K_sub
    config::NO_K2 = true;  // Do not include numerically demanding full NLO part
    
    // Solve up to some rapidity
    BKSolver solver(&dipole);
    solver.SetAlphasScaling(alphas_scaling);
    solver.Solve(maxy);  // Solve up to maxy
    
    // Give solution to the AmplitudeLib object
    AmplitudeLib DipoleAmplitude(solver.GetDipole()->GetData(), solver.GetDipole()->GetYvals(), solver.GetDipole()->GetRvals());
    
    // Test: print dipole amplitude at rapidities 0,0.5,1
    cout << "# r   N(y=0)  N(y=0.5)   N(y=1)" << endl;
    for (double r=0.0001; r<5; r*=1.1)
    {
        cout << r<< " " << DipoleAmplitude.N(r, 0.01*exp(-0)) << " " << DipoleAmplitude.N(r, 0.01*exp(-0.5)) << " " << DipoleAmplitude.N(r,0.01*exp(-1)) << endl;
    }
    
    
    return 0;
}
