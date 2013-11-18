/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013
 */

#include "dipole.hpp"
#include "ic.hpp"
#include <iostream>

using std::cerr; using std::cout; using std::endl;


/*
 * Constructor
 * Parameter: Initial condition
 */

Dipole::Dipole(InitialCondition* ic)
{
    // Initialize rvals
    for (double r=1e-7; r<=50; r*=1.1)
    {
        rvals.push_back(r);
    }
    
    yvals.push_back(0);

    dipole_interp=NULL;

}

/*
 * Creates 1D interpolator which gives dipole amplitude as a function
 * of dipole size at given rapidity yvals[yind]
 *
 */

int Dipole::InitializeInterpolation(int yind)
{
    if (yind >= yvals.size())
    {
        cerr << "Rapidity index " << yind << " is too large, maxindex is " << yvals.size()-1 << LINEINFO << endl;
        throw "Too large rapidity index!";
        return -1;
    }

    if (dipole_interp != NULL) // free old interpolator
    {
        dipole_interp -> Clear();
        delete dipole_interp;
    }

    dipole_interp = new Interpolator(rvals, amplitude[yind]);
    dipole_interp->Initialize();

    return 0;


}

/*
 * Evaluates the initialized interpolator for the given dipole size
 * Trhows an execption and returns -1 if interpolator is not
 * initialized
 */
double Dipole::N(double r)
{
    if (dipole_interp==NULL)
    {
        throw "No dipole interpolator...";
        return -1;
    }

    return dipole_interp -> Evaluate(r);

}

/*
 * Add new rapidity
 * The DE solver in BKSolver class obtains the new dipole amplitude
 * at different dipole sizes and saves it
 * Here we get the table containing the dipole amplitudes and the rapidity value
 *
 * Returns -1 if an error occurs, the corresponding rapidity index if
 * succesfull (that is, YVal[returned_vlaue)=y
 */
int Dipole::AddRapidity(double y, double rgrid[])
{

    if (y<=yvals[yvals.size()-1])
    {
        cerr << "Too small rapidity " << y << " at " << LINEINFO << endl;
        return -1;
    }

    yvals.push_back(y);
    std::vector<double> amp;
    for (unsigned int i=0; i<RPoints(); i++)
    {
        amp.push_back(rgrid[i]);
    }

    amplitude.push_back(amp);

    return yvals.size()-1;

}


double Dipole::MinR()
{
    return rvals[0];
}
double Dipole::MaxR()
{
    return rvals[rvals.size()-1];
}

unsigned int Dipole::RPoints()
{
    return rvals.size()-1;
}

double Dipole::RVal(unsigned int rind)
{
    return rvals[rind];
}

double Dipole::YVal(unsigned int yind)
{
    return yvals[yind];
}

Dipole::~Dipole()
{
    if (dipole_interp != NULL)
    {
        dipole_interp->Clear();
        delete dipole_interp;
    }
}

