/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013
 */

#include "dipole.hpp"
#include "ic.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>

using std::cerr; using std::cout; using std::endl;


/*
 * Constructor
 * Parameter: Initial condition
 */

Dipole::Dipole(InitialCondition* ic)
{
    // Initialize rvals
    std::vector<double> initial_amplitude;
    for (double r=1e-4; r<=40; r*=1.6)
    {
        rvals.push_back(r);
        initial_amplitude.push_back( ic->DipoleAmplitude(r) );
    }
    amplitude.push_back(initial_amplitude);
    
    yvals.push_back(0);

    // Initialize initial condition interpolator
    dipole_interp=NULL;
    InitializeInterpolation(0);
    cout << "Dipole is ready!" << endl;
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
        return -1;
    }

    if (dipole_interp != NULL) // free old interpolator
    {
        delete dipole_interp;
    }

    dipole_interp = new Interpolator(rvals, amplitude[yind]);
    dipole_interp->SetUnderflow(0);
    dipole_interp->SetOverflow(1.0);
    dipole_interp->SetFreeze(true);
    dipole_interp->Initialize();
    interpolator_yind = yind;

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
        return -1;
    }

    return dipole_interp -> Evaluate(r);
}

double Dipole::S(double r)
{
    return 1.0-N(r);
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

/*
 * Save amplitude to the given file
 *
 * Syntax is the same as in rbk and AmplitudeLib
 */
int Dipole::Save(std::string filename)
{
    std::ofstream out;
    out.open(filename.c_str());

    out << "###" << std::scientific << std::setprecision(15) << MinR() << endl;
    out << "###" << std::scientific << std::setprecision(15) <<
        1.1  << endl;   // rmultiplier
    out << "###" << RPoints() << endl;
    out << "###0.01" << endl;   // x0

    for (int yind=0; yind<YPoints(); yind++)
    {
        out << "###" << std::scientific << std::setprecision(15)
            << yvals[yind] << endl;
        for (int rind=0; rind<RPoints(); rind++)
        {
            out << std::scientific << std::setprecision(15)
                << amplitude[yind][rind] << endl;
        }
    }
    out.close();
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

unsigned int Dipole::YPoints()
{
    return yvals.size()-1;
}

double Dipole::RVal(unsigned int rind)
{
    if (rind < 0 or rind>rvals.size()-1)
    {
        cerr << "Invalid dipole size index " << rind << " at " << LINEINFO << endl;
        return 0;
    }
    return rvals[rind];
}

double Dipole::YVal(unsigned int yind)
{
    return yvals[yind];
}

Dipole::~Dipole()
{
    cout << "Destroying dipole..." << endl;
    if (dipole_interp != NULL)
    {
        delete dipole_interp;
    }
}

