/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013
 */

#ifndef _NLOBK_DIPOLE_H
#define _NLOBK_DIPOLE_H

#include "ic.hpp"
#include "nlobk_config.hpp"
#include <vector>
#include <tools/interpolation.hpp>
#include <string>

/*
 * Dipole amplitude
 * Initially initialize this by InitialCondition, the solver
 * fills at larger rapidities
 */


class Dipole
{
    public:
        Dipole(InitialCondition* ic_);             // Initialize from a given ic
        ~Dipole();
        int InitializeInterpolation(int yind);
            // Create interpolator of dipole amplitude values at rapidity yvals[i]

        double N(double r);         // Evaluates the initialized interpolator at r
        double S(double r);         // 1 - N

        int AddRapidity(double y, double rgrid[]);

        double MinR();
        double MaxR();
        unsigned int RPoints();
        unsigned int YPoints();
        double RVal(unsigned int rind);
        double YVal(unsigned int yind);

        int Save(std::string filename);

    private:
        // amplitude[i][j] is vector of dipole amplitude values at rapidity yvals[i]
        // at dipole size rvals[j] 
        std::vector< std::vector<double > > amplitude;   
        std::vector< double > yvals; 
        std::vector< double > rvals;
        InitialCondition* ic;

        Interpolator* dipole_interp;        // Initialized interpolator to evaluate N(r)
        unsigned int interpolator_yind;     // Rapidity index at which the interpolator is initialized
};

const std::string VERSION = "0.0-dev";
#endif
