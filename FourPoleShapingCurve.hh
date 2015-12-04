#ifndef FOUR_POLE_SHAPING_CURVE_HH
#define FOUR_POLE_SHAPING_CURVE_HH
#include "ShapingCurve.hh"

class FourPoleShapingCurve : public ShapingCurve
{
	public:
		FourPoleShapingCurve(double tp1, double tp2);
		virtual ~FourPoleShapingCurve();
	private:
		double tp1,tp2;
        double norm; //peak amplitude of the integral-normalized function
        double getHeight_intnorm(double time);
	public:
		double getHeight(double time); //normalized so peak height = 1
		//double getSlope(double time); //used for fitting
		//virtual double getPeak(); //used for fitting
		//virtual double getSignal(double time, int n, double *par);
};

#endif
