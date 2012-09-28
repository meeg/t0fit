#ifndef SMOOTH_SHAPING_CURVE_HH
#define SMOOTH_SHAPING_CURVE_HH
#include "ShapingCurve.hh"
#include "TSpline.h"

class SmoothShapingCurve : public ShapingCurve
{
	private:
		TSpline3 *theSpline;
	public:
		SmoothShapingCurve(double sTime);
		SmoothShapingCurve(int ni, double *ti, double *yi);
		~SmoothShapingCurve();
		double getHeight(double time); //normalized so peak height = 1
		double getSlope(double time); //used for fitting
};

#endif
