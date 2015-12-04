#include "FourPoleShapingCurve.hh"

#include <math.h>

FourPoleShapingCurve::FourPoleShapingCurve(double tp1, double tp2) : ShapingCurve(3.0*pow(tp1*pow(tp2,3),0.25))
{
	this->tp1 = tp1;
	this->tp2 = tp2;
    norm = getHeight_intnorm(shapingTime);
}

FourPoleShapingCurve::~FourPoleShapingCurve()
{
}

double FourPoleShapingCurve::getHeight_intnorm(double time)
{
	if (time<0) return 0;
	else return (pow(tp1,2)/(pow(tp1-tp2,3)))*(
                exp((-time)/tp1)-
                exp((-time)/tp2)*(1+
                    (time)*(tp1-tp2)/(tp1*tp2)+
                    pow(((time)*(tp1-tp2)/(tp1*tp2)),2)/2));
}

double FourPoleShapingCurve::getHeight(double time)
{
    return getHeight_intnorm(time)/norm;
}

//double FourShapingCurve::getSlope(double time)
//{
    //return 0.0;
//}
