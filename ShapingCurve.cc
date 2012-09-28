#include "ShapingCurve.hh"

#include <math.h>

ShapingCurve::ShapingCurve(double sTime)
{
	shapingTime = sTime;
}

ShapingCurve::~ShapingCurve()
{
}

double ShapingCurve::getHeight(double time)
{
	if (time<0) return 0;
	else return (time/shapingTime)*exp(1.0-time/shapingTime);
}

double ShapingCurve::getSlope(double time)
{
	if (time<0) return 0;
	//if (time<-0.2*shapingTime) return 0;
	//else if (time<0) return (1.0+time/(0.2*shapingTime))*exp(1.0)/shapingTime;
	else return (1.0-time/shapingTime)*exp(1.0-time/shapingTime)/shapingTime;
}

double ShapingCurve::getPeak()
{
	return shapingTime;
}

double ShapingCurve::getSignal(double time, int n, double *par)
{
	double signal = 0;
	int i;	

	for (i=0; i<n; i++)
	{
		signal += par[2*i+1]*getHeight(time-par[2*i]);
	}

	return signal;
}

