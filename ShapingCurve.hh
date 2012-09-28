#ifndef SHAPING_CURVE_HH
#define SHAPING_CURVE_HH

class ShapingCurve
{
	public:
		ShapingCurve(double sTime);
		virtual ~ShapingCurve();
	protected:
		double shapingTime;
	public:
		virtual double getHeight(double time); //normalized so peak height = 1
		virtual double getSlope(double time); //used for fitting
		virtual double getPeak(); //used for fitting
		virtual double getSignal(double time, int n, double *par);
};

#endif
