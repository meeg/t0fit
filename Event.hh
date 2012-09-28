#ifndef EVENT_HH
#define EVENT_HH


#include <gsl/gsl_rng.h>
#include "ShapingCurve.hh"

class Event
{
	public:
		Event(ShapingCurve *s, gsl_rng *r, double noise);
		~Event();
	private:
		int numHits;
		double noiseLevel;
		ShapingCurve *shape;
		double times[10], heights[10];
		gsl_rng *rng;
	public:
		void addHit(double time, double height);
		double getSignal(double time); //return a sample at the given time, with Gaussian noise (not deterministic)
		double Evaluate(double *x, double *p); //no-noise getSignal for use with TF1 constructor
		int getNumHits();
		void clear(); //clear all hits
		void getTimes(double *times);
		void getHeights(double *heights);
		void sortHits(); //sorts hits by increasing time
		void print(); //print times and heights of hits
};

#endif
