#ifndef SAMPLES_HH
#define SAMPLES_HH


#include "Event.hh"
#include <list>

class Samples
{
	public:
		Samples(int nSamples, double sInterval);
		~Samples();
	private:
		int numSamples;
		double startTime, sampleInterval;
		double *samples;
	public:
		void readEvent(Event *evt, double t);
		double getSample(int i);
		int getNumSamples();
		double getStartTime();
		double getSampleInterval();
		void print();
		void readEvent(double *a, double t);
};

#endif
