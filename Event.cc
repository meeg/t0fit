#include "Event.hh"

#include <gsl/gsl_randist.h>

Event::Event(ShapingCurve *s, gsl_rng *r, double noise)
{
	noiseLevel = noise;
	rng = r;
	shape = s;
	numHits = 0;
}

Event::~Event()
{
}

void Event::addHit(double time, double height)
{
	times[numHits] = time;
	heights[numHits] = height;
	numHits++;
}

double Event::getSignal(double time)
{
	double signal = 0;
	int i;	

	for (i=0; i<numHits; i++)
	{
		signal += heights[i]*shape->getHeight(time-times[i]);
	}

	signal += gsl_ran_gaussian(rng,noiseLevel);
	return signal;
}

double Event::Evaluate(double *x, double *p)
{
	double signal = 0;
	int i;	

	for (i=0; i<numHits; i++)
	{
		signal += heights[i]*shape->getHeight(x[0]-times[i]);
	}
	return signal;
}

int Event::getNumHits()
{
	return numHits;
}

void Event::clear()
{
	numHits = 0;
}

void Event::getTimes(double *t)
{
	int i;
	for (i=0;i<numHits;i++)
		t[i] = times[i];
}

void Event::getHeights(double *h)
{
	int i;
	for (i=0;i<numHits;i++)
		h[i] = heights[i];
}

void Event::sortHits()
{
	int i,mark,min;
	double tmp_t, tmp_a;
	for (mark=0; mark<numHits; mark++)
	{
		min = mark;
		for (i=min+1;i<numHits;i++)
			if (times[i]<times[min]) min=i;
		if (min!=mark)
		{
			tmp_t = times[min];
			tmp_a = heights[min];
			times[min] = times[mark];
			heights[min] = heights[mark];
			times[mark] = tmp_t;
			heights[mark] = tmp_a;
		}
	}
}

void Event::print()
{
	int i;
	printf("Times:\t\t");
	for (i=0;i<numHits;i++)
		printf("%lf\t",times[i]);
	printf("\n");
	printf("Heights:\t");
	for (i=0;i<numHits;i++)
		printf("%lf\t",heights[i]);
	printf("\n");
}
