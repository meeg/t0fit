#include "Samples.hh"

Samples::Samples(int nSamples, double sInterval)
{
	numSamples = nSamples;
	sampleInterval = sInterval;
	samples = (double *) calloc(nSamples, sizeof(double));
}

Samples::~Samples()
{
	free(samples);
}

void Samples::readEvent(Event *evt, double t)
{
	int i;
	startTime = t;
	for (i=0;i<numSamples;i++)
	{
		samples[i] = evt->getSignal(t+sampleInterval*i);
	}
}

double Samples::getSample(int i)
{
	return samples[i];
}

int Samples::getNumSamples()
{
	return numSamples;
}

double Samples::getStartTime()
{
	return startTime;
}

double Samples::getSampleInterval()
{
	return sampleInterval;
}

void Samples::print()
{
	int i;
	for (i=0; i<numSamples; i++)
		printf("%lf\t",samples[i]);
	printf("\n");
}

void Samples::readEvent(double *a, double t)
{
	int i;
	startTime = t;
	for (i=0;i<numSamples;i++)
	{
		samples[i] = a[i];
	}
}
