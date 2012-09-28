#include "Fitter.hh"
#include <math.h>


TMinuit *gMinuit = NULL;
Fitter *currentFitter;

Fitter::Fitter(ShapingCurve *sc, int sCount, int pCount, double noise)
{
	nPeaks = pCount;
	shape = sc;
	sigma_noise = noise;



	nSamples = sCount;
	y = (double *) calloc(sCount, sizeof(double));
	sigma = (double *) calloc(sCount, sizeof(double));
	shape = sc;

	fit_par = (double *) calloc(2*pCount, sizeof(double));
	fit_err = (double *) calloc(2*pCount, sizeof(double));
	useSample = (bool *) calloc(sCount, sizeof(bool));
	currentFitter = this;
	verbose = -1;
}

Fitter::~Fitter()
{
	free(y);
	free(sigma);
	free(fit_par);
	free(fit_err);
}

void Fitter::print_fit()
{
	int i;
	printf("Times:\t\t");
	for (i=0;i<nPeaks;i++)
		printf("%lf\t",fit_par[2*i]);
	printf("\n");
	printf("Heights:\t");
	for (i=0;i<nPeaks;i++)
		printf("%lf\t",fit_par[2*i+1]);
	printf("\n");
}

void Fitter::readSamples(Samples *input)
{
	int i;

	for (i=0;i<nSamples;i++)
	{
		y[i] = input->getSample(i);
		sigma[i] = sigma_noise;
		useSample[i] = true;
		//useSample[i] = (y[i]>3.0*sigma_noise);
	}

	startTime = input->getStartTime();
	sampleInterval= input->getSampleInterval();
}

void Fitter::getFitTimes(double *times)
{
	int i;
	for (i=0;i<nPeaks;i++)
		times[i] = fit_par[2*i];
}

void Fitter::getFitHeights(double *heights)
{
	int i;
	for (i=0;i<nPeaks;i++)
		heights[i] = fit_par[2*i+1];
}

void Fitter::getFitPar(double *par)
{
	int i;
	for (i=0;i<2*nPeaks;i++)
		par[i] = fit_par[i];
}

void Fitter::getFitErr(double *err)
{
	int i;
	for (i=0;i<2*nPeaks;i++)
		err[i] = fit_err[i];
}

double Fitter::getChisq(double *par, bool use_all)
{
	int i;
	double Yi, t, delta;
	double chisq = 0;
	for (i = 0; i < nSamples; i++)
	{
		if (!use_all && !useSample[i]) continue;
		Yi = 0;
		t = startTime+sampleInterval*i;
		Yi = shape->getSignal(t,nPeaks,par);
		delta = (Yi - y[i])/sigma[i];
		chisq += delta*delta;
	}
	return chisq;
}

int Fitter::getNumUsedSamples(int start, int end)
{
	int i;
	int count=0;
	for (i=start;i<(end==-1 ? nSamples : end+1);i++)
		if (useSample[i]) count++;
		//if (useSample[i] && startTime+sampleInterval*i>getFirstPeak()) count++;
	return count;
}

double Fitter::getFirstPeak()
{
	double firstPeak = fit_par[0];
	for (int i=1;i<nPeaks;i++)
		if (fit_par[2*i]<firstPeak) firstPeak = fit_par[2*i];
	return firstPeak;
}

double Fitter::getFitChisq()
{
	return getChisq(fit_par);
}

ShapingCurve * Fitter::getShape()
{
	return shape;
}

int Fitter::getDOF()
{
	return dof;
}

int Fitter::getStatus()
{
	return status;
}

void Fitter::sortFit()
{
	int i,mark,min;
	double tmp_t, tmp_a;
	for (mark=0; mark<nPeaks; mark++)
	{
		min = mark;
		for (i=min+1;i<nPeaks;i++)
			if (fit_par[2*i]<fit_par[2*min]) min=i;
		if (min!=mark)
		{
			tmp_t = fit_par[2*min];
			tmp_a = fit_par[2*min+1];
			fit_par[2*min] = fit_par[2*mark];
			fit_par[2*min+1] = fit_par[2*mark+1];
			fit_par[2*mark] = tmp_t;
			fit_par[2*mark+1] = tmp_a;
		}
	}
}

void Fitter::printResiduals()
{
	int i;

	for (i = 0; i < nSamples; i++)
	{
		double t = startTime+sampleInterval*i;
		printf("%lf\t", getSignal(t) - y[i]);
	}
	printf("\n");
}

double Fitter::getSignal(double time)
{
	int i;
	double Yi = 0;

	for (i=0; i<nPeaks; i++)
	{
		Yi += fit_par[2*i+1] * shape->getHeight(time-fit_par[2*i]);
	}
	return Yi;
}

void Fitter::setVerbosity(int verbosity)
{
	verbose = verbosity;
}

int Fitter::getNumPeaks()
{
	return nPeaks;
}

void Fitter::setSigmaNoise(double sigma)
{
	sigma_noise = sigma;
}

Double_t evtCurve(Double_t *x, Double_t *par)
{
	Fitter *theFitter = currentFitter;
	return (theFitter->getShape())->getSignal(x[0],theFitter->getNumPeaks(),par);
}

