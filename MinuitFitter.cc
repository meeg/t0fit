#include "MinuitFitter.hh"

#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH2F.h"
#include <math.h>

MinuitFitter::MinuitFitter(ShapingCurve *sc, int sCount, int pCount, double noise) : Fitter(sc, sCount, pCount, noise)
{
	if (gMinuit==NULL) gMinuit = new TMinuit();
	guess_par = (double *) calloc(2*pCount, sizeof(double));
	tweakFCN = true;
}

MinuitFitter::~MinuitFitter()
{
	free(guess_par);
}

void MinuitFitter::doMinuitFit()
{
	initMinuit();
	double arglist[10];
	int ierflg = 0;
	int i;
	for (i=0; i<nPeaks; i++)
	{
		gMinuit->DefineParameter(2*i,"time",guess_par[2*i],25.0,-100.0,200.0);
		gMinuit->DefineParameter(2*i+1,"height",guess_par[2*i+1],1.0*guess_par[2*i+1],0.0,50.0);
	}

	arglist[0] = 1000;
	gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);
	status = ierflg;
	if (false && status == 0)
	{
		arglist[0] = 1000;
		gMinuit->mnexcm("IMPROVE", arglist ,1,ierflg);
	}
	for (i=0; i<2*nPeaks; i++)
	{
		gMinuit->GetParameter(i,fit_par[i],fit_err[i]);
	}
	dof = nSamples-2*nPeaks;

	if (tweakFCN)
	{
		tweakFCN = false;
		for (i=0; i<nPeaks; i++)
		{
			gMinuit->DefineParameter(2*i,"time",fit_par[2*i],5.0,-100.0,200.0);
			gMinuit->DefineParameter(2*i+1,"height",fit_par[2*i+1],0.05*fit_par[2*i+1],0.0,50.0);
		}
		arglist[0] = 1000;
		gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);
		status = ierflg;
		if (status == 0)
		{
			arglist[0] = 1000;
			gMinuit->mnexcm("IMPROVE", arglist ,1,ierflg);
		}
		for (i=0; i<2*nPeaks; i++)
		{
			gMinuit->GetParameter(i,fit_par[i],fit_err[i]);
		}

		tweakFCN = true;
	}
}

void MinuitFitter::doFit()
{
	makeGuess();
	doMinuitFit();
}

double MinuitFitter::memberFCN(double *par)
{
	size_t i,j;
	double delta, delta2;
	double chisq = 0;
	double Yi;
	double b;
	int firstPeak;

	for (i = 0; i < nSamples; i++)
	{
		if (!useSample[i]) continue;
		Yi = 0;
		double t = startTime+sampleInterval*i;
		Yi = shape->getSignal(t,nPeaks,par);
		delta = (Yi - y[i]);
		if (tweakFCN && Yi<2*sigma[i]) //if sample is before all the t0's, use the distance to the leading edge instead of the residual
		{
			firstPeak = -1;
			for (j=0;j<nPeaks;j++)
			{
				if (par[2*j]+10.0>t && par[2*j]-20.0<t && (firstPeak==-1 || par[2*j]<par[2*firstPeak]))
					firstPeak = j;
			}
			if (firstPeak!=-1)
			{
				b = par[2*firstPeak+1]*shape->getSlope(0); //slope of leading edge
				delta2 = (y[i]-b*(t-par[2*firstPeak]))/sqrt(1+b*b);
				if (delta2>0 && fabs(delta2) < fabs(delta))
					delta = (0.0*fabs(delta)+1.0*fabs(delta2));
			}
		}
		delta /= sigma[i];

		chisq += delta*delta;
	}

	return chisq/(nSamples-2*nPeaks);
}

void MinuitFitter::memberGradFCN(double *par, double *gin)
{
	size_t i,j;
	double delta;
	double Yi;
	double dchisq = 0;

	for (i=0;i<2*nPeaks;i++)
		gin[i] = 0;

	for (i = 0; i < nSamples; i++)
	{
		if (!useSample[i]) continue;
		double t = startTime+sampleInterval*i;
		Yi = 0;
		for (j=0; j<nPeaks; j++)
		{
			Yi += par[2*j+1] * shape->getHeight(t-par[2*j]);
		}
		for (j=0; j<nPeaks; j++)
		{
			delta = (Yi - y[i])/(sigma[i]*sigma[i]);
			gin[2*j] += -2.0*delta*par[2*j+1]*shape->getSlope(t-par[2*j]);
			gin[2*j+1] += 2.0*delta*shape->getHeight(t-par[2*j]);
		}
	}
}

void MinuitFitter::initMinuit()
{
	currentFitter = this;
	gMinuit->SetFCN(FCN);
	//gMinuit->SetPrintLevel(-1);
	gMinuit->Command("CLE");

	gMinuit->Command("SET NOW");
	//gMinuit->Command("SET WAR");
	gMinuit->Command("SET STR 2");
	//gMinuit->Command("SET GRA 1");
	gMinuit->Command("SET ERR 1");
}

void MinuitFitter::setTweakFCN(bool x)
{
	tweakFCN = x;
}

bool MinuitFitter::getTweakFCN()
{
	return tweakFCN;
}

void MinuitFitter::plotFit(Event *evt, char *name)
{
	currentFitter = this;
	double *t, *h;
	int i,j;
	TGraph *graph;
	TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,750);
	TF1 *func;

	TH2F *hpx = new TH2F("hpx","",10,-25.0,175.0,10,-5.0,30.0*nPeaks); // axis range
	hpx->SetStats(kFALSE); // no statistics
	hpx->DrawCopy();
	delete hpx;

	t = (double *) calloc(nSamples, sizeof(double));
	h = (double *) calloc(nSamples, sizeof(double));
	for (i=0;i<nSamples;i++)
		t[i] = startTime+i*sampleInterval;
	graph = new TGraph(nSamples,t,y);
	graph->DrawClone("*");
	delete graph;

	j=0;
	for (i=0;i<nSamples;i++)
	{
		if (useSample[i])
		{
			t[j] = startTime+i*sampleInterval;
			h[j] = y[i];
			j++;
		}
	}
	graph = new TGraph(j,t,h);
	graph->SetMarkerColor(2);
	graph->DrawClone("*");
	delete graph;
	free(t);
	free(h);

	func = new TF1("curve",evt,&Event::Evaluate,-100.0,200.0,0,"Event");
	func->SetLineWidth(1);
	func->SetNpx(1000);
	func->DrawCopy("CSAME");
	delete func;

	func = new TF1("curve",evtCurve,-100.0,200.0,2*nPeaks);
	t = (double *) calloc(nPeaks, sizeof(double));
	h = (double *) calloc(nPeaks, sizeof(double));
	getFitTimes(t);
	getFitHeights(h);
	for (i=0;i<nPeaks;i++)
	{
		func->FixParameter(2*i,t[i]);
		func->FixParameter(2*i+1,h[i]);
	}
	func->SetLineColor(2);
	func->DrawCopy("CSAME");
	delete func;

	char funcname[100];

	for (j=0;j<nPeaks;j++)
	{
		sprintf(funcname,"curve %d",2*j);
		func = new TF1(funcname,FCNCurve,-100.0,200.0,2*nPeaks+1);
		func->SetLineWidth(1);
		func->SetNpx(1000);
		for (i=0;i<2*nPeaks;i++)
		{
			func->FixParameter(i,fit_par[i]);
		}
		func->FixParameter(2*nPeaks,(double)j);
		func->SetLineColor(2*j+3);
		func->DrawCopy("CSAME");
		delete func;

		sprintf(funcname,"curve %d",2*j+1);
		func = new TF1(funcname,FCNCurve2,-100.0,200.0,2*nPeaks+1);
		func->SetLineWidth(1);
		func->SetNpx(1000);
		for (i=0;i<2*nPeaks;i++)
		{
			func->FixParameter(i,fit_par[i]);
		}
		func->FixParameter(2*nPeaks,(double)j);
		func->SetLineColor(2*j+4);
		func->DrawCopy("CSAME");
		delete func;
	}

	c1->SaveAs(name);

	delete c1;
}

void MinuitFitter::setVerbosity(int verbosity)
{
	gMinuit->SetPrintLevel(verbosity);
	verbose = verbosity;
}

void MinuitFitter::makeGuess()
{
	int i,j;

	double min, max;

	double *dy = (double *) calloc(nSamples, sizeof(double));

	dy[0] = y[0];
	for (i=0;i<nSamples-1;i++) dy[i+1]=y[i+1]-y[i];

	min = -1;
	for (j=0;j<nPeaks;j++)
	{
		max = 0;
		for (i=0;i<nSamples;i++)
		{
			if ((min<0 || dy[i]<min) && dy[i]>max)
			{
				max = dy[i];
			}
		}
		min = max;
	}

	j=0;
	for (i=0;i<nSamples;i++)
	{
		if (dy[i]>=min)
		{
			if (j==nPeaks)
			{
				if (verbose>2) printf("bad makeGuess; min = %f\n",min);
				break;
			}
			guess_par[2*j] = startTime+sampleInterval*i;
			guess_par[2*j+1] = dy[i];
			j++;
		}
	}
	free(dy);
}

void MinuitFitter::print_guess()
{
	int i;
	printf("Times:\t\t");
	for (i=0;i<nPeaks;i++)
		printf("%lf\t",guess_par[2*i]);
	printf("\n");
	printf("Heights:\t");
	for (i=0;i<nPeaks;i++)
		printf("%lf\t",guess_par[2*i+1]);
	printf("\n");
}

void MinuitFitter::setGuess(int i, double time, double height)
{
	guess_par[2*i] = time;
	guess_par[2*i+1] = height;
}

double MinuitFitter::guessChisq()
{
	return getChisq(guess_par);
}

void FCN(Int_t&npar, Double_t*gin, Double_t&f, Double_t*par, Int_t flag)
{
	MinuitFitter *theFitter = (MinuitFitter *)currentFitter;
	switch(flag)
	{
		case 1:
			// Initialization
			break;
		case 2:
			//	theFitter->memberGradFCN(par,gin);
			// Compute derivatives
			// store them in gin
			break;
		case 3:
			// after the fit is finished
			break;
		default:
			f = theFitter->memberFCN(par);
			break;
	}
}

Double_t FCNCurve(Double_t *x, Double_t *par)
{
	MinuitFitter *theFitter = (MinuitFitter *)currentFitter;
	bool temp = theFitter->getTweakFCN();
	theFitter->setTweakFCN(true);
	int i, n;
	n = theFitter->getNumPeaks();
	double *newpar = (double *) calloc(2*n, sizeof(double));
	int k = (int)par[2*n];
	for (i=0;i<2*n;i++)
		newpar[i] = par[i];
	newpar[2*k] = x[0];
	double y = theFitter->memberFCN(newpar);
	free(newpar);
	theFitter->setTweakFCN(temp);
	return y;
}

Double_t FCNCurve2(Double_t *x, Double_t *par)
{
	MinuitFitter *theFitter = (MinuitFitter *)currentFitter;
	bool temp = theFitter->getTweakFCN();
	theFitter->setTweakFCN(false);
	int i, n;
	n = theFitter->getNumPeaks();
	double *newpar = (double *) calloc(2*n, sizeof(double));
	int k = (int)par[2*n];
	for (i=0;i<2*n;i++)
		newpar[i] = par[i];
	newpar[2*k] = x[0];
	double y = theFitter->memberFCN(newpar);
	free(newpar);
	theFitter->setTweakFCN(temp);
	return y;
}

