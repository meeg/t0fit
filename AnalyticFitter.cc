#include "AnalyticFitter.hh"
#include <math.h>
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH2F.h"


AnalyticFitter::AnalyticFitter(ShapingCurve *sc, int sCount, int pCount, double noise) : Fitter(sc, sCount, pCount, noise)
{
	shapingTime = sc->getPeak();
	sectionStart = (int *) calloc(pCount, sizeof(int));
	sectionEnd = (int *) calloc(pCount, sizeof(int));
	section_par = (double *) calloc(2*pCount, sizeof(double));
	section_err = (double *) calloc(2*pCount, sizeof(double));
	section_covar = (double *) calloc(pCount, sizeof(double));
}

AnalyticFitter::~AnalyticFitter()
{
	free(sectionStart);
	free(sectionEnd);
	free(section_par);
	free(section_err);
	free(section_covar);
}

int AnalyticFitter::doAnalyticFit(int start, int end, double &time, double &time_err, double &height, double &height_err, double &covar)
{
	int i;
	int length = 0;
	int *k = (int *) calloc(nSamples, sizeof(int));
	for (i=start;i<end+1;i++) if (useSample[i])
	{
		k[length] = i;
		length++;
	}

	if (length<1) return -1; //not enough samples
	if (length == 1) //TODO: do something more intelligent here
	{
		time = startTime+sampleInterval*(k[0]-1);
		height = y[k[0]]*shapingTime*exp(sampleInterval/shapingTime-1.0)/sampleInterval;
		time_err = sampleInterval;
		height_err = height;
		covar = 0;
		if (verbose>2) printf("start %d, end %d, time %f, time_err %f, height %f, height_err %f\n",start, end, time, time_err, height, height_err);
		return 0;
	}
	double *a = (double *) calloc(length, sizeof(double));
	double *t = (double *) calloc(length, sizeof(double));
	double *p = (double *) calloc(length, sizeof(double));
	for (i=0;i<length;i++)
	{
		p[i] = y[k[i]]/sigma[k[i]];
		t[i] = startTime+sampleInterval*(k[i]);
		a[i] = exp(1-t[i]/shapingTime)/(shapingTime*sigma[k[i]]);
	}
	double pa, aatt, pat, aat, aa;
	pa = 0;
	aatt = 0;
	pat = 0;
	aat = 0;
	aa = 0;
	for (i=0;i<length;i++)
	{
		pa += p[i]*a[i];
		aatt += a[i]*a[i]*t[i]*t[i];
		pat += p[i]*a[i]*t[i];
		aat += a[i]*a[i]*t[i];
		aa += a[i]*a[i];
	}

	time = (pa*aatt-pat*aat)/(pa*aat-aa*pat);
	height = pa/(exp(time/shapingTime)*(aat-time*aa));
	double time_var = 0;
	double height_var = 0;
	covar = 0;
	double *dt_dp = (double *) calloc(length, sizeof(double));
	double *dh_dp = (double *) calloc(length, sizeof(double));
	for (i=0;i<length;i++)
	{
		dt_dp[i] = a[i]*(aatt-t[i]*aat-time*(aat-t[i]*aa))/(pa*aat-aa*pat);
		dh_dp[i] = (a[i]*exp(-1.0*time/shapingTime)+height*dt_dp[i]*aa)/(aat-time*aa)-height*dt_dp[i]/shapingTime;
		time_var += dt_dp[i]*dt_dp[i];
		height_var += dh_dp[i]*dh_dp[i];
		covar += dt_dp[i]*dh_dp[i];
	}
	time_err = sqrt(time_var);
	height_err = sqrt(height_var);

	if (verbose>2) printf("start %d, end %d, time %f, time_err %f, height %f, height_err %f\n",start, end, time, time_err, height, height_err);
	
	free(a);
	free(t);
	free(p);
	free(k);
	free(dt_dp);
	free(dh_dp);

	return 0;
}

void AnalyticFitter::doFitSingle()
{
	int i;
	int start, end;
	int best_start, best_end;
	double chisq, best_chisq;
	best_chisq = -1;
	if (nPeaks>1) for (i=2;i<2*nPeaks;i++)
	{
		fit_par[i] = 0;
		fit_err[i] = 0;
	}
	nSections = 1;
	for (start=0;start<nSamples-1;start++)
		for (end=start+1;end<nSamples;end++)
		{
			status = doAnalyticFit(start, end, section_par[0], section_err[0], section_par[1], section_err[1], section_covar[0]);
			if (status!=0) continue;
			mergeSections();
			chisq = getChisq(fit_par);
			if (verbose>2) printf("start %d, end %d, total chisq %f\n",start,end,chisq);
			if (best_chisq == -1 || chisq<best_chisq)
			{
				best_start = start;
				best_end = end;
				best_chisq = chisq;
			}
		}

	/*
	for (i=n-3;i>=0;i--)
	{
		status = doAnalyticFit(i, n-1, fit_par[0], fit_err[0], fit_par[1], fit_err[1]);
		if (status!=0) continue;
		chisq = memberFCN(fit_par);
		if (verbose>2) printf("start %d, total chisq %f\n",i,chisq);
		if (best_chisq == -1 || chisq<best_chisq)
		{
			best = i;
			best_chisq = chisq;
		}

		   //if (chisq/(n-i-2)>10)
		   //{
		   //doAnalyticFit(i+1, n-1, fit_par[0], fit_err[0], fit_par[1], fit_err[1], chisq);
		   //break;
		   //}
	}
	*/

	//status = chisq/(n-i-2)<10.0 ? 0 : -1;
	if (verbose>0) printf("Good range is %d to %d\n",best_start,best_end);

	if (best_chisq != -1)
	{
		status = 0;
		doAnalyticFit(best_start, best_end, fit_par[0], fit_err[0], fit_par[1], fit_err[1], section_covar[0]);
		dof = getNumUsedSamples(best_start, best_end)-2;
	}
	else status = 1;
}

void AnalyticFitter::doFit()
{
	int i,j;
	int *start = (int *) calloc(nPeaks, sizeof(int));
	int *best_start = (int *) calloc(nPeaks, sizeof(int));
	double chisq, best_chisq;
	status=0;
	best_chisq = -1;
	nSections = nPeaks;
	if (nPeaks>1) for (i=2;i<2*nPeaks;i++)
	{
		fit_par[i] = 0;
		fit_err[i] = 0;
	}

	for (i=0;i<nPeaks;i++)
	{
		start[i] = i;
	}
	while (start[0]<nSamples-nPeaks+1)
	{
		for (i=0;i<nPeaks-1;i++)
		{
			status = doAnalyticFit(start[i], start[i+1]-1, section_par[2*i], section_err[2*i], section_par[2*i+1], section_err[2*i+1], section_covar[i]);
			if (status!=0) break;
		}
		if (status==0) status = doAnalyticFit(start[nPeaks-1], nSamples-1, section_par[2*(nPeaks-1)], section_err[2*(nPeaks-1)], section_par[2*(nPeaks-1)+1], section_err[2*(nPeaks-1)+1], section_covar[nPeaks-1]);
		if (status==0)
		{
			mergeSections();
			chisq = getChisq(fit_par);
			if (verbose>2)
			{
				for (j=0;j<nPeaks;j++) printf("%d, ",start[j]);
				printf("total chisq %f\n",chisq);
			}
			if (best_chisq<0 || chisq<best_chisq)
			{
				memcpy(best_start,start,(nPeaks)*sizeof(int));
				best_chisq = chisq;
			}
		}

		start[nPeaks-1]++;
		for (i=nPeaks-1;i>0 && start[i]==nSamples;i--)
		{
			start[i-1]++;
			start[i]=start[i-1]+1;
		}
	}

	if (best_chisq != -1)
	{
		if (verbose>2)
		{
			printf("Best: ");
			for (j=0;j<nPeaks;j++) printf("%d, ",best_start[j]);
			printf("total chisq %f\n",best_chisq);
		}
		for (i=0;i<nPeaks-1;i++)
		{
			doAnalyticFit(best_start[i], best_start[i+1]-1, section_par[2*i], section_err[2*i], section_par[2*i+1], section_err[2*i+1], section_covar[i]);
		}
		doAnalyticFit(best_start[nPeaks-1], nSamples-1, section_par[2*(nPeaks-1)], section_err[2*(nPeaks-1)], section_par[2*(nPeaks-1)+1], section_err[2*(nPeaks-1)+1], section_covar[nPeaks-1]);
		mergeSections();
		dof = getNumUsedSamples()-2*nPeaks;
		if (dof<0) dof = 0;
	}
	else status = 1;
}

void AnalyticFitter::mergeSections()
{
	int i;
	double exp_1, exp_s;
	double t_1, t_2, t_s;
	double h_1, h_2, h_s;

	fit_par[0] = section_par[0];
	fit_par[1] = section_par[1];
	fit_err[0] = section_err[0];
	fit_err[1] = section_err[1];

	for (i=1;i<nSections;i++)
	{
		t_1 = fit_par[2*(i-1)];
		h_1 = fit_par[2*(i-1)+1];
		t_s = section_par[2*i];
		h_s = section_par[2*i+1];
		exp_1 = exp(t_1/shapingTime);
		exp_s = exp(t_s/shapingTime);
		t_2 = (h_s*t_s*exp_s-h_1*t_1*exp_1)/(h_s*exp_s-h_1*exp_1);
		h_2 = (h_s*(t_2-t_s+shapingTime)*exp_s-h_1*(t_2-t_1+shapingTime)*exp_1)*exp(-1.0*t_2/shapingTime)/shapingTime;

		//printf("%f %f\n",h_s*(t_2-t_s)*exp_s,h_1*(t_2-t_1)*exp_1);
		//printf("%f %f\n",h_s*shape->getHeight(t_2-t_s),h_1*shape->getHeight(t_2-t_1));

		fit_par[2*i] = t_2;
		fit_par[2*i+1] = h_2;
		if (verbose>2) printf("merged: peak %d has time %f, height %f\n",i,t_2, h_2);
		fit_err[2*i] = 0; //TODO: actually propagate these errors
		fit_err[2*i+1] = 0;
	}

}

void AnalyticFitter::plotFit(Event *evt, char *name)
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
	for (i=0;i<nPeaks;i++)
	{
		func->FixParameter(2*i,fit_par[2*i]);
		func->FixParameter(2*i+1,fit_par[2*i+1]);
	}

	func->SetLineColor(2);
	func->SetLineWidth(1);
	func->DrawCopy("CSAME");

	/*
	   double t_err, h_err, covar;
	   func->SetLineWidth(1);

	   for (i=n-3;i>=0;i--)
	   {
	   printf("%d\n",i);
	   doAnalyticFit(i, n-1, t[0], t_err, h[0], h_err, covar);
	   func->FixParameter(0,t[0]);
	   func->FixParameter(1,h[0]);
	   func->SetLineColor(n-i-1);
	   func->DrawCopy("CSAME");
	   }
	   */


	delete func;


	c1->SaveAs(name);

	delete c1;
}
