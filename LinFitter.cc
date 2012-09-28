#include "LinFitter.hh"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TH2F.h"
#include "TStyle.h"
#include <math.h>


LinFitter::LinFitter(ShapingCurve *sc, int sCount, int pCount, double noise) : Fitter(sc, sCount, pCount, noise)
{
	if (gMinuit==NULL) gMinuit = new TMinuit();
	guess_t = (double *) calloc(pCount, sizeof(double));
	coeff_mat = gsl_matrix_alloc(pCount,pCount);
	a_vec = gsl_vector_alloc(pCount);
	fit_t = (double *) calloc(pCount, sizeof(double));
	fit_t_err = (double *) calloc(pCount, sizeof(double));
	/*
	v = 
	TensorProd(v,v);
	*/
}

LinFitter::~LinFitter()
{
	free(guess_t);
	gsl_matrix_free(coeff_mat);
	gsl_vector_free(a_vec);
	free(fit_t);
	free(fit_t_err);
}

double LinFitter::doLinFit(double *times, double *par)
{
	int i,j;
	int status;
	double t;
	double norm;


	int length = 0;
	int *k = (int *) calloc(nSamples, sizeof(int));
	for (i=0;i<nSamples;i++) if (useSample[i])
	{
		k[length] = i;
		length++;
	}

	y_vec = gsl_vector_alloc(length);
	sc_mat = gsl_matrix_alloc(nPeaks,length);

	for (i=0;i<length;i++)
		gsl_vector_set(y_vec,i,y[k[i]]/sigma[k[i]]);
	//printf("y_vec:\n");
	//gsl_vector_fprintf(stdout,y_vec,"%f");
	for (i=0;i<nPeaks;i++) for (j=0;j<length;j++)
	{
		t = startTime+sampleInterval*k[j];
		gsl_matrix_set(sc_mat,i,j,shape->getHeight(t-times[i])/sigma[k[j]]);
	}
	//printf("sc_mat:\n");
	//gsl_matrix_fprintf(stdout,sc_mat,"%f");
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,sc_mat,sc_mat,0.0,coeff_mat); //coeff_mat = sc_mat*transpose(sc_mat)
	//printf("coeff_mat:\n");
	//gsl_matrix_fprintf(stdout,coeff_mat,"%f");
	gsl_blas_dgemv(CblasNoTrans,1.0,sc_mat,y_vec,0.0,a_vec); //a_vec = sc_mat*y_vec
	//printf("a_vec:\n");
	//gsl_vector_fprintf(stdout,a_vec,"%f");
	gsl_set_error_handler_off();
	status = gsl_linalg_HH_svx(coeff_mat,a_vec); //solve the matrix equation coeff_mat*x=a_vec, and write the result in a_vec
	if (status==GSL_ESING || gsl_vector_min(a_vec)<0) //if we couldn't solve the matrix equation, or one of the peaks has negative height, give up
	{
		//	printf("fail!\n");
		//	for (i=0;i<nPeaks;i++)
		//		printf("%f\t",times[i]);
		//	printf("\n");
		gsl_vector_set_zero(a_vec); //set all peak heights to 0
	}
	if (par!=NULL) for (i=0;i<nPeaks;i++) //fill par with fitted values
	{
		par[2*i] = times[i];
		par[2*i+1] = gsl_vector_get(a_vec,i);
	}
	gsl_blas_dgemv(CblasTrans,-1.0,sc_mat,a_vec,1.0,y_vec); //subtract fit from samples
	norm = gsl_blas_dnrm2(y_vec);
	gsl_vector_free(y_vec);
	gsl_matrix_free(sc_mat);
	free(k);
	return norm*norm;
}

void LinFitter::initMinuit()
{
	currentFitter = this;
	gMinuit->SetFCN(LinFCN);
	//gMinuit->SetPrintLevel(-1);
	gMinuit->Command("CLE");

	gMinuit->Command("SET NOW");
	//gMinuit->Command("SET WAR");
	gMinuit->Command("SET STR 2");
	gMinuit->Command("SET ERR 1");
}

void LinFitter::doMinuitFit()
{
	doGridFit();
	initMinuit();
	double arglist[10];
	int ierflg = 0;
	int i;
	for (i=0; i<nPeaks; i++)
	{
		gMinuit->DefineParameter(i,"time",guess_t[i],25.0,-200.0,120.0);
	}

	arglist[0] = 1000;
	gMinuit->mnexcm("SIMPLEX", arglist ,1,ierflg);
	status = ierflg;
	if (status == 0)
	{
		arglist[0] = 1000;
		gMinuit->mnexcm("IMPROVE", arglist ,1,ierflg);
	}
	for (i=0; i<nPeaks; i++)
	{
		gMinuit->GetParameter(i,fit_t[i],fit_t_err[i]);
	}

	doLinFit(fit_t,fit_par);
}

void LinFitter::plotFit(Event *evt, const char *name)
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
	func->SetLineWidth(1);
	func->SetNpx(1000);
	func->SetLineColor(2);
	func->DrawCopy("CSAME");
	delete func;

	char funcname[100];

	for (j=0;j<nPeaks;j++)
	{
		sprintf(funcname,"curve %d",j);
		func = new TF1(funcname,LinFCNCurve,-100.0,200.0,nPeaks+1);
		func->SetLineWidth(1);
		func->SetNpx(1000);
		for (i=0;i<nPeaks;i++)
		{
			func->FixParameter(i,fit_par[2*i]);
		}
		func->FixParameter(nPeaks,(double)j);
		func->SetLineColor(j+3);
		func->DrawCopy("CSAME");
		delete func;
	}

	c1->SaveAs(name);

	delete c1;
}

void LinFitter::makeGuess()
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
			guess_t[j] = startTime+sampleInterval*(i-1)-shape->getPeak();
			j++;
		}
	}
	free(dy);
}

void LinFitter::print_guess()
{
	int i;
	printf("Times:\t\t");
	for (i=0;i<nPeaks;i++)
		printf("%lf\t",guess_t[i]);
	printf("\n");
}

void LinFitter::setVerbosity(int verbosity)
{
	gMinuit->SetPrintLevel(verbosity);
	verbose = verbosity;
}

double LinFitter::Evaluate(double *x, double *p)
{
	//return sqrt(doLinFit(x, NULL));
	return doLinFit(x, NULL);
}

void LinFitter::plotFCN(Event *evt, const char *name, int res, const char *style, double x_min, double x_max, double y_min, double y_max)
{
	if (nPeaks<2) return;

	TGraph *graph;

	TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,750);
	c1->SetLogz();

	TF2 *func = new TF2("FCN",this,&LinFitter::Evaluate,x_min,x_max,y_min,y_max,0,"LinFitter","LinFit");
	func->SetNpy(res);
	func->SetNpx(res);
	func->SetLineWidth(1);
	//func->SetContour(20,z);

	//func->DrawCopy("surf4"); //makes a fancy gouraud-shaded surface plot
	gStyle->SetPalette(1);
	//func->DrawCopy("cont1z");
	//func->DrawCopy("colz");
	func->DrawCopy(style);
	delete func;

	double *t = (double *) calloc(nPeaks, sizeof(double));
	double *t2 = (double *) calloc(nPeaks, sizeof(double));
	evt->getTimes(t);
	t2[0] = t[1];
	graph = new TGraph(1,t,t2);
	graph->SetMarkerColor(1);
	graph->DrawClone("*");
	delete graph;

	graph = new TGraph(1,t2,t);
	graph->SetMarkerColor(1);
	graph->DrawClone("*");
	delete graph;

	t2[0] = fit_t[1];
	graph = new TGraph(1,fit_t,t2);
	graph->SetMarkerColor(2);
	graph->DrawClone("*");
	delete graph;

	c1->SaveAs(name);
	delete c1;
	free(t);
	free(t2);
}

void LinFitter::doGridFit()
{
	double gridStep = 8.0;
	//double scanStep = 2.0;
	double *min, *max;
	min = (double *) calloc(nPeaks, sizeof(double));
	max = (double *) calloc(nPeaks, sizeof(double));
	int i,j;
	double chisq, best_chisq, last_chisq;
	for (i=0;i<nPeaks;i++)
	{
		min[i] = -96.0;
		max[i] = 120.0;
	}

	best_chisq=-1;
	for (j=0;j<7;j++)
	{
		for (i=0;i<nPeaks;i++)
		{
			guess_t[i] = max[i];
		}
		last_chisq = best_chisq;
		best_chisq = -1;
		if (verbose>2)
		{
			printf("Grid fit, step size %lf\n",gridStep);
			for (i=0;i<nPeaks;i++)
			{
				printf("Peak %d: min %lf, max %lf\n",i,min[i],max[i]);
			}
		}
		while (guess_t[nPeaks-1]>min[nPeaks-1])
		{
			guess_t[0] -= gridStep;
			for (i=0;i<nPeaks-1 && guess_t[i]<min[i];i++)
			{
				guess_t[i+1] -= gridStep;
				guess_t[i] = guess_t[i+1]<max[i] ? guess_t[i+1] : max[i];
			}
			chisq = doLinFit(guess_t);
			if (best_chisq<0 || chisq<best_chisq)
			{
				best_chisq = chisq;
				for (i=0;i<nPeaks;i++)
				{
					fit_t[i] = guess_t[i];
				}
			}
		}
		if (verbose>2)
		{
			printf("Step %d: Best chisq %lf, improvement %f, times:",j,best_chisq,last_chisq-best_chisq);
			for (i=0;i<nPeaks;i++)
			{
				printf("\t%lf",fit_t[i]);
			}
			printf("\n");
		}
		if (last_chisq!=-1 && last_chisq!=best_chisq && last_chisq-best_chisq<0.01) break;
		/*
		if (nPeaks>1)
		{
			for (i=0;i<nPeaks;i++)
			{
				guess_t[i] = fit_t[i];
			}
			max[1] = fit_t[1]+gridStep;
			guess_t[1] = fit_t[1]-gridStep;
			while (guess_t[1]<max[1])
			{
				guess_t[1]+=gridStep/8;
				chisq = doLinFit(guess_t);
				//printf("guess chisq: %f\n",chisq);
				if (best_chisq<0 || chisq<best_chisq)
				{
					best_chisq = chisq;
					fit_t[1] = guess_t[1];
				}
			}
			guess_t[1] = fit_t[1];
			guess_t[0] = min[0];
			while (guess_t[0]<guess_t[1])
			{
				guess_t[0]+=scanStep;
				chisq = doLinFit(guess_t);
				//printf("guess chisq: %f\n",chisq);
				if (best_chisq<0 || chisq<best_chisq)
				{
					best_chisq = chisq;
					fit_t[0] = guess_t[0];
				}
			}
		}
		*/
		for (i=0;i<nPeaks;i++)
		{
			min[i] = fit_t[i]-2*gridStep;
			max[i] = fit_t[i]+2*gridStep;
		}
		gridStep /= 4.0;
	}
	for (i=0;i<nPeaks;i++)
	{
		guess_t[i] = fit_t[i];
	}
	doLinFit(fit_t,fit_par);
	status = 0;
	dof = getNumUsedSamples()-2*nPeaks;
	free(min);
	free(max);
}

void LinFitter::doFit()
{
	doGridFit();
}

void LinFCN(Int_t&npar, Double_t*gin, Double_t&f, Double_t*par, Int_t flag)
{
	LinFitter *theFitter = (LinFitter *)currentFitter;
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
			f = theFitter->doLinFit(par,NULL);
			break;
	}
}

Double_t LinFCNCurve(Double_t *x, Double_t *par)
{
	LinFitter *theFitter = (LinFitter *)currentFitter;
	int i, n;
	n = theFitter->getNumPeaks();
	double *newpar = (double *) calloc(n, sizeof(double));
	int k = (int)par[n];
	for (i=0;i<n;i++)
		newpar[i] = par[i];
	newpar[k] = x[0];
	double y = theFitter->doLinFit(newpar,NULL);
	free(newpar);
	return y;
}

