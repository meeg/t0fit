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
#include <algorithm>
#include <cmath>
using namespace std;


LinFitter::LinFitter(ShapingCurve *sc, int sCount, int pCount, double noise) : Fitter(sc, sCount, pCount, noise)
{
	if (gMinuit==NULL) gMinuit = new TMinuit();
	guess_t = (double *) calloc(pCount, sizeof(double));
	coeff_mat = gsl_matrix_alloc(pCount,pCount);
	a_vec = gsl_vector_alloc(pCount);
	fit_t = (double *) calloc(pCount, sizeof(double));
	fit_t_err = (double *) calloc(pCount, sizeof(double));
	ework = gsl_eigen_symmv_alloc(nPeaks);
	hess_mat = gsl_matrix_alloc(nPeaks,nPeaks);
	evec = gsl_matrix_alloc(nPeaks,nPeaks);
	proj_mat = gsl_matrix_alloc(nPeaks,nPeaks);
	eval = gsl_vector_alloc(nPeaks);
	grad = gsl_vector_alloc(nPeaks);
	eig_mat = gsl_matrix_alloc(nPeaks,nPeaks);
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
	gsl_eigen_symmv_free(ework);
	gsl_matrix_free(hess_mat);
	gsl_matrix_free(evec);
	gsl_matrix_free(proj_mat);
	gsl_vector_free(eval);
	gsl_vector_free(grad);
	gsl_matrix_free(eig_mat);
}

void LinFitter::doFit()
{
	//doGridFit();
	//doMinuitFit();
	doRecursiveFit();
	//doScanFit();
}

double LinFitter::doLinFit(double *times, double *par)
{
	int i,j;
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
	//gsl_vector_set(a_vec,gsl_vector_min_index(a_vec),0.0);
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
	gMinuit->SetPrintLevel(max(-1,verbose-3));
	gMinuit->Command("CLE");

	gMinuit->Command("SET NOW");
	//gMinuit->Command("SET WAR");
	gMinuit->Command("SET STR 2");
	gMinuit->Command("SET ERR 1");
}

void LinFitter::doMinuitFit()
{
	//makeGuess();
	for (int j=0;j<nPeaks;j++)
	{
		guess_t[j] = sampleInterval*(j);
	}
	if (nPeaks==1) {
		int numUsedSamples = 0;
		int numPositiveSamples = 0;
		int numBigSamples = 0;
		int firstUsedSample = nSamples;
		int lastUsedSample;
		int firstBigSample = nSamples;
		for (int i=0;i<nSamples;i++) if (useSample[i]) {
			numUsedSamples++;
			lastUsedSample = i;
			if (i<firstUsedSample) firstUsedSample = i;
			if (y[i]>0) {
				numPositiveSamples++;
				if (y[i]>3.0*sigma_noise) {
					numBigSamples++;
					if (i<firstBigSample) firstBigSample = i;
				}
			}
		}
		bool made_guess = false;
		bool made_bestfit = false;
		if (numUsedSamples==1) {
			if (useSample[0]) {
				fit_t[0] = -500.0;
				made_bestfit = true;
			} else {
				fit_t[0] = sampleInterval*(firstUsedSample-0.1);
				made_bestfit = true;
			}
		}
		else if (numPositiveSamples==1 && y[lastUsedSample]>0) {
			fit_t[0] = sampleInterval*(lastUsedSample-0.1);
			made_bestfit = true;
		}
		else if (numBigSamples==1 && y[lastUsedSample]>3.0*sigma_noise && useSample[lastUsedSample-1] && y[lastUsedSample-1]<0) {
			fit_t[0] = sampleInterval*(lastUsedSample-0.1);
			made_bestfit = true;
		}
		else if (numUsedSamples==2) {
			guess_t[0] = sampleInterval*(firstUsedSample-0.1);
			made_guess = true;
		}
		/*
		if (!made_bestfit && !made_guess && firstBigSample<nSamples) {
			if (firstBigSample+1<nSamples && useSample[firstBigSample+1]) {
				guess_t[0] = sampleInterval*(firstBigSample-0.1);
			} else {
				guess_t[0] = sampleInterval*(firstBigSample-1.1);
			}
			made_guess = true;
		}*/
		/*
		   else for (int i=0;i<nSamples;i++) if (useSample[i] && y[i]>3.0*sigma_noise) {
		   if (i>0 && useSample[i-1]) {
		   if (y[i-1]>0) {
		   guess_t[0] = sampleInterval*(i-1.1);
		   } else {
		   }
		   }
		   else
		   guess_t[0] = sampleInterval*(i-0.1);
		   made_guess = true;
		   break;
		   }*/
		if (made_bestfit) {
			if (verbose>1)
				printf("fit is degenerate; chose a best-fit value: %f\n",fit_t[0]);
			doLinFit(fit_t,fit_par);
			return;
		}
	}
	if (verbose>1) {
		printf("guesses:\t");
		for (int i=0; i<nPeaks; i++) printf("%f\t",guess_t[i]);
		printf("\n");
	}
	//doGridFit(2);
	initMinuit();
	int i;
	for (i=0; i<nPeaks; i++)
	{
		gMinuit->DefineParameter(i,"time",guess_t[i],25.0,-500.0,120.0);
	}

	arglist[0] = 1000;
	//gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);
	gMinuit->mnexcm("SIMPLEX", arglist ,1,ierflg);
	//gMinuit->mnexcm("MINIMIZE", arglist ,1,ierflg);

	status = ierflg;
	/*
	   if (nPeaks==2) {
	   double new_fval;
	   double *init_par = new double[nPeaks];
	   double *best_par = new double[nPeaks];
	   getPar(init_par);
	   double best_fval = doLinFit(init_par);
	   if (verbose > 0)
	   printf("Before valley: chisq %f\n",best_fval);
	   memcpy(best_par,init_par,nPeaks*sizeof(double));

	   double *move = new double[nPeaks];
	   move[0] = 1.0;
	   for (int i=1;i<nPeaks;i++) move[i] = 0.0;
	   new_fval = valleyImprove(move);
	   if (verbose > 0)
	   printf("Valley right: new chisq %f\n",new_fval);
	   if (new_fval < best_fval) {
	   best_fval = new_fval;
	   getPar(best_par);
	   }

	   setPar(init_par);
	   move[0] = -1.0;
	   for (int i=1;i<nPeaks;i++) move[i] = 0.0;
	   new_fval = valleyImprove(move);
	   if (verbose > 0)
	   printf("Valley left: new chisq %f\n",new_fval);
	   if (new_fval < best_fval) {
	   best_fval = new_fval;
	   getPar(best_par);
	   }
	   setPar(best_par);
	   arglist[0] = 1000;
	   gMinuit->mnexcm("SIMPLEX", arglist ,1,ierflg);
	   status = ierflg;
	   }
	   */
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

void LinFitter::doRecursiveFit()
{
	if (nPeaks==1) {
		//doGridFit(2);
		//lineScan(0,fit_t[0],5.0);
		doMinuitFit();
		return;
	}

	double *best_times = NULL;
	//double *times = new double[nPeaks];
	double best_fval, fval;
	int best_split;

	LinFitter *fitter1 = new LinFitter(shape,nSamples,1,sigma_noise);
	LinFitter *fitter2 = new LinFitter(shape,nSamples,nPeaks-1,sigma_noise);
	fitter1->setVerbosity(verbose);
	fitter2->setVerbosity(verbose);
	Samples *mySamples = new Samples(6,24.0);
	mySamples->readEvent(y,startTime);
	fitter1->readSamples(mySamples);
	for (int split=1;split<nSamples;split++) {
		for (int i=0;i<nSamples;i++) {
			fitter1->setUseSample(i,i<split);
		}
		fitter1->doRecursiveFit();
		fitter1->getFitTimes(fit_t);
		double *y_sub = new double[nSamples];
		double total = 0;
		//int firstBigSample=nSamples-1;
		for (int i=0;i<nSamples;i++) {
			y_sub[i] = y[i]-fitter1->getSignal(startTime+i*sampleInterval);
			total+=y_sub[i];
			//if (i<firstBigSample && y_sub[i]>sigma_noise) firstBigSample = i;
		}
		if (total>0) {
			mySamples->readEvent(y_sub,startTime);
			if (verbose>1) {
				printf("samples with first peak subtracted:\t");
				mySamples->print();
			}
			fitter2->readSamples(mySamples);
			for (int i=0;i<nSamples;i++) {
				fitter2->setUseSample(i,true);
				//fitter2->setUseSample(i,split-1<=i||firstBigSample<=i);
			}
			fitter2->doRecursiveFit();
			fitter2->getFitTimes(fit_t+1);
		}
		delete[] y_sub;

		if (verbose>1) {
			fval = doLinFit(fit_t);
			printf("split %d, before refit: chisq = %f\n",split,fval);
			for (int i=0; i<nPeaks; i++) printf("time %d: %f\t",i,fit_t[i]);
			printf("\n");
		}

		initMinuit();
		for (int i=0; i<nPeaks; i++)
		{
			gMinuit->DefineParameter(i,"time",fit_t[i],25.0,-500.0,120.0);
		}
		arglist[0] = 1000;
		gMinuit->mnexcm("SIMPLEX", arglist ,1,ierflg);
		status = ierflg;
		if (status == 0)
		{
			arglist[0] = 1000;
			gMinuit->mnexcm("IMPROVE", arglist ,1,ierflg);
		}
		getPar(fit_t);

		fval = doLinFit(fit_t);
		if (verbose>1) {
			printf("split %d, after refit: chisq = %f\n",split,fval);
			for (int i=0; i<nPeaks; i++) printf("time %d: %f\t",i,fit_t[i]);
			printf("\n");
		}

		if (best_times == NULL) {
			best_fval = fval;
			best_split = split;
			best_times = new double[nPeaks];
			memcpy(best_times,fit_t,nPeaks*sizeof(double));
		} else if (fval<best_fval) {
			best_fval = fval;
			best_split = split;
			memcpy(best_times,fit_t,nPeaks*sizeof(double));
		}
	}

	delete fitter1;
	delete fitter2;
	delete mySamples;

	if (verbose>1) {
		printf("Best split: %d, chisq = %f\n",best_split,best_fval);
		printf("Times:\t\t");
		for (int i=0;i<nPeaks;i++) printf("%lf\t",best_times[i]);
		printf("\n");
	}

	memcpy(fit_t,best_times,nPeaks*sizeof(double));
	delete[] best_times;

	/*
	initMinuit();
	for (int i=0; i<nPeaks; i++)
	{
		gMinuit->DefineParameter(i,"time",fit_t[i],25.0,-500.0,120.0);
	}
	arglist[0] = 1000;
	gMinuit->mnexcm("SIMPLEX", arglist ,1,ierflg);
	status = ierflg;
	if (status == 0)
	{
		arglist[0] = 1000;
		gMinuit->mnexcm("IMPROVE", arglist ,1,ierflg);
	}
	getPar(fit_t);
	*/

	fval = doLinFit(fit_t,fit_par);

	/*
	if (verbose>1) {
		printf("After refit: chisq = %f\n",fval);
		printf("Times:\t\t");
		for (int i=0;i<nPeaks;i++) printf("%lf\t",fit_t[i]);
		printf("\n");
	}
	*/
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

	func = new TF1("curve",evt,&Event::Evaluate,-100.0,200.0,0);
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
	printf("Times:\t\t");
	for (int i=0;i<nPeaks;i++)
		printf("%lf\t",guess_t[i]);
	printf("\n");
}

void LinFitter::setVerbosity(int verbosity)
{
	verbose = verbosity;
}

double LinFitter::Evaluate(double *x, double *p)
{
	//return sqrt(doLinFit(x));
	return doLinFit(x);
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

void LinFitter::lineScan(int i, double start, double step) {
	fit_t[i] = start;
	double fval, old_fval;
	bool success, old_success;
	int k=0;

	fval = doLinFit(fit_t);
	success = false;

	while (k<1000) {
		old_fval = fval;
		fit_t[i]+=step;
		fval = doLinFit(fit_t);
		old_success = success;
		success = fval<=old_fval;
		if (success) {
			step *= 3.0;
		} else {
			fit_t[i]-=step;
			if (old_success && fval-old_fval < 0.01 && fval!=old_fval) return;
			fval = doLinFit(fit_t);
			step *= -0.4;
			//if (old_success) return;
		}
		k++;
	}
}

void LinFitter::doScanFit()
{
	double *best_t = NULL;
	double best_fval, fval;

	initMinuit();
	for (int i=0; i<nPeaks; i++)
	{
		gMinuit->DefineParameter(i,"time",fit_t[i],5.0,0,0);
	}

	if (nPeaks==1) {
		lineScan(0, 0.0, 5.0);
		setPar(fit_t);
	}

	if (nPeaks==2) {
		double *move = new double[nPeaks];

		fit_t[0] = 110.0;
		lineScan(1, 0.0, 5.0);

		if (verbose>0) {
			printf("Times:\t\t");
			for (int i=0;i<nPeaks;i++)
				printf("%lf\t",fit_t[i]);
			printf("\n");
		}

		/*
		   setPar(fit_t);
		   arglist[0] = 1000;
		   gMinuit->mnexcm("SIMPLEX", arglist ,1,ierflg);
		   getPar(fit_t);
		   fval = doLinFit(fit_t);
		   */

		move[0] = -1.0;
		for (int i=1;i<nPeaks;i++) move[i] = 0.0;
		setPar(fit_t);
		fval = valleyImprove(move);
		getPar(fit_t);

		if (best_t == NULL) {
			best_fval = fval;
			best_t = new double[nPeaks];
			memcpy(best_t,fit_t,nPeaks*sizeof(double));
		} else if (fval<best_fval) {
			best_fval = fval;
			memcpy(best_t,fit_t,nPeaks*sizeof(double));
		}
		fit_t[0] = -200.0;
		lineScan(1, 0.0, 5.0);

		if (verbose>0) {
			printf("Times:\t\t");
			for (int i=0;i<nPeaks;i++)
				printf("%lf\t",fit_t[i]);
			printf("\n");
		}

		/*
		   setPar(fit_t);
		   arglist[0] = 1000;
		   gMinuit->mnexcm("SIMPLEX", arglist ,1,ierflg);
		   getPar(fit_t);
		   fval = doLinFit(fit_t);
		   */

		move[0] = 1.0;
		for (int i=1;i<nPeaks;i++) move[i] = 0.0;
		setPar(fit_t);
		fval = valleyImprove(move);
		getPar(fit_t);

		if (best_t == NULL) {
			best_fval = fval;
			best_t = new double[nPeaks];
			memcpy(best_t,fit_t,nPeaks*sizeof(double));
		} else if (fval<best_fval) {
			best_fval = fval;
			memcpy(best_t,fit_t,nPeaks*sizeof(double));
		}

		setPar(best_t);
	}

	arglist[0] = 1000;
	gMinuit->mnexcm("SIMPLEX", arglist ,1,ierflg);
	getPar(fit_t);
	doLinFit(fit_t,fit_par);
}

void LinFitter::doGridFit(int nIterations)
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
	for (j=0;j<nIterations;j++)
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
			//printf("chisq %f: time %f\n",chisq,guess_t[0]);
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

double LinFitter::valleyImprove(double *move)
{
	double *best_par = new double[nPeaks];
	double best_fval;

	getPar(best_par);

	double init_fval = doLinFit(best_par);
	best_fval = init_fval;

	double old_fval = init_fval;
	double step = 1.0;
	int k = 0;
	while(true) {
		double new_fval = valleyStep(move, step, 0.0);
		if (new_fval<best_fval) {
			//printf("New best: %f\n",new_fval);
			best_fval = new_fval;
			getPar(best_par);
			//printf("Best par:\t");
			//for (int i=0;i<nPeaks;i++)
			//	printf("%f\t",best_par[i]);
			//printf("\n");
		}

		//if (new_fval > old_fval + 0.1 || new_fval>init_fval + 10) {
		if (new_fval>init_fval + 10) {
			if (verbose > 0)
				printf("%d valley steps\n",k);
			break;
		}
		old_fval = new_fval;
		if (++k>1000) {
			if (verbose > 0)
				printf("%d valley steps\n",k);
			break;
		}
	}

	/*
	   setPar(best_par);
	   arglist[0] = 1000;
	   gMinuit->mnexcm("SIMPLEX", arglist ,1,ierflg);
	   getPar(best_par);
	   */

	setPar(best_par);
	return doLinFit(best_par);
	}

	double LinFitter::valleyStep(double *dir, double &step, double mindot)
	{
		arglist[0] = 1000;
		gMinuit->mnexcm("HESSE", arglist ,1,ierflg);
		int maxPar = gMinuit->fMaxpar;

		for (int i = 0;i<nPeaks;i++) {
			for (int j = 0;j<nPeaks;j++) {
				gsl_matrix_set(hess_mat,i,j,gMinuit->fP[i+j*maxPar]);
				if (verbose>1) printf("%f\t",gMinuit->fP[i+j*maxPar]);
			}
			gsl_vector_set(grad,i,gMinuit->fGrd[i]);
			if (verbose>1) printf("\t%f\n",gMinuit->fGrd[i]);
		}

		//gsl_matrix_set_identity(proj_mat);
		//gsl_matrix_set_zero(eig_mat);
		//gsl_blas_dsyr(CblasLower,1.0/gsl_blas_dnrm2(grad),grad,proj_mat);

		//gsl_blas_dsymm(CblasLeft,CblasLower,1.0,proj_mat,hess_mat,0.0,eig_mat);

		//gsl_eigen_symmv(eig_mat,eval,evec,ework);
		gsl_eigen_symmv(hess_mat,eval,evec,ework);
		gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_ASC);
		/*
		   printf("Eigenvalues:\n");
		   gsl_vector_fprintf(stdout,eval,"%f");
		   printf("Eigenvectors:\n");
		   gsl_matrix_fprintf(stdout,evec,"%f");
		   */
		gsl_vector_view evec0 = gsl_matrix_column(evec, 0);
		gsl_vector_view evec1 = gsl_matrix_column(evec, 1);

		double *move = evec0.vector.data;
		double *correct = evec1.vector.data;

		double dotproduct = 0;
		for (int i=0;i<nPeaks;i++) dotproduct+=move[i]*dir[i];

		/*
		   if (abs(dotproduct) < mindot) {
		   return valleyStep(dir,step/2,mindot);
		   }
		   */

		if (verbose > 1) {
			printf("dot: %f\n",dotproduct);
		}

		if (dotproduct<0) {
			for (int i=0;i<nPeaks;i++) move[i]*=-1;
		}

		/*
		   if (abs(dotproduct) > 0.95) {
		   step = min(2*step,5.0);
		   } else {
		   step = max(0.5*step,1.0);
		   }
		   */

		memcpy(dir,move,nPeaks*sizeof(double));

		double *par = new double[nPeaks];
		getPar(par);

		if (verbose > 1) {
			printf("old par:\t");
			for (int i=0;i<nPeaks;i++)
				printf("%f\t",par[i]);
			printf("\n");
			printf("step:\t");
			for (int i=0;i<nPeaks;i++)
				printf("%f\t",step*move[i]);
			printf("\n");
		}

		for (int i=0;i<nPeaks;i++) {
			//	if (par[i]<-150.0) return 100000;
			par[i]+=step*move[i];
		}
		setPar(par);

		if (verbose > 1) {
			printf("stepped par:\t");
			for (int i=0;i<nPeaks;i++)
				printf("%f\t",par[i]);
			printf("\n");
		}

		gMinuit->mnhes1();
		double slope = 0;
		for (int i=0;i<nPeaks;i++) slope+=correct[i]*gMinuit->fGrd[i];
		if (slope > 0) {
			slope*=-1;
			for (int i=0;i<nPeaks;i++) correct[i]*=-1;
		}
		double fval = doLinFit(par);
		gMinuit->mnline(par,fval,correct,slope,0.05);

		if (verbose > 1) {
			getPar(par);
			printf("mnline par:\t");
			for (int i=0;i<nPeaks;i++)
				printf("%f\t",par[i]);
			printf("\n");
		}

		return fval;
	}

	void LinFitter::getPar(double *par) {
		double dummy;
		for (int i=0; i<nPeaks; i++) gMinuit->GetParameter(i,par[i],dummy);
	}

	void LinFitter::setPar(double *par) {
		for (int i=0; i<nPeaks; i++) {
			arglist[0] = i+1;
			arglist[1] = par[i];
			gMinuit->mnexcm("SET PARAM",arglist,2,ierflg);
		}
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
				/*
				   for (int i=0;i<npar;i++) if (par[i]<-200.0) {
				   f = 10000000000.0;
				   return;
				   }
				   */
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
