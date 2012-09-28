#include "apvtime.hh"
#include "Event.hh"
#include "ShapingCurve.hh"
#include "Samples.hh"
#include "Fitter.hh"
#include "MinuitFitter.hh"
#include "LinFitter.hh"
#include "AnalyticFitter.hh"

#include <stdio.h>
#include <stdarg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>

#include "TFile.h"
#include "TTree.h"
#include "SmoothShapingCurve.hh"
#include "TMath.h"

int main(int argc,char** argv) //fitter type (1: Minuit, 2: linear, 3: analytic); nEvents; nPeaks (1 or 2)
{
	int c;

	int fitterType = 3;
	int nPeaks = 1;
	int nEvents = 100;

	while ((c = getopt(argc,argv,"hf:n:p:")) !=-1)
		switch (c)
		{
			case 'h':
				printf("-h: print this help\n");
				printf("-f: fitter type\n");
				printf("-n: number of events\n");
				printf("-p: number of peaks per event\n");
				return(0);
				break;
			case 'f':
				fitterType = atoi(optarg);
				break;
			case 'n':
				nEvents = atoi(optarg);
				break;
			case 'p':
				nPeaks = atoi(optarg);
				break;
			case '?':
				printf("Invalid option or missing option argument; -h to list options\n");
				return(1);
			default:
				abort();
		}

	printf("type %d, nEvents %d, nPeaks %d\n", fitterType, nEvents, nPeaks);

	//initialize the RNG
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();

	T = gsl_rng_mt19937;
	r = gsl_rng_alloc (T);


	ShapingCurve *myShape = new ShapingCurve(35.0);

	Fitter *myFitter, *myFitter2;
	switch (fitterType)
	{
		case 1:
			myFitter = new MinuitFitter(myShape,6,1,1.0);
			myFitter2 = new MinuitFitter(myShape,6,2,1.0);
			break;
		case 2:
			myFitter = new LinFitter(myShape,6,1,1.0);
			myFitter2 = new LinFitter(myShape,6,2,1.0);
			break;
		case 3:
			myFitter = new AnalyticFitter(myShape,6,1,1.0);
			myFitter2 = new AnalyticFitter(myShape,6,2,1.0);
			break;
	}

	myFitter->setVerbosity(-1);
	myFitter2->setVerbosity(-1);
	Event *myEvent = new Event(myShape,r,1.0);

	Samples *mySamples = new Samples(6,24.0);

	TFile *hfile = new TFile("stuff.root","RECREATE","Stuff");
	TTree *tree = new TTree("myTree","A ROOT tree");


	int event;
	double time, height;
	int status;
	double fit_times[2], fit_heights[2];
	double times[2], heights[2];
	double *temp;
	int peaks;
	int badfits=0;
	char name[100];
	int start,end;

	int status2;
	double fit_times2[2], fit_heights2[2];

	tree->Branch("true_times",times,"true_times[2]/D");
	tree->Branch("true_heights",heights,"true_heights[2]/D");

	tree->Branch("status",&status,"status/I");
	tree->Branch("fit_times",fit_times,"fit_times[2]/D");
	tree->Branch("fit_heights",fit_heights,"fit_heights[2]/D");

	tree->Branch("status2",&status2,"status2/I");
	tree->Branch("fit_times2",fit_times2,"fit_times2[2]/D");
	tree->Branch("fit_heights2",fit_heights2,"fit_heights2[2]/D");

	double true_par[4], fit_par[4], fit_err[4];
	double true_chisq, fit_chisq;
	double prob;
	int dof;
	double fit_par2[4], fit_err2[4];
	double fit_chisq2;
	double prob2;
	int dof2;
	int samples;
	tree->Branch("samples",&samples,"samples/I");
	tree->Branch("true_chisq",&true_chisq,"true_chisq/D");
	//tree->Branch("true_par",true_par,"true_par[4]/D");

	tree->Branch("dof",&dof,"dof/I");
	//tree->Branch("fit_par",fit_par,"fit_par[4]/D");
	tree->Branch("fit_err",fit_err,"fit_err[4]/D");
	tree->Branch("fit_chisq",&fit_chisq,"fit_chisq/D");
	tree->Branch("prob",&prob,"prob/D");

	tree->Branch("dof2",&dof2,"dof2/I");
	//tree->Branch("fit_par2",fit_par2,"fit_par2[4]/D");
	tree->Branch("fit_err2",fit_err2,"fit_err2[4]/D");
	tree->Branch("fit_chisq2",&fit_chisq2,"fit_chisq2/D");
	tree->Branch("prob2",&prob2,"prob2/D");

	true_par[2] = 0;
	true_par[3] = 0;

	for (event=0;event<nEvents;event++)
	{
		//printf("event %d\n",event);
		gsl_rng_set (r,event); //seed = event


		time = gsl_ran_flat(r, 24.0, 48.0);
		//time = gsl_ran_flat(r, -0.0, 90.0);
		height = 25.0;
		myEvent->addHit(time,height);
		true_par[0] = time;
		true_par[1] = height;


		if (nPeaks>1)
		{
			time = gsl_ran_flat(r, -100.0, 120.0);
			//time = gsl_ran_flat(r, -100.0, 0.0);
			height = 25.0;
			myEvent->addHit(time,height);
			true_par[2] = time;
			true_par[3] = height;
		}

		//myEvent->sortHits();
		myEvent->getTimes(times);
		myEvent->getHeights(heights);
		mySamples->readEvent(myEvent, 0.0);


		myFitter->readSamples(mySamples);
		myFitter2->readSamples(mySamples);



		myFitter->doFit();
		myFitter2->doFit();

		status = myFitter->getStatus();
		myFitter->getFitPar(fit_par);
		myFitter->getFitErr(fit_err);
		dof = myFitter->getDOF();
		samples = myFitter->getNumUsedSamples();

		myFitter->getFitTimes(fit_times);
		myFitter->getFitHeights(fit_heights);

		fit_chisq = myFitter->getChisq(fit_par);
		true_chisq = myFitter->getChisq(true_par);
		prob = TMath::Prob(fit_chisq,dof);

		status2 = myFitter2->getStatus();
		myFitter2->getFitPar(fit_par2);
		myFitter2->getFitErr(fit_err2);
		dof2 = myFitter2->getDOF();

		myFitter2->getFitTimes(fit_times2);
		myFitter2->getFitHeights(fit_heights2);

		fit_chisq2 = myFitter2->getChisq(fit_par2);
		prob2 = TMath::Prob(fit_chisq2,dof2);

		tree->Fill();

		/*
		   printf("Fit %d (status %d, chisq %lf)\n",event,status,chisq);
		   printf("True params (chisq %lf), samples:\n",myFitter->memberFCN(par));
		   myEvent->print();
		   mySamples->print();
		   printf("Guessed params (chisq %lf):\n",myFitter->guessChisq());
		   myFitter->print_guess();
		   printf("Fitted params:\n");
		   myFitter->print_fit();
		   */

		if (status == 0 && fit_chisq>true_chisq)
			status = -1;

		if (status2 == 0 && fit_chisq2>true_chisq)
			status2 = -1;

		//printf("Fit %d: chisq1 %f, chisq2 %f, true chisq %f\n",i,fit_chisq, fit_chisq2, true_chisq);
		//myEvent->print();

		if (status != 0 && status2 != 0)
		{
			badfits++;
			printf("Rerunning bad fit %d (status %d, chisq %lf)\n",event,status,fit_chisq);
			printf("(status2 %d, chisq2 %lf)\n",status2,fit_chisq2);
			printf("True params (chisq %lf), samples:\n",true_chisq);
			myEvent->print();
			mySamples->print();
			printf("Fitted params:\n");
			myFitter->print_fit();
			myFitter2->print_fit();

			if (fitterType==1)
			{
				//MinuitFitter-specific
				printf("Guessed params (chisq %lf):\n",myFitter->guessChisq());
				myFitter->print_guess();
				printf("Guessed params (chisq %lf):\n",myFitter2->guessChisq());
				myFitter2->print_guess();
			}

			if (fitterType==2)
			{
				//LinFitter-specific
				sprintf(name,"fit_%06d_fcn.png",event);
				myFitter2->plotFCN(myEvent,name,300,"cont1z",-100,120,-100,120);
				sprintf(name,"fit_%06d_fcn_zoom.png",event);
				myFitter2->plotFCN(myEvent,name,100,"cont1z",fit_par2[0]-5.0,fit_par2[0]+5.0,fit_par2[2]-5.0,fit_par2[2]+5.0);
				sprintf(name,"fit_%06d_fcn_zoom_true.png",event);
				myFitter2->plotFCN(myEvent,name,100,"cont1z",true_par[0]-5.0,true_par[0]+5.0,true_par[2]-5.0,true_par[2]+5.0);
			}

			/*
			   sprintf(name,"fit_%06d.png",event);
			   myFitter->plotFit(myEvent,name);
			   sprintf(name,"fit_%06d_2.png",event);
			   myFitter2->plotFit(myEvent,name);

			   myFitter->setVerbosity(3);
			   myFitter2->setVerbosity(3);
			   myFitter->doFit();
			   myFitter2->doFit();
			   myFitter->setVerbosity(0);
			   myFitter2->setVerbosity(0);
			   */
		}

		/*
		   if (status == 0)
		   {
		   sprintf(name,"goodfit_%6d.png",i);
		   myFitter->plotFit(myEvent,name);
		   }
		   */

		myEvent->clear();
	}
	printf("Bad fits: %d\n",badfits);
	hfile->Write();

	//clean up
	gsl_rng_free (r);
	return 0;
}

