#ifndef FITTER_HH
#define FITTER_HH

#include "Event.hh"
#include "Samples.hh"
#include "ShapingCurve.hh"
#include "TObject.h"
#include "TMinuit.h"

class Fitter : public TObject
{
	public:
		Fitter(ShapingCurve *sc, int sCount, int pCount, double noise);
		~Fitter();
	protected:
		int nPeaks; //number of peaks to fit
		int nSamples; //number of samples
		double * y; //samples
		double sigma_noise; //noise level
		double * sigma; //sample error
		double startTime, sampleInterval; //sample start time, interval
		int status; //0 for good fit
		ShapingCurve *shape;
		double *fit_par, *fit_err;
		bool * useSample;
		int verbose;
		int dof;

	public:
		void print_fit();
		void setVerbosity(int verbosity);
		void getFitTimes(double *times);
		void getFitHeights(double *heights);
		void getFitPar(double *par);
		void getFitErr(double *err);
		void readSamples(Samples *input);
		double getChisq(double *par, bool use_all=false);
		double getFitChisq();
		int getDOF();
		int getStatus();
		ShapingCurve * getShape();
		int getNumPeaks();
		void printResiduals();
		double getSignal(double time);
		void setSigmaNoise(double sigma);
		double getFirstPeak();
		//void setUsedSamples();

		void sortFit();
		void setUseSample(int i, bool use);
		int getNumUsedSamples(int start=0, int end=-1);

		virtual void doFit() {}
		virtual void plotFit(Event *evt, const char *name) =0;
		virtual void plotFCN(Event *evt, const char *name, int res, const char *style, double x_min, double x_max, double y_min, double y_max) {}
		virtual void print_guess() {}
		virtual double guessChisq() { return -1;}
};


extern TMinuit *gMinuit;
extern Fitter *currentFitter;

Double_t evtCurve(Double_t *x, Double_t *par);
#endif

