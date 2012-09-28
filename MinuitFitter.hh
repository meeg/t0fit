#ifndef MINUIT_FITTER_HH
#define MINUIT_FITTER_HH
#include "Fitter.hh"

#include "Event.hh"
#include "Samples.hh"
#include "ShapingCurve.hh"
#include "TObject.h"
#include "TMinuit.h"

class MinuitFitter : public Fitter
{
	public:
		MinuitFitter(ShapingCurve *sc, int sCount, int pCount, double noise);
		~MinuitFitter();
	private:
		bool tweakFCN;
		double *guess_par;
	public:
		void doMinuitFit();
		void doFit();
		double memberFCN(double *par);
		void memberGradFCN(double *par, double *gin);
		void initMinuit();
		void setTweakFCN(bool x);
		bool getTweakFCN();
		void plotFit(Event *evt,char *name);
		void setVerbosity(int verbosity);
		void makeGuess();
		void setGuess(int i, double time, double height);
		void print_guess();
		double guessChisq();
};

void FCN(Int_t&npar, Double_t*gin, Double_t&f, Double_t*par, Int_t flag);
Double_t FCNCurve(Double_t *x, Double_t *par); //with tweak
Double_t FCNCurve2(Double_t *x, Double_t *par); //no tweak

#endif
