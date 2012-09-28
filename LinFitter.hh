#ifndef LIN_FITTER_HH
#define LIN_FITTER_HH
#include "Fitter.hh"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

class LinFitter : public Fitter
{
	public:
		LinFitter(ShapingCurve *sc, int sCount, int pCount, double noise);
		~LinFitter();
	private:
		double *guess_t;
		gsl_vector *y_vec, *a_vec;
		gsl_matrix *sc_mat, *coeff_mat;
		double *fit_t, *fit_t_err;
	public:
		double doLinFit(double *times, double *par=NULL); // does linear fit and returns chisq; writes fitted parameters to *par
		void initMinuit();
		void doMinuitFit();
		void plotFit(Event *evt, const char *name); //plot event, fit and chisq dependence on each fitted t_0
		void makeGuess();//guesses initial conditions based on steepest slopes between samples
		void print_guess();
		void setVerbosity(int verbosity);
		double Evaluate(double *x, double *p); //wrapped doLinFit for use with TF2 constructor
		void plotFCN(Event *evt, const char *name, int res, const char *style, double x_min, double x_max, double y_min, double y_max); //makes a contour plot of chisq; assumes 2 hits and 2-peak fit
		void doGridFit(); //does a grid minimization of doLinFit; also writes fit result to guess_t
		void doFit();
};

void LinFCN(Int_t&npar, Double_t*gin, Double_t&f, Double_t*par, Int_t flag); //doLinFit wrapped as a Minuit FCN
Double_t LinFCNCurve(Double_t *x, Double_t *par); //function used to make TF1 for plotFit
#endif
