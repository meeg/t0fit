#ifndef ANALYTIC_FITTER_HH
#define ANALYTIC_FITTER_HH
#include "Fitter.hh"

class AnalyticFitter : public Fitter
{
	public:
		AnalyticFitter(ShapingCurve *sc, int sCount, int pCount, double noise);
		~AnalyticFitter();
	private:
		double shapingTime;
		int *sectionStart, *sectionEnd;
		double *section_par, *section_err;
		double *section_covar;
		int nSections;
	public:
		int doAnalyticFit(int start, int end, double &time, double &time_err, double &height, double &height_err, double &covar); //fits to a range of samples
		void mergeSections(); //takes the section fits and subtracts pileup
		void doFitSingle(); //fits to a single peak (deprecated)
		void doFit(); //fits multiple peaks
		void plotFit(Event *evt, const char *name);
};

#endif
