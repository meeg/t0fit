CPP = g++

CPPFLAGS = -g -I$(ROOTSYS)/include -I$(ROOTSYS)/include/root
CXXFLAGS =

ROOTLIBS = `$(ROOTSYS)/bin/root-config --glibs`

LIBS = -lgsl -lgslcblas $(ROOTLIBS) -lMinuit
SOURCES = apvtime.cc Event.cc ShapingCurve.cc Samples.cc Fitter.cc SmoothShapingCurve.cc LinFitter.cc AnalyticFitter.cc MinuitFitter.cc
OBJS = $(SOURCES:.cc=.o)
#DEPS = $(SOURCES:.cpp=.d)


%.o: %.cc
	@echo "[CPP]   $<"
	@$(CPP) -c $< $(CPPFLAGS) $(CXXFLAGS) -o $@


all: apvtime

apvtime: $(OBJS)
	@echo "[LINK]  apvtime"
	@$(CPP) $(LIBS) -o apvtime $(OBJS)

clean:
	@rm -f $(OBJS) apvtime
