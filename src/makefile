objects = CalculateFlux.o CramSolver.o ErrorMessage.o GaussLaguerreQuadrature.o \
GaussLegendreQuadrature.o ModecInitialization.o ModecOutput.o ModecProcedure.o \
NuclList.o ReadCouple.o ReadDepthLib.o ReadMainLib.o RunMODEC.o SetOnlineProcess.o \
tinyxml2.o TTASolver.o

CXX = g++
CXXFLAGS = -std=c++11 -O3

MODEC : $(objects)
	$(CXX) $(CXXFLAGS) -o MODEC $(objects)

CalculateFlux.o : NuclList.h
CramSolver.o : SolveTrans.h
ErrorMessage.o : ErrorMessage.h
GaussLaguerreQuadrature.o : IntegralMethods.h
GaussLegendreQuadrature.o : IntegralMethods.h
ModecInitialization.o : ModecClass.h
ModecOutput.o : ModecClass.h
ModecProcedure.o : ModecClass.h
NuclList.o : NuclList.h
ReadCouple.o : ModecClass.h
ReadDepthLib.o : ModecClass.h
ReadMainLib.o : ModecClass.h
RunMODEC.o : ModecClass.h
SetOnlineProcess.o : ModecClass.h
tinyxml2.o : tinyxml2.h
TTASolver.o : SolveTrans.h

.PHONY : clean
clean :
	-rm $(objects)