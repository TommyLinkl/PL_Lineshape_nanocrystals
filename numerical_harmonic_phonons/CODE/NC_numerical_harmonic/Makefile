CXX = gcc
CXXFLAGS:= -O3
OBJ = main.o init.o read.o rand.o symplectic.o force.o energy.o correlation.o rand_r.o rand_v.o nerror.o
LIBS = -lm

%.o : %.c
	$(CXX) -c -o $@ $< $(CXXFLAGS)

md_run : $(OBJ)
	$(CXX) -o md_run.exe $^ $(CXXFLAGS) $(LIBS)

clean :
	rm -r $(OBJ) *.exe 
