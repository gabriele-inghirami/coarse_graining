# compiler to be used (GNU gcc or Intel icc)
comp=gcc
#comp=icc

# linker flags
ldflags=-lm

# gsl linker flags for all target, which includes to_landau.exe
gslflags=-lgsl -lblas

# optimization or debugging flags

# safe optimization flag, valid for all processors
#opt=-O2 

#gcc 8.3.0 optimization flag for the CSC GOETHE cluster in Frankfurt am Main, partition general1, or CSC Puhti cluster in Finland (skylake architecture)
#opt=-O3 -march=skylake

#gcc 8.3.0 optimization flag for the CSC GOETHE cluster in Frankfurt am Main, partition general2 (broadwell architecture)

#gcc 8.2.0 optimization flag for the CSC GOETHE cluster in Frankfurt am Main, partition fuchs
#opt=-O3 -march=ivybridge

#gcc 4.8.5 optimization flag for the CSC GOETHE cluster in Frankfurt am Main, partition general2 (broadwell architecture)
#opt=-O3 -march=avx2
#opt=-O3 -mavx

#icc optimization flag for the CSC GOETHE cluster in Frankfurt am Main, partition general1, or CSC Puhti cluster in Finland (skylake architecture)
#opt=-O3 -xCORE-AVX512 -ipo

#icc optimization flag for the CSC GOETHE cluster in Frankfurt am Main, partition general2 (broadwell architecture)
#opt=-O3 -xAVX -ipo

#gcc optimization flag that assumes that the computers used to compile and run the program have the same architecture
opt=-O3 -march=native

# debugging flags for GNU gdb (and also ddd, Affinic Debugger, Allinea DDT,...)
#opt=-g -ggdb -O0

# debugging flag for Absoft FX3 debugger
# opt=-g -gdwarf-3 -O0

#additional debugging flags to detect floating point exceptions
#opt+=-fsignaling-nans

#additional flag to ignore warning for unused returning values
#opt+=-Wno-unused-result

#all: cg.exe to_text.exe to_2D.exe to_landau.exe
all: cg.exe to_text.exe to_2D.exe 

cg.exe: main.c calculate.c tools.c io.c particles.c
	$(comp) $? -o $@ $(ldflags) $(opt)
#$(comp) $? -o $@ -L /home/hireaction/inghirami/mylib/lib $(ldflags) $(opt)

to_text.exe: to_text.c
	$(comp) $? -o $@ $(ldflags) $(opt)

to_2D.exe: to_2D.c 
	$(comp) $? -o $@ $(ldflags) $(opt)

# target: to delete the products of the compilation (executable, object files, modules...)
clean:
		  \rm *.exe


