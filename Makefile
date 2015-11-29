#===============================================================================
# User Options
#===============================================================================

PROGRAM  = 3DKBA
CC       = g++
DEBUG    = no
OPTIMIZE = yes
PROFILE  = no

#===============================================================================
# Object Files
#===============================================================================

OBJECTS = \
  main.o \
  auxiliary_function.o \
  miniapp.o \
  mesh.o \
  eas.o \
  esa.o \
  sae.o \
  sea.o \
  ase.o \
  aes.o \
  eas_mod.o

#===============================================================================
# COMPILER FLAGS
#===============================================================================

ifeq ($(CC), icpc)
	CCFLAGS += -openmp
else
	CCFLAGS += -fopenmp
endif

ifeq ($(OPTIMIZE),yes)
  CCFLAGS += -O3
else
  CCFLAGS += -O0
endif

ifeq ($(DEBUG),yes)
  CCFLAGS += -g
endif

#===============================================================================
# Targets
#===============================================================================

$(PROGRAM): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(CCFLAGS) $(LDFLAGS) $(LIB) $(INCLUDE)

clean:
	@rm -f *.o $(PROGRAM)

neat:
	@rm -f *.o

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .cc .o
.PHONY: clean neat

%.o: %.cc
	$(CC) $(CCFLAGS) $(INCLUDE) -c $<

#===============================================================================
# Dependencies
#===============================================================================

main.o: auxiliary_function.o
main.o: miniapp.o 
#main.o: eas.o 
#main.o: esa.o 
#main.o: sae.o 
#main.o: sea.o 
#main.o: ase.o 
#main.o: aes.o 
#main.o: eas_mod.o

miniapp.o: auxiliary_function.o 
miniapp.o: eas.o 
miniapp.o: esa.o 
miniapp.o: sae.o 
miniapp.o: sea.o 
miniapp.o: ase.o 
miniapp.o: aes.o 
miniapp.o: eas_mod.o

eas.o: auxiliary_function.o
esa.o: auxiliary_function.o
aes.o: auxiliary_function.o
ase.o: auxiliary_function.o
sae.o: auxiliary_function.o
sea.o: auxiliary_function.o
eas_mod.o: auxiliary_function.o