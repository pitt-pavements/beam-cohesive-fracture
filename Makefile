#Makefile for ILLITHERM-FV

PROGRAM = secant-base-3
SOURCE = misc.f PPR_model_educational.f matrices.f crack-element-x.f main.f inverse1.f inverse2.f
OBJECTS = misc.o PPR_model_educational.o matrices.o crack-element-x.o main.o inverse1.o inverse2.o
ARCHIVE = secant-base.tar

CC = gfortran
DEPS = ABA_PARAM.INC
CFLAGS = -fbacktrace
OPTIONS = -g 

%.o: %.f $(DEPS)
	$(CC) $(OPTIONS) -c -o $@ $< $(CFLAGS)

$(PROGRAM):	$(OBJECTS)
	$(CC) $(OPTIONS) -o $(PROGRAM).out $(OBJECTS) $(CFLAGS)
	
clean: 
	rm -f $(OBJECTS) $(PROGRAM).out *.tar *.plt
	
rebuild:
	clean $(PROGRAM)
	
archive:
	@echo Creating code archive tar-file $(ARCHIVE) ...
	tar cf $(ARCHIVE) $(SOURCE) Makefile
	@ls -l $(ARCHIVE)

help:
	@echo Try:
	@echo make $(PROGRAM) .... to build the program named $(PROGRAM)
	@echo make clean .... to clean up, removing object files and program $(PROGRAM)
	@echo make archive .... to make an archive tar file you can transfer or submit