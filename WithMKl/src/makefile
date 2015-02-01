# New Makefile that automatically depends itself 
#--------------------------------------------------------------
# -warn with no options disables all warning messages
# -WB turns a compile time bounds check into warning
# LD specifies that a program should be linked as a dynamic link # library
# libs tells the compiler which type of run-time library to link # to
# -mp enables improved floating point consistency
# O0 disables all optimizations. Also implies -mp. So #intermediate floating-point results are evaluated at extended #precision
# O1 enables optimizations for speed
# O2 further enables optimizations for speed
# O3 is for more aggressive optimizations
# fast maximises speed across the entire program
# check (-C or -CB) checks for certain conditions at run time
# -C is check all
# -CB is check bounds, error is reported if expression is outside 
# dimension of array or length of string. Once debugged, omit  #this option to improve run time
# traceback option, when severe error occurs, source file, #routine name and line number correlation information is #displayed alongwith call stack hexadecimal addresses (program #counter trace), this option increaes size of executable program, #but has no impact on run-time execution speeds
# -fp-stack-check causes significant code generation after every #function/subroutine call to ensure floating point is maintained #in the correct state. Best used only as a debugging aid to find #floating point stack underflow/overflow problems, which can be #otherwise hard to find
# fpe0, program will print an error message if any divide by zero # exceptions occur
#--------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) -c $(F90FLAGS) $(INCLUDE) $< 
F90 = ifort -warn -WB -C -mkl
LD  = $(F90)
INCLUDE =
LIBS =  
PROFILE = -O0
OPTFLAGS = -g
CFLAGS = $(OPTFLAGS) $(PROFILE)
F90FLAGS = $(CFLAGS)
LDFLAGS =   

PROGRAM = fadd2d
F90SRCS = global.f90 main.f90 nrtype.f90 nrutil.f90 remesh_module.f90 \
          prep.f90 elem.f90 \
          assemb.f90 formakf.f90 gauss.f90 shape.f90 solver.f90 \
          post_remesh.f90 proc.f90 volume_solver.f90 pressure_solver.f90 \
          rearge.f90 sif.f90 sifaniso.f90 root.f90 kernel1_c.f90 kernel2.f90 \
          tstress.f90 tstress1.f90 tstress2.f90 tstress2a.f90 tstress3.f90 \
          assemb_tstress.f90 rigidinner4.f90

SRCS = $(F90SRCS) 

OBJS = $(F90SRCS:.f90=.o)

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(LD) -o $@ $(LDFLAGS) $(OBJS) $(LIBS)

elem.o: elem.f90 ekc_2.f90 ekd_2.f90 ek_b2.f90 ek_at2.f90 \
        eatts2.f90 eattip.f90 eat3.f90 ebts2.f90 ebtip.f90 eb3.f90 \
        eat_rigid.f90 eb_rigid.f90
	$(F90) -c $(F90FLAGS) $(INCLUDE) elem.f90

#Dec10_09: add prep.f90 to following line to recompile prep.f90
prep.o: prep.f90 setint.f90
	$(F90) -c $(F90FLAGS) $(INCLUDE) prep.f90

check:
	../test/test.sh

clean:
	rm -f $(OBJS) $(PROGRAM)
	rm *.mod *genmod.f90
tidy:
	rm -f *.BAK *.bak *.CKP *~

undepend:
	rm -f $(OBJS:%.o=.%.d) 

spotless: tidy clean undepend

