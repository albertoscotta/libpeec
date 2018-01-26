LIB = libpeec.a
AR = ar
AR_OPTS = r
CC = gcc
CC_OPTS = -g -Wall --std=c99 -fopenmp -D DEBUG #-D MACKENZIE -D RECT_APPROX  
OBJS = coeff_mats.o mesh.o constants.o utils/lin_algebra.o \
	utils/integrals.o coeff_mats/rp_mat.o coeff_mats/pp_mat.o \
	coeff_mats/lp_mat.o
MAKEFILE = Makefile
INCLUDE = coeff_mats.h mesh.h constants.h utils/lin_algebra.h \
	utils/integrals.h coeff_mats/rp_mat.h coeff_mats/lp_mat.h \
	coeff_mats/pp_mat.h
LIBDIR = /usr/local/lib
INCLUDEDIR = /usr/local/include/peec
DOCDIR = ../doc

../$(LIB): $(OBJS)
	$(AR) $(AR_OPTS) $@ $(OBJS)

%.o: %.c $(MAKEFILE) $(INCLUDE)
	$(CC) $(CC_OPTS) -c $< -o $@

install:
	cp ../$(LIB) $(LIBDIR)
	mkdir $(INCLUDEDIR)
	cp --parents $(INCLUDE) $(INCLUDEDIR)

uninstall:
	rm $(LIBDIR)/$(LIB)
	rm -r $(INCLUDEDIR)

doc:
	mkdir -p $(DOCDIR)
	doxygen

clean:
	rm -f ../$(LIB) $(OBJS)
	rm -rf $(DOCDIR)

.PHONY: install uninstall doc clean
