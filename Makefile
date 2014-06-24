###########################################################################
# MUPAGE Makefile
# G.Carminati - carminati@bo.infn.it
###########################################################################
SHELL = /bin/csh

# Needed libraries and include file folders
package = mupage
# include ROOT libraries
ROOTLIBS  = `${ROOTSYS}/bin/root-config --libs`
ROOTINCL  = `${ROOTSYS}/bin/root-config --cflags`

C++    = g++
C++OPT = -g -O2 -Wall -Wno-deprecated
RM = rm

################################## MUPAGE #################################
# folders list
SRC     = src
INCLUDE = inc
BUILD   = Linux
EVENTS  = evt
INFO    = livetime

# include file list
INCFILES  = -I$(INCLUDE) $(ROOTINCL)
LINKFILES = $(ROOTLIBS)

#----------------------------- Suffix Rules -------------------------------
# set up C++ suffixes and relationship between file.cc and file.o

 .cc.o:
	$(C++) -c $(C++OPT) $(INCFILES) -o $(BUILD)/$@ $<

 .o :
	$(C++) $(C++OPT) $(INCFILES) $< -o $@

#--------------------------- File Dependencies ----------------------------
SOURCES:= $(wildcard $(SRC)/*.cc)
OBJECTS:= $(notdir $(SOURCES:.cc=.o))

vpath %.cc $(SRC)
vpath %.o $(BUILD)

all: config ${package}.exe

${package}.exe: $(OBJECTS)
	$(C++) $(C++OPT) -o ${package}.exe $(INCFILES) \
	$(addprefix $(BUILD)/,$(OBJECTS)) $(LINKFILES)

clean:
	$(RM) -rf $(BUILD) ${package}.exe

config:
	mkdir -p $(BUILD) $(EVENTS) $(INFO)
