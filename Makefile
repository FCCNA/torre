SHELL = /bin/bash

SRCSUF:=cpp
OBJSUF:=o
DLLSUF:=so

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

CC ?= $(shell root-config --cc)
CXX ?= $(shell root-config --cxx)
CFLAGS += $(shell root-config --libs --cflags) -Wall

SRCDIR = src
BINDIR = bin
LIBDIR = lib

TRACKING_OBJECTS = tdata_track.$(OBJSUF) tracking.$(OBJSUF)
#PROCESS_OBJECTS = tdata.$(OBJSUF) process.$(OBJSUF)
PARSE_OBJECTS = dynamictree.$(OBJSUF)
#SIM_OBJECTS = montecarlo.$(OBJSUF)

TRACK_OBJ = $(patsubst %,$(BINDIR)/%,$(TRACKING_OBJECTS))
#PROC_OBJ = $(patsubst %,$(BINDIR)/%,$(PROCESS_OBJECTS))
PARS_OBJ = $(patsubst %,$(BINDIR)/%,$(PARSE_OBJECTS))
#SIM_OBJ = $(patsubst %,$(BINDIR)/%,$(SIM_OBJECTS))

all: $(PROC_OBJ) $(PARS_OBJ) $(SIM_OBJ) $(TRACK_OBJ)
	$(CXX) $(CFLAGS) -o tracking $(TRACK_OBJ)
#	$(CXX) $(CFLAGS) -o process $(PROC_OBJ)
	$(CXX) $(CFLAGS) -o parse $(PARS_OBJ)
#	$(CXX) $(CFLAGS) -o simulate $(SIM_OBJ)

$(BINDIR)/%.$(OBJSUF): $(SRCDIR)/%.$(SRCSUF) | $(BINDIR)
	$(CXX)  $(CFLAGS) -I$(SRCDIR) -c -o $@ $<

$(BINDIR):
	@[[ -d "$(BINDIR)" ]] || mkdir -p "$(BINDIR)"

.PHONY: clean
clean:
	rm -f $(TRACK_OBJ)
	rm -f $(PARS_OBJ)
#	rm -f $(PROC_OBJ)
#	rm -f $(SIM_OBJ)
	rmdir $(BINDIR)
	rm -f parse
#	rm -f process
#	rm -f simulate
	rm -f tracking

