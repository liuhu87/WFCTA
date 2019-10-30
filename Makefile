SRCDIR := CC
INCDIR := include
LIBDIR := lib
OBJDIR := obj

CXX  :=`root-config --cxx`
LD   :=`root-config --ld`
ARCH :=`root-config --arch`

HEADERS  := $(INCDIR)/common.h
HEADERS  += $(INCDIR)/CorsikaIO.h $(INCDIR)/CorsikaChain.h $(INCDIR)/CorsikaEvent.h
HEADERS  += $(INCDIR)/WFTelescope.h $(INCDIR)/WFMirror.h $(INCDIR)/WFCone.h $(INCDIR)/WFCamera.h $(INCDIR)/Readparam.h $(INCDIR)/WFCTAMCEvent.h  
HEADERS  += $(INCDIR)/WFCTALedEvent.h
HEADERS  += $(INCDIR)/WFCTALaserEvent.h
HEADERS  += $(INCDIR)/WFCTAEvent.h
HEADERS  += $(INCDIR)/WFCTADecode.h
HEADERS  += $(INCDIR)/EventNtuple.h 
HEADERS  += $(INCDIR)/FluxModel.h
HEADERS  += $(INCDIR)/Cloud.h
HEADERS  += $(INCDIR)/LHChain.h
HEADERS  += $(INCDIR)/Laser.h
HEADERS  += $(INCDIR)/StatusDB.h
HEADERS  += $(INCDIR)/ReadTrack.h
HEADERS  += $(INCDIR)/ShowerPlot.h

SOURCES  := $(SRCDIR)/common.C
SOURCES  += $(SRCDIR)/CorsikaIO.C $(SRCDIR)/CorsikaChain.C $(SRCDIR)/CorsikaEvent.C 
SOURCES  += $(SRCDIR)/WFTelescope.C $(SRCDIR)/WFMirror.C $(SRCDIR)/WFCone.C $(SRCDIR)/WFCamera.C $(SRCDIR)/Readparam.C $(SRCDIR)/WFCTAMCEvent.C
SOURCES  += $(SRCDIR)/WFCTALedEvent.C
SOURCES  += $(SRCDIR)/WFCTALaserEvent.C
SOURCES  += $(SRCDIR)/WFCTAEvent.C
SOURCES  += $(SRCDIR)/WFCTADecode.C
SOURCES  += $(SRCDIR)/EventNtuple.C 
SOURCES  += $(SRCDIR)/FluxModel.C
SOURCES  += $(SRCDIR)/Cloud.C
SOURCES  += $(SRCDIR)/LHChain.C
SOURCES  += $(SRCDIR)/Laser.C
SOURCES  += $(SRCDIR)/StatusDB.C
SOURCES  += $(SRCDIR)/ReadTrack.C
SOURCES  += $(SRCDIR)/ShowerPlot.C
SOURCES  += event.C
SOURCES  += eventSort.C
SOURCES  += eventSortMerge.C
SOURCES  += eventYMJ.C
SOURCES  += status.C
SOURCES  += read.C
SOURCES  += dosim.C
SOURCES  += showcloudmap.C
SOURCES  += dolasersim.C

OBJS     := $(OBJDIR)/common.o
OBJS     += $(OBJDIR)/CorsikaIO.o $(OBJDIR)/CorsikaChain.o $(OBJDIR)/CorsikaEvent.o 
OBJS     += $(OBJDIR)/WFTelescope.o $(OBJDIR)/WFMirror.o $(OBJDIR)/WFCone.o $(OBJDIR)/WFCamera.o $(OBJDIR)/Readparam.o $(OBJDIR)/WFCTAMCEvent.o
OBJS     += $(OBJDIR)/WFCTALedEvent.o
OBJS     += $(OBJDIR)/WFCTALaserEvent.o
OBJS     += $(OBJDIR)/WFCTAEvent.o
OBJS     += $(OBJDIR)/WFCTADecode.o
OBJS     += $(OBJDIR)/EventNtuple.o 
OBJS     += $(OBJDIR)/FluxModel.o
OBJS     += $(OBJDIR)/Cloud.o
OBJS     += $(OBJDIR)/LHChain.o
OBJS     += $(OBJDIR)/Laser.o
OBJS     += $(OBJDIR)/StatusDB.o
OBJS     += $(OBJDIR)/ReadTrack.o
OBJS     += $(OBJDIR)/ShowerPlot.o
OBJS     += $(OBJDIR)/dictionary.o

DEFINES  := -I. -I$(INCDIR) -I$(OBJDIR) `root-config --cflags`

#CXXFLAGS := -O3 -fPIC -qopenmp
CXXFLAGS := -O3 -fPIC
#CXXFLAGS += -D_THIN_

LDFLAGS  := `root-config --libs`
# -lRGL -lEve -lGeom -lMinuit -lTMVA -lXMLIO -lMLP -lTreePlayer -lXrdClient -lGpad -lNet -lHist -lHistPainter -lGraf -lMatrix -lRooFit

all:event.exe eventSort.exe eventSortMerge.exe eventYMJ.exe status.exe read.exe dosim.exe showcloudmap.exe dolasersim.exe $(LIBDIR)/lib.so

event.exe: $(OBJDIR)/event.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

eventSort.exe: $(OBJDIR)/eventSort.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

eventSortMerge.exe: $(OBJDIR)/eventSortMerge.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

eventYMJ.exe: $(OBJDIR)/eventYMJ.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

status.exe: $(OBJDIR)/status.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

read.exe: $(OBJDIR)/read.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

dosim.exe: $(OBJDIR)/dosim.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

showcloudmap.exe: $(OBJDIR)/showcloudmap.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

dolasersim.exe: $(OBJDIR)/dolasersim.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(LIBDIR)/lib.so: $(OBJS)
	mkdir -p $(LIBDIR)
	$(LD) -shared -o $@ $^ $(LDFLAGS)

$(OBJDIR)/event.o: event.C
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/eventSort.o: eventSort.C
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/eventSortMerge.o: eventSortMerge.C
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/eventYMJ.o: eventYMJ.C
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/status.o: status.C
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/read.o: read.C
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/dosim.o: dosim.C
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/showcloudmap.o: showcloudmap.C
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/dolasersim.o: dolasersim.C
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.C
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/%.o: $(OBJDIR)/%.C
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/dictionary.C: $(HEADERS) $(INCDIR)/linkdef.h
	rootcint -f $@ -c $(DEFINES) -p $^

clean:
	rm -rfv *.exe $(OBJDIR)/* $(LIBDIR)/lib.so $(OBJDIR)
