SRCDIR := CC
INCDIR := include
LIBDIR := lib
OBJDIR := obj

CXX  :=`root-config --cxx`
LD   :=`root-config --ld`
ARCH :=`root-config --arch`

HEADERS  := $(INCDIR)/CorsikaIO.h $(INCDIR)/CorsikaChain.h $(INCDIR)/CorsikaEvent.h
HEADERS  += $(INCDIR)/WFTelescope.h $(INCDIR)/WFMirror.h $(INCDIR)/WFCone.h $(INCDIR)/WFCamera.h $(INCDIR)/Readparam.h $(INCDIR)/WFCTAMCEvent.h  
HEADERS  += $(INCDIR)/WFCTALedEvent.h
HEADERS  += $(INCDIR)/WFCTALaserEvent.h
HEADERS  += $(INCDIR)/WFCTAEvent.h
HEADERS  += $(INCDIR)/WFCTADecode.h
HEADERS  += $(INCDIR)/EventNtuple.h 
HEADERS  += $(INCDIR)/FluxModel.h

SOURCES  := $(SRCDIR)/CorsikaIO.C $(SRCDIR)/CorsikaChain.C $(SRCDIR)/CorsikaEvent.C 
SOURCES  += $(SRCDIR)/WFTelescope.C $(SRCDIR)/WFMirror.C $(SRCDIR)/WFCone.C $(SRCDIR)/WFCamera.C $(SRCDIR)/Readparam.C $(SRCDIR)/WFCTAMCEvent.C
SOURCES  += $(SRCDIR)/WFCTALedEvent.C
SOURCES  += $(SRCDIR)/WFCTALaserEvent.C
SOURCES  += $(SRCDIR)/WFCTAEvent.C
SOURCES  += $(SRCDIR)/WFCTADecode.C
SOURCES  += $(SRCDIR)/EventNtuple.C 
SOURCES  += $(SRCDIR)/FluxModel.C
SOURCES  += event.C
SOURCES  += status.C
SOURCES  += read.C

OBJS     := $(OBJDIR)/CorsikaIO.o $(OBJDIR)/CorsikaChain.o $(OBJDIR)/CorsikaEvent.o 
OBJS     += $(OBJDIR)/WFTelescope.o $(OBJDIR)/WFMirror.o $(OBJDIR)/WFCone.o $(OBJDIR)/WFCamera.o $(OBJDIR)/Readparam.o $(OBJDIR)/WFCTAMCEvent.o
OBJS     += $(OBJDIR)/WFCTALedEvent.o
OBJS     += $(OBJDIR)/WFCTALaserEvent.o
OBJS     += $(OBJDIR)/WFCTAEvent.o
OBJS     += $(OBJDIR)/WFCTADecode.o
OBJS     += $(OBJDIR)/EventNtuple.o 
OBJS     += $(OBJDIR)/FluxModel.o
OBJS     += $(OBJDIR)/dictionary.o

DEFINES  := -I. -I$(INCDIR) -I$(OBJDIR) `root-config --cflags`

#CXXFLAGS := -O3 -fPIC -qopenmp
CXXFLAGS := -O3 -fPIC
#CXXFLAGS += -D_THIN_

LDFLAGS  := `root-config --libs`
# -lRGL -lEve -lGeom -lMinuit -lTMVA -lXMLIO -lMLP -lTreePlayer -lXrdClient -lGpad -lNet -lHist -lHistPainter -lGraf -lMatrix -lRooFit

all:event.exe status.exe read.exe $(LIBDIR)/lib.so

event.exe: $(OBJDIR)/event.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

status.exe: $(OBJDIR)/status.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

read.exe: $(OBJDIR)/read.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(LIBDIR)/lib.so: $(OBJS)
	mkdir -p $(LIBDIR)
	$(LD) -shared -o $@ $^ $(LDFLAGS)

$(OBJDIR)/event.o: event.C
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/status.o: status.C
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/read.o: read.C
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
