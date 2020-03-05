SRCDIR := CC
INCDIR := include
LIBDIR := lib
OBJDIR := obj

CXX  :=`root-config --cxx`
LD   :=`root-config --ld`
ARCH :=`root-config --arch`

HEADERS  += $(INCDIR)/WFCTALedEvent.h
HEADERS  += $(INCDIR)/WFCTAEvent.h
HEADERS  += $(INCDIR)/WFCTAMerge.h
HEADERS  += $(INCDIR)/WFCTADecode.h
HEADERS  += $(INCDIR)/LHChain.h

SOURCES  += $(SRCDIR)/WFCTALedEvent.C
SOURCES  += $(SRCDIR)/WFCTAEvent.C
SOURCES  += $(SRCDIR)/WFCTAMerge.C
SOURCES  += $(SRCDIR)/WFCTADecode.C
SOURCES  += $(SRCDIR)/LHChain.C
SOURCES  += event.C
SOURCES  += eventSort.C
SOURCES  += eventSortMerge.C
SOURCES  += eventYMJ.C
SOURCES  += status.C
SOURCES  += read.C

OBJS     += $(OBJDIR)/WFCTALedEvent.o
OBJS     += $(OBJDIR)/WFCTAEvent.o
OBJS     += $(OBJDIR)/WFCTAMerge.o
OBJS     += $(OBJDIR)/WFCTADecode.o
OBJS     += $(OBJDIR)/LHChain.o
OBJS     += $(OBJDIR)/dictionary.o

DEFINES  := -I. -I$(INCDIR) -I$(OBJDIR) `root-config --cflags`

#CXXFLAGS := -O3 -fPIC -qopenmp
CXXFLAGS := -O3 -fPIC #-shared
#CXXFLAGS += -D_FCNTL_H
#CXXFLAGS += -D_THIN_

LDFLAGS  := `root-config --libs`
# -lRGL -lEve -lGeom -lMinuit -lTMVA -lXMLIO -lMLP -lTreePlayer -lXrdClient -lGpad -lNet -lHist -lHistPainter -lGraf -lMatrix -lRooFit

all:event.exe eventSort.exe eventSortMerge.exe eventYMJ.exe status.exe read.exe $(LIBDIR)/lib.so

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

$(LIBDIR)/lib.so: $(OBJS)
	mkdir -p $(LIBDIR)
	cp $(OBJDIR)/dictionary_rdict.pcm $(LIBDIR)/
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
