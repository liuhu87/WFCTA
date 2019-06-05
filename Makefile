SRCDIR := CC
INCDIR := include
OBJDIR := obj

HEADERS  :=  $(INCDIR)/WFCTAEvent.h $(INCDIR)/WFCTADecode.h $(INCDIR)/camera.h $(INCDIR)/dumpPack.h
SOURCES  :=  $(SRCDIR)/WFCTAEvent.C $(SRCDIR)/WFCTADecode.C event.C status.C
OBJS     :=  $(OBJDIR)/WFCTAEvent.o $(OBJDIR)/WFCTADecode.o $(OBJDIR)/dictionary.o

DEFINES  := -I. -I$(INCDIR) -I$(OBJDIR) `root-config --cflags`
CXXFLAGS := -O3 -fPIC

all:event.exe status.exe

event.exe: $(OBJDIR)/event.o $(OBJS)
	echo "event.exe"
	g++ -g -O2 -o $@ $^ `root-config --libs`

status.exe: $(OBJDIR)/status.o $(OBJS)
	echo "status.exe"
	g++ -g -O2 -o $@ $^  `root-config --libs`

$(OBJDIR)/event.o: event.C
	echo "event.C"
	mkdir -p $(OBJDIR)
	g++ $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/status.o: status.C
	echo "status.C"
	mkdir -p $(OBJDIR)
	g++ $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.C
	echo "src/.C"
	mkdir -p $(OBJDIR)
	g++ $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/%.o:$(OBJDIR)/%.C
	echo "obj/.C"
	mkdir -p $(OBJDIR)
	g++ $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/dictionary.C: $(HEADERS) $(INCDIR)/linkdef.h
	echo "cint"
	rootcint -f $@ -c $(DEFINES) -p $^

clean:
	rm -rfv *.exe $(OBJDIR)/*
