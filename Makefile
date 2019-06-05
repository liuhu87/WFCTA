SRCDIR := CC
INCDIR := include
OBJDIR := obj

HEADERS  :=  $(INCDIR)/WFCTAEvent.h $(INCDIR)/WFCTADecode.h
SOURCES  :=  $(SRCDIR)/WFCTAEvent.C $(SRCDIR)/WFCTADecode.C event.C status.C
OBJS     :=  $(OBJDIR)/WFCTAEvent.o $(OBJDIR)/WFCTADecode.o $(OBJDIR)/dictionary.o

CXX      :=`root-config --cxx`
DEFINES  := -I. -I$(INCDIR) -I$(OBJDIR) `root-config --cflags`
CXXFLAGS := -O3 -fPIC
LDFLAGS  := `root-config --libs`

all:event.exe status.exe

event.exe: $(OBJDIR)/event.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

status.exe: $(OBJDIR)/status.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/event.o: event.C
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/status.o: status.C
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.C
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/%.o:$(OBJDIR)/%.C
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/dictionary.C: $(HEADERS) $(INCDIR)/linkdef.h
	rootcint -f $@ -c $(DEFINES) -p $^

clean:
	rm -rfv *.exe $(OBJDIR)/*
