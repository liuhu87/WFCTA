SRCDIR := CC
INCDIR := include
OBJDIR := obj

HEADERS  :=  $(INCDIR)/WFCTAEvent.h $(INCDIR)/WFCTADecode.h $(INCDIR)/camera.h $(INCDIR)/dumpPack.h
SOURCES  :=  $(SRCDIR)/WFCTAEvent.C $(SRCDIR)/WFCTADecode.C
OBJS     :=  $(OBJDIR)/WFCTAEvent.o $(OBJDIR)/WFCTADecode.o WFCTAEventDict.o WFCTADecodeDict.o

DEFINES  := -I. -I$(INCDIR) -I$(OBJDIR) `root-config --cflags`
CXXFLAGS := -O3 -fPIC

all:event.exe status.exe

event.exe: event.o $(OBJS)
	g++ -g -O2 -o $@ $^  `root-config --cflags --libs`

status.exe: status.o $(OBJS)
	g++ -g -O2 -o $@ $^  `root-config --cflags --libs`

event.o:event.C
	g++ -g -O2 -c $< `root-config --cflags`

status.o:status.C
	g++ -g -O2 -c $< `root-config --cflags`

WFCTADecodeDict.o:$(INCDIR)/WFCTADecodeDict.C
	g++ -g -O2 -c $< `root-config --cflags --libs`

WFCTAEventDict.o:$(INCDIR)/WFCTAEventDict.C
	g++ -g -O2 -c $< `root-config --cflags --libs`

$(OBJDIR)/%.o: $(SRCDIR)/%.C
	mkdir -p $(OBJDIR)
	g++ $(CXXFLAGS) $(DEFINES) -c $^ -o $@

$(OBJDIR)/dictionary.C

#WFCTADecode.o:$(SRCDIR)/WFCTADecode.C $(INCDIR)/WFCTADecode.h
#	g++ -g -O2 -c $^ `root-config --cflags --libs`
#
#WFCTAEvent.o:$(SOURCES) $(HEADERS)
#	g++ -g -O2 -c $< `root-config --cflags --libs` 

clean:
	rm *.exe *.o
