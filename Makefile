SRCDIR := CC
INCDIR := include
#OBJDIR := obj

HEADERS  :=  $(INCDIR)/WFCTAEvent.h $(INCDIR)/camera.h $(INCDIR)/dumpPack.h
SOURCES  :=  $(SRCDIR)/WFCTAEvent.cc 
OBJS     :=  WFCTAEvent.o WFCTADecode.o WFCTAEventDict.o WFCTADecodeDict.o

all:event.exe status.exe

event.exe: event.o $(OBJS)
	g++ -g -O2 -o $@ $^  `root-config --cflags --libs`

status.exe: status.o $(OBJS)
	g++ -g -O2 -o $@ $^  `root-config --cflags --libs`

event.o:event.cc
	g++ -g -O2 -c $< `root-config --cflags`

status.o:status.cc
	g++ -g -O2 -c $< `root-config --cflags`

WFCTADecodeDict.o:$(INCDIR)/WFCTADecodeDict.cc
	g++ -g -O2 -c $< `root-config --cflags --libs`

WFCTAEventDict.o:$(INCDIR)/WFCTAEventDict.cc
	g++ -g -O2 -c $< `root-config --cflags --libs`

WFCTADecode.o:$(SRCDIR)/WFCTADecode.cc $(INCDIR)/WFCTADecode.h
	g++ -g -O2 -c $^ `root-config --cflags --libs`

WFCTAEvent.o:$(SOURCES) $(HEADERS)
	g++ -g -O2 -c $< `root-config --cflags --libs` 

clean:
	rm *.exe *.o
