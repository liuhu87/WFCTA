read.exe: read.o WFCTAEvent.o WFCTAEventDict.o
	g++ -g -O2 -o read.exe read.o WFCTAEventDict.o WFCTAEvent.o  `root-config --cflags --libs`

read.o:read.cc
	g++ -g -O2 -c read.cc `root-config --cflags`

WFCTAEventDict.o:WFCTAEventDict.cc
	g++ -g -O2 -c WFCTAEventDict.cc `root-config --cflags --libs` 

WFCTAEvent.o:WFCTAEvent.cc
	g++ -g -O2 -c WFCTAEvent.cc `root-config --cflags --libs` 

clean:
	rm *.exe *.o
