CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

OBJS =		 cell.o SimProcess.o parsimu_beta.o

LIBS =

CXX = mpicxx

TARGET =	parsimu_beta

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)
	
cell.o: cell.cpp
	$(CXX) -c -o cell.o cell.cpp
	
SimProcess.o: SimProcess.cpp
	$(CXX) -c -o SimProcess.o SimProcess.cpp

clean:
	rm -f $(OBJS) $(TARGET)
