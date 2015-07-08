CXXFLAGS =	-O0 -g -Wall -fmessage-length=0

OBJS =		 timers.o cell.o SimProcess.o parsimu_beta.o

LIBS =

CXX = mpicxx

TARGET =	parsimu_beta

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

timers.o: timers.cpp
	$(CXX) -c -o timers.o timers.cpp
	
cell.o: cell.cpp
	$(CXX) -c -o cell.o cell.cpp
	
SimProcess.o: SimProcess.cpp
	$(CXX) -c -o SimProcess.o SimProcess.cpp

clean:
	rm -f $(OBJS) $(TARGET)
