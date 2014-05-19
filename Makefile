CXXFLAGS =	-O2 -g -Wall -fmessage-length=0 -std=c++0x -fopenmp -lpthread

OBJS =		NOSE-rates-stochastic.o

LIBS =

TARGET =	NOSE-rates-stochastic

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS) -fopenmp

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
