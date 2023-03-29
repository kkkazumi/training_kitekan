#Makefile
all:estimate estimate_mood 

#name of program and the project name
TARGET1 = estimate
TARGET2 = estimate_mood
OBJS1 = func_new.o pseudo_inverse.o estimate.o
OBJS2 = func_new.o pseudo_inverse.o est_mental.o estimate_mood.o

#macro
CXX = g++
CXXFLAGS = -std=c++11 -g
CXXFLAGS +=`pkg-config eigen3 --cflags`

#primali target
$(TARGET1): $(OBJS1)
	$(CXX) -o $@ $(CXXFLAGS) $^

$(TARGET2): $(OBJS2)
	$(CXX) -o $@ $(CXXFLAGS) $^

#suffix rule
.cpp .o:
	$(CXX) $(CXXFLAGS) -c $<

#definition except the for suffix rule
.SUFFIXES: .cpp .o

#delete
.PHONY: clean
clean:
	rm -rf $(TARGET1) $(TARGET2) *.o
