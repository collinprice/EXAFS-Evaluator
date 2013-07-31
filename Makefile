CC = gcc
CXX = g++

CXXFLAGS = -c -std=c++0x -pthread -g -Wall -Wextra -pedantic #-DDEBUG=1
CCFLAGS = -c
LDFLAGS = -lpthread

CPP_SRC = \
	main.cpp \
	pdbhelper.cpp \
	pdbatom.cpp \
	genfig/genfig.cpp \
	vmdhelper.cpp \
	ifeffithelper.cpp \
	exafsevaluator.cpp \
	iohelper/iohelper.cpp \
	chromosome.cpp \
	exafsga.cpp

CPP_SRC_OBJS = $(CPP_SRC:.cpp=.o)

CC_SRC = file_read_write.c
CC_SRC_OBJS = $(CC_SRC:.c=.o)

EXECUTABLE = main

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(CPP_SRC_OBJS) $(CC_SRC_OBJS)
	$(CXX) $(CPP_SRC_OBJS) $(CC_SRC_OBJS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

.c.o:
	$(CC) $(CCFLAGS) $< -o $@


clean:
		rm -rf $(CPP_SRC_OBJS) $(CC_SRC_OBJS) $(EXECUTABLE)
