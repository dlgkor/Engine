#all:
#	g++ -I src/include -L src/lib -o main main.cpp -lmingw32 -lSDL2main -lSDL2 -lSDL2_image

CC = g++
Target = main
SRCS = $(wildcard *.cpp) $(wildcard *.c)
OBJS = $(SRCS:.cpp=.o)
INCLUDE_PATH = -I src/include
LIBRARY_PATH = -L src/lib
LINKER_FLAG = -lmingw32 -lSDL2main -lSDL2 -lSDL2_image

STATICLIBS = -static-libgcc -static-libstdc++


all: $(Target)
	$(CC) $(STATICLIBS) -o $(Target) $(OBJS) $(INCLUDE_PATH) $(LIBRARY_PATH) $(LINKER_FLAG) 

$(Target):
	$(CC) -c $(SRCS) $(INCLUDE_PATH) $(LIBRARY_PATH) $(LINKER_FLAG)

clean : 
	-del -fR $(OBJS)