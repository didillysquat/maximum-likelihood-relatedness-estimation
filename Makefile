CXX ?= c++
CFLAGS = -Wall -Werror -fmax-errors=3
CCFLAGS = -pthread -Wall
LIBARG	= -g -std=c++11 -O3
EIGEN = include
INC = -I $(EIGEN)
TARGET	= lcmlkin #Variant split 
FILES = relatedness utils Variant split
SRC	= $(addprefix src/, $(addsuffix .cpp, $(FILES))) #src/filevercmp.c

$(TARGET): $(SRC)
	$(CXX) $(CCFLAGS) $(SRC) $(LIBARG) $(INC) -o $@

clean:
	rm -f $(TARGET)
