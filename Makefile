CC	= g++
CFLAGS = -c -Wall -Werror -fmax-errors=3
LIBARG	= -g -std=c++11 -O3
EIGEN = include
INC = -I $(EIGEN)
TARGET	= relatedness utils Variant split 
SRC	= $(addprefix src/, $(addsuffix .cpp, $(TARGET))) src/filevercmp.c

$(TARGET): $(SRC)
	$(CC) $(SRC) $(LIBARG) $(INC) -o $@

clean:
	rm -f $(TARGET)
