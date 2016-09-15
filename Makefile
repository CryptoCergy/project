CC=gcc
CFLAGS=-W -Wall -Wextra -Wno-unused-parameter -O3 
LDFLAGS=-W -Wall -Wextra -O3
SRC= $(wildcard *.c) 
OBJ= $(SRC:.c=.o)

EXEC=var

all: $(EXEC)

$(EXEC): $(OBJ) 
	$(CC) -o $@ $(OBJ) $(LDFLAGS) 
	rm *.o

%.o: %.c include/common.h include/%.h
	$(CC) $(CFLAGS) -c $<  

clean :
	rm $(OBJ)
	rm $(EXEC)
