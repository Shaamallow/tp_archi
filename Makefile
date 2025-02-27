CC = gcc
CFLAGS = -pthread -mavx -lm -std=c11 -Wall
TARGET = main

all: $(TARGET)

$(TARGET): main.o utils.o
	$(CC) $(CFLAGS) -o $(TARGET) main.o utils.o

main.o: main.c utils.h
	$(CC) $(CFLAGS) -c main.c

utils.o: utils.c utils.h
	$(CC) $(CFLAGS) -c utils.c

clean:
	rm -f $(TARGET) *.o

