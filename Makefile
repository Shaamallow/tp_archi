CC = gcc
CFLAGS = -pthread -mavx -lm -std=c11 -Wall
TARGET = main

all: $(TARGET)

$(TARGET): main.o utils.o operation.o
	$(CC) $(CFLAGS) -o $(TARGET) main.o utils.o operation.o

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

utils.o: utils.c utils.h
	$(CC) $(CFLAGS) -c utils.c

operation.o: operation.c operation.h
	$(CC) $(CFLAGS) -c operation.c

clean:
	rm -f $(TARGET) *.o

