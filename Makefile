CC = gcc
CFLAGS = -pthread -mavx -std=c11
TARGET = parallel_distance

all: $(TARGET)

$(TARGET): parallel_distance.c
	$(CC) $(CFLAGS) -o $(TARGET) parallel_distance.c -lm

clean:
	rm -f $(TARGET)

