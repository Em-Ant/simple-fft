CC = gcc
CFLAGS= -Wall -O3
LDFLAGS= -lm 
OBJS= 	fft_c.o\
		main.o
		
EXEC=	fft_test

	
all: $(OBJS)
	$(CC) -o $(EXEC) $(OBJS) $(LDFLAGS)
	
clean:
	-@rm -f $(EXEC) *~ *.o *.a

exec: all
	./$(EXEC)
	
	 
