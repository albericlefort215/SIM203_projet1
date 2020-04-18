
CC=gcc
COPTS= -O3 -Wall


mesh:	mesh.o main_mesh.o lplib3.o libmeshb7.o 
	$(CC) $(COPTS) -o mesh mesh.o  main_mesh.o libmeshb7.o lplib3.o -lpthread -lm

lplib3.o :	lplib3.c lplib3.h 
	$(CC) -c $(COPTS)  -I. lplib3.c

libmeshb7.o :	libmeshb7.c libmeshb7.h 
	$(CC) -c $(COPTS)  -I. libmeshb7.c  

mesh.o :	mesh.c mesh.h libmeshb7.h 
	$(CC) -c $(COPTS)  -I. mesh.c

main_mesh.o :	main_mesh.c mesh.h
	$(CC) -c $(COPTS)  -I. main_mesh.c  

clean :
	-rm mesh mesh.o main_mesh.o lplib3.o libmeshb7.o 
