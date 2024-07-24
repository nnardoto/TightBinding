#####################################################################################
# 									openmx2xatu										#
#									===========										#
#																				    #
#			N. N. Batista															#
#			nnardoto@gmail.com														#
#			23/Jul./2024															#
#																					#
#			Adaptation of Makefile from openmx@3.9 code								#
#			https://www.openmx-square.org/											#
#																					#
#####################################################################################
CC = gcc -std=c99 -fcommon

openmx2xatu.x: openmx2xatu.o read_scfout.o
	$(CC) read_scfout.o openmx2xatu.o -o openmx2xatu.x

openmx2xatu.o: openmx2xatu.c
	$(CC) -c openmx2xatu.c -o openmx2xatu.o

read_scfout.o: read_scfout.c read_scfout.h 
	$(CC) -c read_scfout.c -o read_scfout.o

clean:
	rm openmx2xatu.x *.o
