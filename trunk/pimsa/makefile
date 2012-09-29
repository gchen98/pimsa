# compile options
#CFLAGS = -Wall -O3  -I/usr/include
CFLAGS = -Wall  -I/usr/include/mysql -g
LINKFLAGS = -L/usr/lib64/mysql -L/usr/lib64 -lm -lgsl -lgslcblas -lmysqlclient 
#LINKFLAGS = -lm -lgsl -lgslcblas

# compiler
CC = g++

# program name
PROGNAME=../bin/sampler

OBJS=db_manager.o pathwaysearch.o ramrepository.o utility.o mysqlrepository.o

$(PROGNAME)  : $(OBJS) 
	$(CC) -o $(PROGNAME) $(OBJS) $(LINKFLAGS)

pathwaysearch.o :pathwaysearch.cpp pathwaysearch.hpp
	$(CC) $(CFLAGS) -c $<

utility.o: utility.cpp utility.hpp
	$(CC) $(CFLAGS) -c $<

db_manager.o: db_manager.cpp db_manager.hpp
	$(CC) $(CFLAGS) -c $<

ramrepository.o: ramrepository.cpp db_manager.hpp ramrepository.hpp
	$(CC) $(CFLAGS) -c $<

mysqlrepository.o: mysqlrepository.cpp mysqlrepository.hpp db_manager.hpp
	$(CC) $(CFLAGS) -c $<


all: $(PROGNAME)

clean :
	rm -fr *.o *.o $(PROGNAME)
