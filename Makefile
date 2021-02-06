########################
# Makefile for denada  #
########################

LIB = -lm
INC = -Iik

OBJECTS = \
	denada.o\
	ik/sequence.o\
	ik/toolbox.o\

APP = denada
SRC = main.c
OBJ = main.o

DATE = $(shell date +\%Y-\%m-\%d)

###########
# Targets #
###########

default:
	make gcc

$(APP): $(OBJ) $(OBJECTS)
	$(CC) -o $(APP) $(CFLAGS) $(OBJ) $(OBJECTS) $(LIB)

clean:
	rm -f *.o $(APP)
	cd ik; make clean

depend: $(OBJECTS:.o=.c)
	gcc $(INC) -MM $^ > $@

tar:
	rm -rf /tmp/$(APP)
	mkdir /tmp/$(APP)
	cp -r * /tmp/$(APP)
	cd /tmp/$(APP); make clean; rm -rf CVS */CVS
	cd /tmp; tar -zcvf $(APP)-$(DATE).tar.gz $(APP)
	rm -rf /tmp/$(APP)


#################
# Architectures #
#################

gcc:
	make $(APP) CC="gcc" CFLAGS="-O2 -Wall"


###################
# Inference Rules #
###################

%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

################
# Dependancies #
################

include depend

