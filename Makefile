CC = /usr/bin/g++	
DEBUGFLAGS = -Wall
OPTFLAGS = -o3

RUN_DIR=./run
RAT_OUTPUT=$(RUN_DIR)/rateliq.x

#folder to include 

INC_LIST= -I ./inc \
	  -I/home/vishal/eigen3/


# Source Folders

SRC_DIR=./src

# Object folders

OBJ_DIR=./obj

# Library folders

LIB_DIR=./lib

# Object list

OBJ_FILES=$(OBJ_DIR)/main.o \
	  $(OBJ_DIR)/dSFMT.o \
	  $(OBJ_DIR)/mathfunc.o \
	  $(OBJ_DIR)/print.o \
	  $(OBJ_DIR)/dist.o\
	  $(OBJ_DIR)/pairlist.o\
	  $(OBJ_DIR)/CalEnerFor.o\
	  $(OBJ_DIR)/pot_forc.o\
	  $(OBJ_DIR)/minimization.o\
	  $(OBJ_DIR)/montecarlo.o\
	  $(OBJ_DIR)/md.o\
	  $(OBJ_DIR)/CalDens.o\
	  $(OBJ_DIR)/simuAnn.o\

# Make Targets
all:$(OBJ_FILES) output

output:$(RAT_OUTPUT)


# build object files
$(OBJ_DIR)/main.o:$(SRC_DIR)/main.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o $(OBJ_DIR)/main.o $(INC_LIST)
$(OBJ_DIR)/dSFMT.o:$(SRC_DIR)/dSFMT.c
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/dSFMT.o $(INC_LIST)
$(OBJ_DIR)/mathfunc.o:$(SRC_DIR)/mathfunc.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/mathfunc.o $(INC_LIST)
$(OBJ_DIR)/print.o:$(SRC_DIR)/print.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/print.o $(INC_LIST)
$(OBJ_DIR)/dist.o:$(SRC_DIR)/dist.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o   $(OBJ_DIR)/dist.o $(INC_LIST)
$(OBJ_DIR)/pairlist.o:$(SRC_DIR)/pairlist.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/pairlist.o $(INC_LIST)
$(OBJ_DIR)/CalEnerFor.o:$(SRC_DIR)/CalEnerFor.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/CalEnerFor.o $(INC_LIST)
$(OBJ_DIR)/pot_forc.o:$(SRC_DIR)/pot_forc.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/pot_forc.o $(INC_LIST)
$(OBJ_DIR)/montecarlo.o:$(SRC_DIR)/montecarlo.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/montecarlo.o $(INC_LIST)
$(OBJ_DIR)/minimization.o:$(SRC_DIR)/minimization.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/minimization.o $(INC_LIST)
$(OBJ_DIR)/md.o:$(SRC_DIR)/md.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/md.o $(INC_LIST)
$(OBJ_DIR)/CalDens.o:$(SRC_DIR)/CalDens.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/CalDens.o $(INC_LIST)
$(OBJ_DIR)/simuAnn.o:$(SRC_DIR)/simuAnn.C
	$(CC) -DDSFMT_MEXP=19937 -c $^ $(OPTFLAGS) -o  $(OBJ_DIR)/simuAnn.o $(INC_LIST)





$(RAT_OUTPUT):$(OBJ_FILES)
	$(CC) $(OPTFLAGS)  $(INC_LIST) -o  $(RAT_OUTPUT) $(OBJ_FILES)


# Clean objects and library
clean:
	$(RM) $(OBJ_FILES)
