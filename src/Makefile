#***************************************************************************
#*   Copyright (C) 2008 by Antonio Jesús Adsuar Gómez                      *
#*   adsuar@lsi.upc.edu                                                    *
#*                                                                         *
#*   This file is part of CORTEX                                           *
#*   CORTEX is free software; you can redistribute it and/or modify        *
#*   it under the terms of the GNU General Public License as published by  *
#*   the Free Software Foundation; either version 2 of the License, or     *
#*   (at your option) any later version.                                   *
#*                                                                         *
#*   CORTEX is distributed in the hope that it will be useful,             *
#*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
#*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
#*   GNU General Public License for more details.                          *
#*                                                                         *
#*   You should have received a copy of the GNU General Public License     *
#*   along with this program; if not, write to the                         *
#*   Free Software Foundation, Inc.,                                       *
#*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
#***************************************************************************/

CC = c++
CFLAGS_debug= #-Wall -g -Wno-deprecated
CFLAGS= -Wall -Wno-deprecated
OBJ = fsdata/fsdata.o fsaprs/fsaprs.o fspata/fspata.o fspd/fspd.o fsvotes/fsvotes.o fsmatrix/fsmatrix.o fsdistance/fsdistance.o fsevaluation/fsevaluation.o fsvector/fsvector.o fsset/fsset.o main.o
LINKOBJ = fsdata/fsdata.o fsaprs/fsaprs.o fspata/fspata.o fspd/fspd.o fsvotes/fsvotes.o fsmatrix/fsmatrix.o fsdistance/fsdistance.o fsevaluation/fsevaluation.o fsvector/fsvector.o fsset/fsset.o main.o
BIN = ../bin/cortexaprs
RM = rm -rf

all:$(OBJ)
	$(CC) $(LINKOBJ) -o $(BIN)

main.o: main.cc
	$(CC) $(CFLAGS_debug) -c main.cc -o main.o
	
fsset/fsset.o: fsset/fsset.cc
	$(CC) $(CFLAGS_debug) -c fsset/fsset.cc -o fsset/fsset.o
fsvector/fsvector.o: fsvector/fsvector.cc
	$(CC) $(CFLAGS_debug) -c fsvector/fsvector.cc -o fsvector/fsvector.o
fsdata/fsdata.o: fsdata/fsdata.cc
	$(CC) $(CFLAGS_debug) -c fsdata/fsdata.cc -o fsdata/fsdata.o
fsmatrix/fsmatrix.o: fsmatrix/fsmatrix.cc
	$(CC) $(CFLAGS_debug) -c fsmatrix/fsmatrix.cc -o fsmatrix/fsmatrix.o
fsdistance/fsdistance.o: fsdistance/fsdistance.cc
	$(CC) $(CFLAGS_debug) -c fsdistance/fsdistance.cc -o fsdistance/fsdistance.o
fsevaluation/fsevaluation.o: fsevaluation/fsevaluation.cc
	$(CC) $(CFLAGS_debug) -c fsevaluation/fsevaluation.cc -o fsevaluation/fsevaluation.o
fsaprs/fsaprs.o: fsaprs/fsaprs.cc
	$(CC) $(CFLAGS_debug) -c fsaprs/fsaprs.cc -o fsaprs/fsaprs.o
fspata/fspata.o: fspata/fspata.cc
	$(CC) $(CFLAGS_debug) -c fspata/fspata.cc -o fspata/fspata.o
fspd/fspd.o: fspd/fspd.cc
	$(CC) $(CFLAGS_debug) -c fspd/fspd.cc -o fspd/fspd.o
fsvotes/fsvotes.o: fsvotes/fsvotes.cc
	$(CC) $(CFLAGS_debug) -c fsvotes/fsvotes.cc -o fsvotes/fsvotes.o

#------------------------ METODOS AUXILIARES DE MAKEFILE --------------------

tgz:
	tar -zcvf CORTEX_distrib.tgz *

copy:
	cp CORTEX_distrib.tgz $(HOME)/Distrib

clean:
	$(RM) *~ $(LINKOBJ) DEADJOE core $(BIN) *.tgz
 
