HOME = .
DEMO_DIR = $(HOME)/demo
SRC_DIR = $(HOME)/src
LIB_DIR = $(HOME)/lib
LIBS = -lpthread 
INCS = -I$(SRC_DIR) 
CC = mpicxx
CFLAGS=  -O3 -D__x86_64__
#CFLGAS= -D_LINUX_ -O3 -D__x86_64__
LDFLAGS=-O3 #-Wall

#SOURCES = $(wildcard $(SRC_DIR)/*.cc)
#OBJS := $(patsubst %.cc, %.o, $(SOURCES))

.cc.o:
	$(CC) $(CFLAGS) -c $< -o $@

deps:
#	$(CC) -MM $(SOURCES) > Makefile.deps
	echo > Makefile.deps

.PHONY: all clean deps release debug echo doc

all: fgp prep pic pitc ppitc ppic 

prep: deps $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCS) $(LIBS) $(OBJS) $(DEMO_DIR)/demo_prep.cc -o $(DEMO_DIR)/$@
fgp: deps $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCS) $(LIBS) $(OBJS) $(DEMO_DIR)/demo_fgp.cc -o $(DEMO_DIR)/$@
pic: deps $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCS) $(LIBS) $(OBJS) $(DEMO_DIR)/demo_pic.cc -o $(DEMO_DIR)/$@
pitc: deps $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCS) $(LIBS) $(OBJS) $(DEMO_DIR)/demo_pitc.cc -o $(DEMO_DIR)/$@
ppic: deps $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCS) $(LIBS) $(OBJS) $(DEMO_DIR)/demo_ppic.cc -o $(DEMO_DIR)/$@
ppitc: deps $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCS) $(LIBS) $(OBJS) $(DEMO_DIR)/demo_ppitc.cc -o $(DEMO_DIR)/$@

echo:
	@echo CC: $(CC)
	@echo CFLAGS: $(CFLAGS)
	@echo LDFLAGS: $(LDFLAGS)
	@echo SOURCES: $(SOURCES)
	@echo OBJS: $(OBJS)
	@echo LIBS: $(LIBS)

clean:
	cd $(DEMO_DIR)&&rm -fv fgp prep pic pitc ppitc ppic 
	rm -fv $(SRC_DIR)/*.o
	rm -fv $(DEMO_DIR)/*.o
	rm -fv Makefile.deps
clog:
	rm *.log log *.obs
sty:
	astyle --style=linux --indent=spaces=2 -p $(SRC_DIR)/*.h
	astyle --style=linux --indent=spaces=2 -p $(SRC_DIR)/*.cc
doc:
	doxygen Doxyfile
	cd doc/latex && make 
	cp doc/latex/refman.pdf .
