ROOT_DIR= ../..
SOURCE	= main.cc
TARGET	= run.exe

CXX	= mpicxx
CXXFLAGS= -Wall -std=c++17 -O3 -fmax-errors=1 -fopenmp
CPPFLAGS= -I$(ROOT_DIR)/src -Iminim/include -DMPIPARALLEL
LDFLAGS = -L.
LDLIBS	= -lm -lminim
DEPS    = libminim.a

ifdef strip_templates
STRIP_TEMPLATES = -fdiagnostics-color=always 2>&1 | sed 's/ \[with [^]]*\]//g'
endif


all: $(TARGET)

$(TARGET): $(SOURCE) $(DEPS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS) $(STRIP_TEMPLATES)

libminim.a:
	git submodule update --init minim
	make -C minim
	ln -s minim/bin/libminim.a ./

clean:
	rm -f *.o $(TARGET) $(DEPS)
	rm -rf data
	$(MAKE) -C minim clean

