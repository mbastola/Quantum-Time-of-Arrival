CXX ?= g++
CXXFLAGS ?= -std=c++11 -O2 -Wall -Wextra -pedantic

TARGET := qtoa
SRCS := qtoa.cpp common.cpp physics.cpp arrival.cpp plotting.cpp generation.cpp
HEADERS := common.hpp integration.hpp physics.hpp arrival.hpp plotting.hpp generation.hpp

.PHONY: all run quick clean

all: $(TARGET)

$(TARGET): $(SRCS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRCS)

run: $(TARGET)
	./$(TARGET) --out output

quick: $(TARGET)
	./$(TARGET) --quick --out output

clean:
	rm -f $(TARGET)
