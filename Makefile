CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -O3 -march=native
TARGET = pricing

all: $(TARGET)

$(TARGET): main.cpp Options.h Market.h Pricers.h
	$(CXX) $(CXXFLAGS) main.cpp -o $(TARGET) -lm

clean:
	rm -f $(TARGET)

run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run
