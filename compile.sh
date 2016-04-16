#!/bin/sh
g++ -Wall -o "./sequencer_simulator" "./sequencer_simulator.cpp" -std=c++14  -march=native  -O3 -pthread
g++ -Wall -o "./smile_sequencer" "./smile_sequencer.cpp" -std=c++14  -march=native  -O3 -pthread
