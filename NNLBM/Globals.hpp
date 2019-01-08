#pragma once
#include <assert.h>
#include <string>

typedef double field_t;

// The "rest" direction has been given the index 8, as opposed to the more usual 0.
// This was done for compability with the directions indexes starting at east = 0.
// This might not be necessary.
const enum  CellDirection : int
{
	east = 0,
	northEast = 1,
	north = 2,
	northWest = 3,
	west = 4,
	southWest = 5,
	south = 6,
	southEast = 7,
	rest = 8
};

const static int nDirections = 8;
const static int threeHalfDirection = (3 * nDirections) / 2;

// This function gives the opposite direction of what you put in.
// This is also the reason why I chose to not follow the convension
// having the rest direction be the 0 direction. 
inline int reverseDirectionIndex(const int directionIndex) {
	
	return (directionIndex + threeHalfDirection) % nDirections;
}

const int N_RUN = 10000;
const int N_VELOCITY_PRINT_INTERVAL = 100;
const int N_GRID_X_DIM = 10;
const int N_GRID_Y_DIM = 100;
const field_t F_TAU = 0.77;
const field_t F_BODY_FORCE_X = 1.00e-8;
const field_t F_BODY_FORCE_Y = 0;

const field_t N_VELOCITY_MULTIPLICATION_FACTOR = 1;

const std::string S_BC_BASED = "BCBased/";
const std::string S_SHORT_TESTS = "ShortTests/";

const std::string S_FILE_NAME_BASE = "PoiseuilleFlow";
const std::string S_OUTPUT_DIRECTORY_BASE = "SimulationOutput/";
const std::string S_FLOW_TYPE = "Poiseuille/";
const std::string S_FLOW_SUBCATEGORY = S_SHORT_TESTS;
const std::string S_DIRECTORY_VELOCITY_OUTPUT = S_OUTPUT_DIRECTORY_BASE + S_FLOW_TYPE + S_FLOW_SUBCATEGORY;