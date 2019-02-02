#pragma once
#include <assert.h>
#include <string>

typedef double field_t;

static const int nDirections = 8;

//#define N_SOLID_CELL_MARGIN_X 1
//#define N_SOLID_CELL_MARGIN_Y 1

#define B_CONSOL_OUT_DENSITY  0;

#define STREAM_TEST     0
#define COUETTE_TEST    1
#define POISEUILLE_TEST 2

#define TEST_TYPE 2

#if TEST_TYPE == STREAM_TEST
const int N_RUN = 20;
const int N_VELOCITY_PRINT_INTERVAL = 1;
const int N_DENSITY_PRINT_INTERVAL = N_VELOCITY_PRINT_INTERVAL;
const int N_GRID_X_DIM_FLUID = 3;
const int N_GRID_Y_DIM_FLUID = 9;
const field_t F_TAU = 1.00;
const field_t F_BODY_FORCE_X = 1.00e-6;
const field_t F_BODY_FORCE_Y = 0;
const field_t F_TOP_PLATE_VELOCITY = 0;
#endif

#if TEST_TYPE == COUETTE_TEST
const int N_RUN = 1000;
const int N_VELOCITY_PRINT_INTERVAL = 100;
const int N_DENSITY_PRINT_INTERVAL = N_VELOCITY_PRINT_INTERVAL;
const int N_GRID_X_DIM_FLUID = 20;
const int N_GRID_Y_DIM_FLUID = 10;
const field_t F_TAU = 0.70;
//const field_t F_BODY_FORCE_X = 1.00e-8;
const field_t F_BODY_FORCE_X = 0;
const field_t F_BODY_FORCE_Y = 0;
const field_t F_TOP_PLATE_VELOCITY = 0.001;
#endif

#if TEST_TYPE == POISEUILLE_TEST
const int N_RUN = 60000;
const int N_VELOCITY_PRINT_INTERVAL = 2000;// N_RUN / 10;
const int N_DENSITY_PRINT_INTERVAL = N_VELOCITY_PRINT_INTERVAL;
const int N_GRID_X_DIM_FLUID = 10;
const int N_GRID_Y_DIM_FLUID = 50;
const field_t F_TAU = 0.90;
const field_t F_BODY_FORCE_X = 1.00e-8;
const field_t F_BODY_FORCE_Y = 0;
const field_t F_TOP_PLATE_VELOCITY = 0;
#endif

const std::string S_BC_BASED = "BCBased/";
const std::string S_SHORT_TESTS = "ShortTests/";
const std::string S_TAU_TESTS = "TauTests/";
const std::string S_LENGTH_TESTS_2 = "LengthTests2/";
const std::string S_DIAMETER_TESTS = "DiameterTests/";

const std::string S_OUTPUT_DIRECTORY_BASE = "SimulationOutput/";

// Populations related
const std::string S_POPULATION_OUTPUT_PATH_BASE = "Populations/";
const std::string S_POPULATION_OUTPUT_FULL_PATH = S_OUTPUT_DIRECTORY_BASE + S_POPULATION_OUTPUT_PATH_BASE;

// Density related
const std::string S_DENSITY_OUTPUT_DIRECTORY_BASE = "DensityTests/";
const std::string S_DENSITY_OUTPUT_FILE_NAME_BASE = "DensityTest";
const std::string S_DENSITY_OUTPUT_FULL_PATH = S_OUTPUT_DIRECTORY_BASE + S_DENSITY_OUTPUT_DIRECTORY_BASE;

// Velocity related
const std::string S_FILE_NAME_BASE = "Velocity";
const std::string S_FLOW_TYPE = "VelocityTests/";
const std::string S_FLOW_SUBCATEGORY = "";
const std::string S_DIRECTORY_VELOCITY_OUTPUT = S_OUTPUT_DIRECTORY_BASE + S_FLOW_TYPE + S_FLOW_SUBCATEGORY;
//const std::string S_DIRECTORY_VELOCITY_OUTPUT = "";