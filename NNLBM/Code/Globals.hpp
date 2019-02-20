#pragma once
#include <assert.h>
#include <string>

typedef double field_t;

static const int nDirections = 8;

//#define N_SOLID_CELL_MARGIN_X 1
//#define N_SOLID_CELL_MARGIN_Y 1



//#if BOOK_BOUNCE_BACK == 1
//const int N_RUN = 10;
//#endif
//
//#if BOOK_BOUNCE_BACK == 0
//const int N_RUN = 11;
//#endif

#define B_CONSOL_OUT_DENSITY  0;


#define BOOK_BOUNCE_BACK 1

#define STREAM_TEST_PIPE     0x01
#define STREAM_TEST_BOX		 0x02
#define COUETTE_TEST		 0x04
#define CAVITY_TEST			 0x08
#define POISEUILLE_TEST		 0x10

#define TEST_TYPE 0x10

#if TEST_TYPE == STREAM_TEST_PIPE
const int N_RUN = 10;
const int N_VELOCITY_PRINT_INTERVAL = 1;
const int N_DENSITY_PRINT_INTERVAL = N_VELOCITY_PRINT_INTERVAL;
const int N_GRID_X_DIM_FLUID = 3;
const int N_GRID_Y_DIM_FLUID = 3;
const field_t F_TAU = 10;
const field_t F_BODY_FORCE_X = 0;
const field_t F_BODY_FORCE_Y = 0;
const field_t F_TOP_PLATE_VELOCITY = 0.000;
#endif

#if TEST_TYPE == STREAM_TEST_BOX
const int N_RUN = 20;
const int N_VELOCITY_PRINT_INTERVAL = 1;
const int N_DENSITY_PRINT_INTERVAL = N_VELOCITY_PRINT_INTERVAL;
const int N_GRID_X_DIM_FLUID = 5;
const int N_GRID_Y_DIM_FLUID = 5;
const field_t F_TAU = 1.00;
const field_t F_BODY_FORCE_X = 0;
const field_t F_BODY_FORCE_Y = 0;
const field_t F_TOP_PLATE_VELOCITY = 1;
#endif

#if TEST_TYPE == COUETTE_TEST
const int N_RUN = 150000;
const int N_VELOCITY_PRINT_INTERVAL = 15000;
const int N_DENSITY_PRINT_INTERVAL = N_VELOCITY_PRINT_INTERVAL;
const int N_GRID_X_DIM_FLUID = 2;
const int N_GRID_Y_DIM_FLUID = 129;
const field_t F_TAU = 0.6161;
//const field_t F_BODY_FORCE_X = 1.00e-8;
const field_t F_BODY_FORCE_X = 0;
const field_t F_BODY_FORCE_Y = 0;
const field_t F_TOP_PLATE_VELOCITY = 0.003;
#endif

#if TEST_TYPE == CAVITY_TEST
const int N_RUN = 100000;
const int N_VELOCITY_PRINT_INTERVAL = 5000;
const int N_DENSITY_PRINT_INTERVAL = N_VELOCITY_PRINT_INTERVAL;
const int N_GRID_X_DIM_FLUID = 129;
const int N_GRID_Y_DIM_FLUID = 129;
const field_t F_TAU = 0.6161;
//const field_t F_BODY_FORCE_X = 1.00e-8;
const field_t F_BODY_FORCE_X = 0;
const field_t F_BODY_FORCE_Y = 0;
const field_t F_TOP_PLATE_VELOCITY = 0.03;
#endif

#if TEST_TYPE == POISEUILLE_TEST
const int N_RUN = 100000;
const int N_VELOCITY_PRINT_INTERVAL = 10000;// N_RUN / 10;
const int N_DENSITY_PRINT_INTERVAL = N_VELOCITY_PRINT_INTERVAL;
const int N_GRID_X_DIM_FLUID = 2;
const int N_GRID_Y_DIM_FLUID = 58;
const field_t F_TAU = 0.554;
const field_t F_BODY_FORCE_X = 0.000001284185;
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

