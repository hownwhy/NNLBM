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

#define TEST_TYPE 0x08

#if TEST_TYPE == STREAM_TEST_PIPE
const int N_RUN = 5;
const int N_VELOCITY_PRINT_INTERVAL = 1;
const int N_DENSITY_PRINT_INTERVAL = N_VELOCITY_PRINT_INTERVAL;
const int N_GRID_X_DIM_FLUID = 5;
const int N_GRID_Y_DIM_FLUID = 5;
const field_t F_TAU = 100;
const field_t F_BODY_FORCE_X = 0;
const field_t F_BODY_FORCE_Y = 0;
const field_t F_TOP_PLATE_VELOCITY = 0.8;
#endif

#if TEST_TYPE == STREAM_TEST_BOX
const int N_RUN = 30;
const int N_VELOCITY_PRINT_INTERVAL = 1;
const int N_DENSITY_PRINT_INTERVAL = N_VELOCITY_PRINT_INTERVAL;
const int N_GRID_X_DIM_FLUID = 5;
const int N_GRID_Y_DIM_FLUID = 5;
const field_t RE = 1;
const field_t F_LID_VELOCITY = 1;
const field_t F_BODY_FORCE_X = 0.;
const field_t F_BODY_FORCE_Y = 0.;
const field_t F_TAU = 1;
#endif

#if TEST_TYPE == COUETTE_TEST
const int N_RUN = 1<<(17);
const int N_VELOCITY_PRINT_INTERVAL = 2;
const int N_DENSITY_PRINT_INTERVAL = N_VELOCITY_PRINT_INTERVAL;
const int N_GRID_X_DIM_FLUID = 2;
const int N_GRID_Y_DIM_FLUID = 100;
const field_t RE = 60;
const field_t F_LID_VELOCITY = 0.1;
const field_t F_BODY_FORCE_X = 0.;
const field_t F_BODY_FORCE_Y = 0.;
const field_t F_TAU = 0.5 + (3. * F_LID_VELOCITY * N_GRID_Y_DIM_FLUID / RE);
#endif

#if TEST_TYPE == CAVITY_TEST
const int N_RUN						= 1<<17;
const int N_VELOCITY_PRINT_INTERVAL = 2;
const int N_DENSITY_PRINT_INTERVAL	= N_VELOCITY_PRINT_INTERVAL;
const int N_GRID_X_DIM_FLUID		= 100;
const int N_GRID_Y_DIM_FLUID		= 100;
const field_t RE					= 400;
const field_t F_LID_VELOCITY		= 0.1;
const field_t F_BODY_FORCE_X		= 0.;
const field_t F_BODY_FORCE_Y		= 0.;
const field_t F_TAU = 0.5 + (3. * F_LID_VELOCITY * N_GRID_X_DIM_FLUID / RE);

#endif

#if TEST_TYPE == POISEUILLE_TEST
const int N_RUN = 1<<17;
const int N_VELOCITY_PRINT_INTERVAL = 2;
const int N_DENSITY_PRINT_INTERVAL = N_VELOCITY_PRINT_INTERVAL;
const int N_GRID_X_DIM_FLUID = 2;
const int N_GRID_Y_DIM_FLUID = 100;
const field_t RE = 100;
const field_t F_LID_VELOCITY = 0.;
const field_t F_MAX_VELOCITY = 0.1;
const field_t F_BODY_FORCE_X = 4. * F_MAX_VELOCITY * F_MAX_VELOCITY / N_GRID_Y_DIM_FLUID / RE;
const field_t F_BODY_FORCE_Y = 0.;// 4. * F_LID_VELOCITY * F_LID_VELOCITY / N_GRID_Y_DIM_FLUID / RE;
const field_t F_TAU = 0.5 + (3. * F_MAX_VELOCITY * N_GRID_Y_DIM_FLUID / RE);
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
const std::string S_FLOW_SUBCATEGORY = "Cavity Simulations For Poster/";
const std::string S_DIRECTORY_VELOCITY_OUTPUT = S_OUTPUT_DIRECTORY_BASE + S_FLOW_TYPE + S_FLOW_SUBCATEGORY;
//const std::string S_DIRECTORY_VELOCITY_OUTPUT = "";

