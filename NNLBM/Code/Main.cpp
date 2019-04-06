#pragma once
#include "Globals.hpp"
#include "Grid.hpp"
#include "InputOutput.hpp"
#include <sstream>
#include <string>
//#include <exception>
#include <chrono>
using Clock = std::chrono::high_resolution_clock;


// !!!!!*************************!!!!!************************************!!!!!
// The use of hpp files for almost all code in this project is done in order to
// make code development faster. Most code will be moved to cpp files later on.
// Templates will remain in the header files.

//void setSpeedAtPoint(Grid grid, uint_t x, uint_t y) {
//	grid.getCell(x, y)->initializeVelocity(SpatialDirection::x, 0.01);
//	grid.getCell(x, y+1)->initializeVelocity(SpatialDirection::x, -0.01);
//	grid.getCell(x, y+2)->initializeVelocity(SpatialDirection::y, 0.01);
//	grid.getCell(x, y+3)->initializeVelocity(SpatialDirection::y, -0.01);
//
//}

int main() {

	bool runIndex = 0;
	Grid grid;

	std::string populationOutputString = "";
	std::string velocityString = "";
	std::string velocityFileName = "";
	std::string densityString = "";
	std::string densityFileName = "";

#if TEST_TYPE & STREAM_TEST_PIPE + COUETTE_TEST + POISEUILLE_TEST
	grid.makePipeGeometry();
#elif TEST_TYPE & STREAM_TEST_BOX + CAVITY_TEST
	grid.makeBoxGeometry();
#endif

	//if (TEST_TYPE & (STREAM_TEST_PIPE | COUETTE_TEST | POISEUILLE_TEST)) {
	//	grid.makePipeGeometry();
	//}
	//else if (TEST_TYPE & (STREAM_TEST_BOX | CAVITY_TEST)) {
	//	grid.makeBoxGeometry();
	//}

#if 0
	grid.printGeometry();
	system("pause");
#endif

	grid.makeGrid();
#if 0
	grid.printCellType();
	system("pause");
#endif	

	grid.linkNeighbours();
#if 0
	grid.printNeighboursCellType();
	system("pause");
#endif

#if 0
	grid.printCellPopulation(runIndex);
	system("pause");
#endif

#if 0
	grid.propagate(runIndex);
	grid.printCellPopulation(!runIndex);
	system("pause");
#endif

	grid.gridInitialize();
	//grid.gridInitialize(!runIndex);

#if 0
	std::cout << "Initial population for runIndex" << std::endl;
	grid.printCellPopulation(runIndex);
	std::cout << "Initial population for !runIndex" << std::endl;
	grid.printCellPopulation(!runIndex);
	system("pause");
#endif

#if 0
	grid.printCellVelocity(runIndex);
	system("pause");
#endif



#if TEST_TYPE == STREAM_TEST_PIPE | STREAM_TEST_BOX
	//grid.getCell(1, 1)->initializeVelocity(SpatialDirection::x, 0.9);
	//grid.getCell(2, 2)->initializeVelocity(SpatialDirection::x, 0.9);
	//grid.getCell(1, 3)->initializeVelocity(SpatialDirection::x, 0.9);
	//grid.getCell(1, 2)->initializeVelocity(SpatialDirection::y, -0.9);
	//grid.getCell(2, 2)->initializeVelocity(SpatialDirection::y, 0.9);
	//grid.getCell(3, 2)->initializeVelocity(SpatialDirection::y, -0.9);


	//grid.getCell(1, 2)->setPopulation(0, CellDirection::east, 0.8);
	//grid.getCell(1, 2)->setPopulation(1, CellDirection::northEast, 2);
	/*grid.getCell(2, 2)->computeDensity(runIndex);
	grid.getCell(2, 2)->computeVelocity(runIndex);
	grid.getCell(2, 2)->computePopulationsEq();*/

	//grid.getCell(1, 3)->setPopulation(CellDirection::west, 1.0);
	//grid.getCell(2, 2)->setPopulation(CellDirection::east, 1.0);
	//grid.getCell(2, 2)->setPopulation(CellDirection::south, 1.0);
	//grid.getCell(2, 2)->setPopulation(CellDirection::north, 1.0);
	//grid.getCell(2, 2)->setPopulation(CellDirection::west, 1.0);
	//grid.getCell(5, 6)->setPopulation(0, CellDirection::west, 0.9);
#endif


	for (int superRun = 0; superRun < N_TESTRUN_LOOPS; superRun++) {

		const uint_t nRun = N_RUN;
		velocityFileName = getVelocityFileName();
		densityFileName = getDensityFileName();
		uint_t velocityPrintInterval{ N_VELOCITY_PRINT_INTERVAL };
		//uint_t populationPrintInterval{ N_POPULATION_PRINT_INTERVAL };

		auto t1 = Clock::now();

		for (uint_t run = 0; run < nRun; run++) {
			if (run % velocityPrintInterval == 0) {
				std::cout << "\r Processing: " << run << " of " << nRun;
#if !TESTRUN
				grid.appendGridVelocityList(runIndex, velocityString);
				//grid.appendGridDensityList(runIndex, densityString);								
#endif
				velocityPrintInterval = velocityPrintInterval * 2;
			}

			//if (run % populationPrintInterval == 0) {
			//	std::cout << "\r Processing: " << run << " of " << nRun;
			//	grid.appendGridPoplulationsList(runIndex, populationOutputString);
			//	//grid.appendGridNonEqPoplulationsList(runIndex, populationOutputString);
			//}

			grid.collideAndPropagate(runIndex);
			runIndex = !runIndex;
			//system("pause");
		}

		auto t2 = Clock::now();
		auto nodeUpdates = N_RUN * N_GRID_X_DIM_FLUID * N_GRID_Y_DIM_FLUID;
		auto totalTime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
		double mlups = nodeUpdates * 1000. / totalTime;
		std::cout << "\nNode updates = " << nodeUpdates << std::endl;
		std::cout << "Total time = " << totalTime << std::endl;
		std::cout << "MLUPS = " << mlups << std::endl;		
	}

	
	//system("pause");

#if !TESTRUN
	//populationListToFile(populationOutputString, "population.txt");
	velocityListToFile(velocityString, velocityFileName);
	//densityListToFile(densityString, densityFileName);
#endif
	system("pause");
}