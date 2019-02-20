#pragma once
#include "Globals.hpp"
#include "Grid.hpp"
#include "InputOutput.hpp"
#include <sstream>
#include <string>
//#include <exception>

// !!!!!*************************!!!!!************************************!!!!!
// The use of hpp files for almost all code in this project is done in order to
// make code development faster. Most code will be moved to cpp files later on.
// Templates will remain in the header files.

void setSpeedAtPoint(Grid grid, int x, int y) {
	grid.getCell(x, y)->initializeVelocity(SpatialDirection::x, 0.01);
	grid.getCell(x, y+1)->initializeVelocity(SpatialDirection::x, -0.01);
	grid.getCell(x, y+2)->initializeVelocity(SpatialDirection::y, 0.01);
	grid.getCell(x, y+3)->initializeVelocity(SpatialDirection::y, -0.01);

}

int main() {

	bool runIndex = 0;	
	Grid grid;

	std::string populationOutputString = "";
	std::string velocityString = "";
	std::string velocityFileName = "";
	std::string densityString = "";
	std::string densityFileName = "";
	
	
	if (TEST_TYPE & (STREAM_TEST_PIPE | COUETTE_TEST | POISEUILLE_TEST)) {
		grid.makePipeGeometry();
	}
	else if (TEST_TYPE & (STREAM_TEST_BOX | CAVITY_TEST)) {
		grid.makeBoxGeometry();
	}

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

#if B_CONSOL_OUT_DENSITY
	std::cout << "Initial density runIndex = " << runIndex << std::endl;
	grid.printCellRho(runIndex);
	std::cout << "Initial density !runIndex = " << !runIndex << std::endl;
	grid.printCellRho(!runIndex);
	system("pause");
#endif




#if 1	
	const int nRun = N_RUN;
	const int printInterval = N_VELOCITY_PRINT_INTERVAL;

	velocityFileName = getVelocityFileName();
	densityFileName = getDensityFileName();

#if TEST_TYPE == STREAM_TEST_PIPE
	//grid.getCell(1, 1)->initializeVelocity(SpatialDirection::x, 0.9);
	//grid.getCell(2, 2)->initializeVelocity(SpatialDirection::x, 0.9);
	//grid.getCell(1, 3)->initializeVelocity(SpatialDirection::x, 0.9);
	//grid.getCell(1, 2)->initializeVelocity(SpatialDirection::y, -0.9);
	//grid.getCell(2, 2)->initializeVelocity(SpatialDirection::y, 0.9);
	//grid.getCell(3, 2)->initializeVelocity(SpatialDirection::y, -0.9);
	
	
		grid.getCell(1, 2)->setPopulation(0, CellDirection::northEast, 0.8);
		grid.getCell(1, 2)->setPopulation(1, CellDirection::northEast, 0.8);
		grid.getCell(2, 2)->computeDensity(runIndex);
		grid.getCell(2, 2)->computeVelocity(runIndex);
		grid.getCell(2, 2)->computePopulationsEq();

	//grid.getCell(0, 1)->setPopulation(0, CellDirection::west, 0.9);
	//grid.getCell(2, 2)->setPopulation(0, CellDirection::east, 1.0);
	//grid.getCell(5, 6)->setPopulation(0, CellDirection::west, 0.9);
#endif
	
	//setSpeedAtPoint(grid, 2, 4);


	/*grid.printCellVelocity(runIndex);
	system("pause");*/	

	for (int run = 0; run < nRun; run++) {
#if B_CONSOL_OUT_DENSITY
		std::cout << "\n\nDensity after collide" << std::endl;
		grid.printCellRho(runIndex);
		std::cout << std::endl;
		grid.printCellRho(!runIndex);
		std::cout << std::endl;
		//system("pause");
#endif		
		if (run % printInterval == 0) {
			std::cout << "\r Processing: " << run << " of " << nRun;			
			//grid.appendGridPoplulationsList(runIndex, populationOutputString);

			grid.appendGridNonEqPoplulationsList(runIndex, populationOutputString);

			grid.appendGridVelocityList(runIndex, velocityString);	
			grid.appendGridDensityList(runIndex, densityString);
		}
		//grid.propagate(runIndex);
		grid.collide(runIndex);		
		//grid.moveBoundary(runIndex);
		grid.propagate(runIndex);
		

		//grid.collideAndPropagate(runIndex);

#if B_CONSOL_OUT_DENSITY
		std::cout << "\n\nDensity after propagate" << std::endl;
		grid.printCellRho(runIndex);
		std::cout << std::endl;
		grid.printCellRho(!runIndex);
		std::cout << std::endl;
		system("pause");
#endif
		runIndex = !runIndex;
	}
	populationListToFile(populationOutputString, "population.txt");
	velocityListToFile(velocityString, velocityFileName);
	densityListToFile(densityString, densityFileName);

	/*stringToFile(velocityString, "velocity.txt");*/
	/*system("pause");*/

#endif

//#if 0
//	std::array<int, 2> arrayy = { 1,2 };
//	try {
//		arrayy.at(2) = 0;		
//	}
//	catch(std::exception& e){
//		std::cout << "Catch" << e.what();
//		system("pause");
//	}
//	std::cout << arrayy.at(2);
//	system("pause");
//#endif
	
}