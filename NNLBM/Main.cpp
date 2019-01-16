#pragma once
#include "Globals.hpp"
#include <Windows.h>
#include "Grid.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <exception>

// !!!!!*************************!!!!!************************************!!!!!
// The use of hpp files for almost all code in this project is done in order to
// make code development faster. Most code will be moved to cpp files later on.
// Templates will remain in the header files.

void createFolder(std::string path)
{
	if (!CreateDirectory(path.c_str(), NULL))
	{
		return;
	}
}


void stringToFile(const std::string string, const std::string filename){	
	createFolder(S_OUTPUT_DIRECTORY_BASE);
	createFolder(S_OUTPUT_DIRECTORY_BASE + S_FLOW_TYPE);
	createFolder(S_DIRECTORY_VELOCITY_OUTPUT);
	std::ofstream myfile(S_DIRECTORY_VELOCITY_OUTPUT + filename);
	std::string outputString = string;
	outputString.insert(0, "{");
	outputString.append("}");

	if (myfile.is_open())
	{
		myfile << outputString;
		myfile.close();
	}
	else {
		std::cout << "Unable to open file";
		system("pause");
	}
}

void densityListToFile(const std::string string, const std::string filename) {
	createFolder(S_DENSITY_OUTPUT_DIRECTORY_BASE);
	std::ofstream myfile(S_DENSITY_OUTPUT_DIRECTORY_BASE + filename);
	std::string outputString = string;
	outputString.insert(0, "{");
	outputString.append("}");

	if (myfile.is_open())
	{
		myfile << outputString;
		myfile.close();
	}
	else {
		std::cout << "\nUnable to open file: " << S_DENSITY_OUTPUT_FILE_NAME_BASE << filename;
		system("pause");
	}
}


int main() {

	bool runIndex = 0;
	//field_t rho = 1;
	//field_t velocity = 0;
	Grid grid;

	std::string populationOutputString = "";
	//std::string populationFileName = "";
	std::string velocityString = "";
	std::string velocityFileName = "";
	std::string densityString = "";
	std::string densityFileName = "";
	

	grid.makeGeometry();
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

	grid.gridInitialize(runIndex);
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







#if 1
	//// Streaming test
	//const int nRun = 30;
	//const int printInterval = 1;

	//// Collision test
	//const int nRun = 5;
	//const int printInterval = 1;
	
	// Poisuille test
	const int nRun = N_RUN;
	const int printInterval = N_VELOCITY_PRINT_INTERVAL;

	std::ostringstream velocityFileNameStream;
	velocityFileNameStream << std::setprecision(3);
	velocityFileNameStream << S_FILE_NAME_BASE << "_X" << N_GRID_X_DIM
		<< "_Y" << N_GRID_Y_DIM
		<< "_T" << N_RUN
		<< "_Step" << N_VELOCITY_PRINT_INTERVAL
		<< "_Tau" << F_TAU
		<< "_F" << F_BODY_FORCE_X
		<< ".txt";
	velocityFileName = velocityFileNameStream.str();

	std::ostringstream densityFileNameStream;
	densityFileNameStream << std::setprecision(3);
	densityFileNameStream << S_DENSITY_OUTPUT_FILE_NAME_BASE << "_X" << N_GRID_X_DIM
		<< "_Y" << N_GRID_Y_DIM
		<< "_T" << N_RUN
		<< "_Step" << N_DENSITY_PRINT_INTERVAL
		<< "_Tau" << F_TAU
		<< "_F" << F_BODY_FORCE_X
		<< ".txt";
	densityFileName = densityFileNameStream.str();


	//grid.getCell(3, 3)->initializeVelocity(runIndex, SpatialDirection::x, 0.4);
	//grid.getCell(3, 3)->initializeVelocity(runIndex, SpatialDirection::y, -0.4);
	//grid.getCell(3, 3)->initializeVelocity(!runIndex, SpatialDirection::x, 0.4);
	//grid.getCell(3, 3)->initializeVelocity(!runIndex, SpatialDirection::y, -0.4);

	//grid.getCell(4, 1)->initializeVelocity(runIndex, SpatialDirection::x, 0.9);
	//grid.getCell(5, 4)->initializeVelocity(!runIndex, SpatialDirection::x, -0.9);
	//grid.getCell(1, 4)->initializeVelocity(runIndex, SpatialDirection::y, -0.9);
	//grid.getCell(4, 5)->initializeVelocity(!runIndex, SpatialDirection::y, 0.9);
	

	/*std::cout << std::endl;
	std::cout << grid.getCell(3, 3)->getVelocity(runIndex, SpatialDirection::x) << "\n";
	std::cout << grid.getCell(3, 3)->getVelocity(runIndex, SpatialDirection::y) << "\n";
	std::cout << grid.getCell(2, 2)->getVelocity(!runIndex, SpatialDirection::x) << "\n";
	std::cout << grid.getCell(2, 2)->getVelocity(!runIndex, SpatialDirection::y) << "\n";*/

	/*grid.printCellVelocity(runIndex);
	system("pause");*/


	populationOutputString = "";
	velocityString = "";
	densityString = "";
	

	/*grid.appendGridPolulationsList(runIndex, populationOutputString);
	stringToFile(populationOutputString, "population.txt");
	grid.appendGridVelocityList(runIndex, velocityString);
	stringToFile(velocityString, "velocity.txt");*/
	//system("pause");
	//std::cout << "1*10^-8 = " << (1e-8) << std::endl;
	//system("pause");
	for (int run = 0; run < nRun; run++) {
		
		grid.collide(runIndex);
		if (run % printInterval == 0) {
			std::cout << "\r Processing: " << run << " of " << nRun;
			/*grid.appendGridPolulationsList(runIndex, populationOutputString);
			stringToFile(populationOutputString, "population.txt");*/
			/*std::cout << velocityString;
			system("pause");*/
			grid.appendGridVelocityList(runIndex, velocityString);			
			stringToFile(velocityString, velocityFileName);

			grid.appendGridDensityList(runIndex, densityString);
			densityListToFile(densityString, densityFileName);

		}
		grid.propagate(runIndex);
	/*	grid.appendGridPolulationsList(!runIndex, populationOutputString);
		grid.appendGridVelocityList(!runIndex, velocityString);
		stringToFile(populationOutputString, "testfile.txt");
		stringToFile(velocityString, "velocity.txt");*/

		//system("pause");
		
		runIndex = !runIndex;
	}
	/*stringToFile(velocityString, "velocity.txt");*/
	/*system("pause");*/

#endif



#if 0
	std::array<int, 2> arrayy = { 1,2 };
	try {
		arrayy[2] = 0;		
	}
	catch(std::exception& e){
		std::cout << "Catch" << e.what();
		system("pause");
	}
	std::cout << arrayy[2];
	system("pause");
#endif

	
	
}