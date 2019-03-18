#pragma once
#include "Globals.hpp"
#include <Windows.h>
#include <iostream>
#include <fstream>
#include <sstream>
//#include <iomanip>
#include <string>
//#include <exception>

void createFolder(std::string path)
{
	if (!CreateDirectory(path.c_str(), NULL))
	{
		return;
	}
}

void stringListToFile(const std::string string, const std::string path, const std::string filename) {
	std::ofstream myfile(path + filename);
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

void populationListToFile(const std::string string, const std::string filename) {
	createFolder(S_OUTPUT_DIRECTORY_BASE);
	createFolder(S_OUTPUT_DIRECTORY_BASE + S_POPULATION_OUTPUT_PATH_BASE);
	createFolder(S_POPULATION_OUTPUT_FULL_PATH);

	stringListToFile(string, S_POPULATION_OUTPUT_FULL_PATH, filename);
}


void velocityListToFile(const std::string string, const std::string filename) {
	createFolder(S_OUTPUT_DIRECTORY_BASE);
	createFolder(S_OUTPUT_DIRECTORY_BASE + S_FLOW_TYPE);
	createFolder(S_DIRECTORY_VELOCITY_OUTPUT);

	stringListToFile(string, S_DIRECTORY_VELOCITY_OUTPUT, filename);
}

void densityListToFile(const std::string string, const std::string filename) {
	createFolder(S_OUTPUT_DIRECTORY_BASE);
	createFolder(S_DENSITY_OUTPUT_FULL_PATH);
	std::ofstream myfile(S_DENSITY_OUTPUT_FULL_PATH + filename);
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

std::string getVelocityFileName() {
	std::ostringstream velocityFileNameStream;
	velocityFileNameStream << std::setprecision(4);
	velocityFileNameStream << S_FILE_NAME_BASE 
		<< "_X" << N_GRID_X_DIM_FLUID
		<< "_Y" << N_GRID_Y_DIM_FLUID
		<< "_T" << N_RUN
		<< "_Step" << N_VELOCITY_PRINT_INTERVAL
		<< "_Tau" << F_TAU
		<< "_Re" << RE
		<< "_LidV" << F_LID_VELOCITY
		<< "_F" << F_BODY_FORCE_X
		<< ".txt";
	return velocityFileNameStream.str();
}

std::string getDensityFileName() {
	std::ostringstream densityFileNameStream;
	densityFileNameStream << std::setprecision(3);
	densityFileNameStream << S_DENSITY_OUTPUT_FILE_NAME_BASE << "_X" << N_GRID_X_DIM_FLUID
		<< "_Y" << N_GRID_Y_DIM_FLUID
		<< "_T" << N_RUN
		<< "_Step" << N_DENSITY_PRINT_INTERVAL
		<< "_Tau" << F_TAU
		<< "_F" << F_BODY_FORCE_X
		<< ".txt";
	return densityFileNameStream.str();

}