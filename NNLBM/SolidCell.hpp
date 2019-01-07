#pragma once
#include "Globals.hpp"
#include "Cell.hpp"
#include <iostream>
#include <algorithm>


class SolidCell : public Cell
{
private:
	

public:
	
	~SolidCell() = default;
	SolidCell() = default;

	
	//*****************************************************************************************
	// Do functions
		
	//// BounceBack
	//// The current population of this run is reversed and made available for the next run.
	//void collide(const bool runIndex) override {
	//	field_t currentPopulation;		
	//	for (int cellDirection = 0; cellDirection < nDirections / 2; cellDirection++) {
	//		currentPopulation = populations[getArrayIndex(runIndex, cellDirection)];
	//		populations[getArrayIndex(runIndex, cellDirection)] = populations[getArrayIndex(runIndex, reverseDirectionIndex(cellDirection))];
	//		populations[getArrayIndex(runIndex, reverseDirectionIndex(cellDirection))] = currentPopulation;
	//	}
	//}

	//// TODO: Use other relaxation times
	//void collide(const bool runIndex) override {
	//	//std::cout << "collide" << std::endl;
	//	const int dt = 1;
	//	const field_t tau = 10;
	//	computeRho(runIndex);
	//	computeVelocity(runIndex);
	//	computePopulationsEq(runIndex);
	//	for (int cellDirection = 0; cellDirection < nDirections; cellDirection++) {
	//		populations[getArrayIndex(runIndex, cellDirection)]
	//			= populations[getArrayIndex(runIndex, cellDirection)] - dt * (populations[getArrayIndex(runIndex, cellDirection)] - populationsEq[getArrayIndex(runIndex, cellDirection)]) / tau;
	//	}
	//}

	//Half way bounce back
	void setReceived(const bool runIndex, const int populationIndex, const field_t fieldValue) {
		int arrayIndex = getArrayIndex(runIndex, reverseDirectionIndex(populationIndex));
		//std::cout << "setPopulation ARRAY_INDEX : " << arrayIndex << std::endl;;
		assert(("setPopulations: arrayIndex is negative", arrayIndex >= 0));
		assert(("setPopulations: arrayIndex to high", arrayIndex < nFieldDuplicates * nPopulations));
		populations[arrayIndex] = fieldValue;
	}

	// BounceBack and proppagate
	// The current population of this run is reversed and made available in the neighbour cell for the next run.
	void collideAndPropagate(const bool runIndex) override{
		field_t currentPopulation;
		std::shared_ptr<Cell> targetCell;

		for (int cellDirection = 0; cellDirection < nDirections; cellDirection++) {
			currentPopulation = getPolulation(runIndex, cellDirection);
			targetCell = neighbours.getNeighbour(cellDirection);
			targetCell->setPopulation(!runIndex, reverseDirectionIndex(cellDirection), currentPopulation);
		}
	}

	char getCellTypeChar() const override {
		return 'S';
	}
	
	CellType getCellType() const override{
		return CellType::solidCell;
	}
};