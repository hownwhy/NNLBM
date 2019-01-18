#pragma once
#include "Cell.hpp"
#include <iostream>

class BulkCell : public Cell 
{
public:
	~BulkCell() = default;
	BulkCell() = default;


	//*****************************************************************************************
	// Do functions

	//// TODO: Use other relaxation times
	//void collide(const bool runIndex) override{
	//	//std::cout << "collide" << std::endl;
	//	const int dt = 1;
	//	const field_t tau = 10;
	//	computeDensity(runIndex);
	//	computeVelocity(runIndex);
	//	computePopulationsEq();
	//	for (int cellDirection = 0; cellDirection < nDirections; cellDirection++) {			
	//		populations[getArrayIndex(runIndex, cellDirection)]
	//			= populations[getArrayIndex(runIndex, cellDirection)] - dt * (populations[getArrayIndex(runIndex, cellDirection)] - populationsEq[getArrayIndex(runIndex, cellDirection)]) / tau;
	//	}
	//}

	void collide(const bool runIndex) {
		//std::cout << "\n\nDensity before computeDensity = " << density_;
		computeDensity(runIndex);
		//std::cout << "\n\nDensity after computeDensity = " << density_;
		computeVelocity(runIndex, force);
		computePopulationsEq();
		computeForcePopulations();
		int populationIndex = 0;
		for (int cellDirection = 0; cellDirection < nDirections; cellDirection++) {		
			populationIndex = getArrayIndex(runIndex, cellDirection);
			populations_[populationIndex]
				= populations_[populationIndex] - dt * (populations_[populationIndex] - populationsEq_[cellDirection]) / tau;
			// Add force term
			populations_[populationIndex] += ((1.0 - (dt / (2 * tau))) * forcePopulations_[cellDirection]);
		}		
	}

	void setReceived(const bool runIndex, const int populationIndex, const field_t fieldValue) {
		int arrayIndex = getArrayIndex(runIndex, populationIndex);
		populations_[arrayIndex] = fieldValue;
	}

	// The only real difference of collideAndPropagate compared with collide is where the result is stored:
	// TODO: put correct source term in call to "comouteVelocity"
	void collideAndPropagate(const bool runIndex) override{
		const int dt = 1;
		const field_t tau = 1;
		computeDensity(runIndex);
		computeVelocity(runIndex, force);
		computePopulationsEq();
		field_t currentPopulation;
		std::shared_ptr<Cell> targetCell;

		for (int cellDirection = 0; cellDirection < nDirections; cellDirection++) {			
			int arrayIndex = getArrayIndex(runIndex, cellDirection);
			assert(("collideAndPropagate: arrayIndex is negaive:" , arrayIndex >= 0));
			assert(("collideAndPropagate: arrayIndex to high", arrayIndex < nFieldDuplicates * nPopulations));
			// Imediate relaxation f[i] = f_eq[i]
			currentPopulation = populationsEq_[cellDirection];
			targetCell = neighbours_.getNeighbour(cellDirection);			
			targetCell->setPopulation(!runIndex, cellDirection, currentPopulation);

			/*populations[getArrayIndex(!runIndex, cellDirection)]
				= populations[getArrayIndex(runIndex, cellDirection)] - dt * (populations[getArrayIndex(runIndex, cellDirection)] - populationsEq[getArrayIndex(runIndex, cellDirection)]) / tau;*/
		}
	}
	
	char getCellTypeChar() const override{
		return 'B';
	}
	
	CellType getCellType() const override {
		return CellType::bulkCell;
	}
};