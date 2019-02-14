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

	//void addMovingBoundaryTerm(const bool runIndex) override{
	//	std::cout << "\nBulkCell::addMovingBoundaryTerm()";
	//	system("pause");
	//}
	   
	void propageteTo(const bool runIndex) override{
		field_t propagationPopulation;
		std::shared_ptr<Cell> targetCell;

		for (int direction = 0; direction < nDirections; direction++) {
			propagationPopulation = populations_.at(getArrayIndex(runIndex, direction));
			targetCell = neighbours_.getNeighbour(direction);
			targetCell->setReceived(!runIndex, direction, propagationPopulation);
		}
	}

	void collideAndPropagate(const bool runIndex) override {
		//std::cout << "\n\nDensity before computeDensity = " << density_;
		computeDensity(runIndex);
		//std::cout << "\n\nDensity after computeDensity = " << density_;
		computeVelocity(runIndex);
		computePopulationsEq();
		computeForcePopulations();
		int populationIndex = 0;
		int populationIndex2 = 0;
		field_t propagationPopulation = 0;
		std::shared_ptr<Cell> targetCell;
		for (int cellDirection = 0; cellDirection < nDirections; cellDirection++) {
			populationIndex = getArrayIndex(runIndex, cellDirection);
			populationIndex2 = getArrayIndex(!runIndex, cellDirection);
			targetCell = neighbours_.getNeighbour(cellDirection);
			propagationPopulation = populations_.at(populationIndex) - dt * (populations_.at(populationIndex) - populationsEq_.at(cellDirection)) / tau;
			// Add force term
			propagationPopulation += ((1.0 - (dt / (2 * tau))) * forcePopulations_.at(cellDirection));

			targetCell->setReceived(!runIndex, cellDirection, propagationPopulation);
		}
	}

	char getCellTypeChar() const override {
		return 'B';
	}

	CellType getCellType() const override {
		return CellType::bulk;
	}
};