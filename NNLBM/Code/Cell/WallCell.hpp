#pragma once
#include "../Globals.hpp"
#include "Cell.hpp"
#include <iostream>
#include <algorithm>


class WallCell : public Cell
{
private:


public:

	~WallCell() = default;
	WallCell() = default;


	//*****************************************************************************************
	// Do functions

	// Propagate populations to neighbouring cells (notice the: !runIndex) used
	// so that it will not overwrite populations that will be used in the current run (loop).	
	// TODO: The local variables used for readability may (or may not) decrease performance. Find out. 
	void propageteTo(const bool runIndex) override {
		field_t propagationPopulation;
		std::shared_ptr<Cell> targetCell;

		for (int cellDirection = 0; cellDirection < nDirections; cellDirection++) {
			targetCell = neighbours_.getNeighbour(cellDirection);
			int populationIndex = getArrayIndex(runIndex, cellDirection);
			switch (targetCell->getCellType()) {
			case CellType::bulk:
			case CellType::wall:				
				targetCell->setReceived(!runIndex, cellDirection, populations_[populationIndex]);
				break;
			case CellType::solid:
				// Add velocity from moving walls
				const field_t csSqr = 1. / 3;
				// ! Excludes rest population in the following definitions of c and w. !
				const std::array<int, nDirections * nDimensions> c = { 1,1,0,-1,-1,-1,0,1, 0,1,1,1,0,-1,-1,-1 };
				const std::array<int, nDirections * nDimensions> cMask = { 0,1,0,1,0,0,0,0, 0,0,0,0,0,0,0,0 };
				const std::array<field_t, nDirections> w = { (1. / 9), (1. / 36), (1. / 9), (1. / 36), (1. / 9), (1. / 36), (1. / 9), (1. / 36) };				
				int index = (SpatialDirection::x * nDirections) + cellDirection;
				populations_[populationIndex] -= (2 * w[cellDirection] * density_ / csSqr)
					* ((c[index] * cMask[index] * wallVelocityX)
						+ (c[index] * cMask[index] * wallVelocityY));
				
				// Bounce back
				setReceived(!runIndex, reverseDirectionIndex(cellDirection), populations_[populationIndex]);
				break;
			}
		}
	}

	// BounceBack and proppagate
	// The current population of this run is reversed and made available in the neighbour cell for the next run.
	void collideAndPropagate(const bool runIndex) override {
		computeDensity(runIndex);
		computeVelocity(runIndex);
		computePopulationsEq();
		computeForcePopulations();
		int populationIndex = 0;
		field_t propagationPopulation = 0;
		std::shared_ptr<Cell> targetCell;

		for (int cellDirection = 0; cellDirection < nDirections; cellDirection++) {
			populationIndex = getArrayIndex(runIndex, cellDirection);
			targetCell = neighbours_.getNeighbour(cellDirection);
			propagationPopulation = populations_[populationIndex] - dt * (populations_[populationIndex] - populationsEq_[cellDirection]) / tau;
			// Add force term
			propagationPopulation += ((1.0 - (dt / (2 * tau))) * forcePopulations_[cellDirection]);

			switch (targetCell->getCellType()) {
			case CellType::bulk:
			case CellType::wall:
				targetCell->setReceived(!runIndex, cellDirection, propagationPopulation);
				break;
			case CellType::solid:
				setReceived(!runIndex, reverseDirectionIndex(cellDirection), propagationPopulation);
				break;
			}
		}
	}

	////Half way bounce back
	//void setReceived(const bool runIndex, const int populationIndex, const field_t fieldValue) {
	//	int arrayIndex = getArrayIndex(runIndex, reverseDirectionIndex(populationIndex));
	//	populations_[arrayIndex] = fieldValue;
	//}

	// BounceBack and proppagate
	// The current population of this run is reversed and made available in the neighbour cell for the next run.
	/*void collideAndPropagate(const bool runIndex) override {
		field_t currentPopulation;
		std::shared_ptr<Cell> targetCell;

		for (int cellDirection = 0; cellDirection < nDirections; cellDirection++) {
			currentPopulation = getPolulation(runIndex, cellDirection);
			targetCell = neighbours_.getNeighbour(cellDirection);
			targetCell->setPopulation(!runIndex, reverseDirectionIndex(cellDirection), currentPopulation);
		}
	}*/

	char getCellTypeChar() const override {
		return 'W';
	}

	CellType getCellType() const override {
		return CellType::wall;
	}
};