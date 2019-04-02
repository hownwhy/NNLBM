#pragma once
#include "Cell.hpp"
#include <iostream>
#include <chrono>
using Clock = std::chrono::high_resolution_clock;

class BulkCell : public Cell
{
public:
	~BulkCell() = default;
	BulkCell() = default;


	//*****************************************************************************************
	// Do functions	


	void collideAndPropagate() override {				
		computeDensity();
		computeVelocity();
		computePopulationsEq();
		//computeForcePopulations();

#if DEBUG
		nodeCalcCounter++;
		auto t1 = Clock::now();
#endif

		for (int cellDirection = 0; cellDirection < nCellDirections; cellDirection++) {
			populations_[currentPopulationIndexOffset_ + cellDirection] = (1 - (dt / tau)) * populations_[currentPopulationIndexOffset_ + cellDirection] + (dt / tau) * populationsEq_[cellDirection] // Collision			
				;// +forcePopulations_[cellDirection]; // Adding force term
			
			neighbours_[cellDirection]->setReceived(populations_, (CellDirection)cellDirection);
		}
		// Updating rest population (no propagation)
		populations_[getArrayIndex(!runIndex_, CellDirection::rest)] = (1 - (dt / tau)) * populations_[currentPopulationIndexOffset_ + CellDirection::rest] + (dt / tau) * populationsEq_[CellDirection::rest] // Collision			
			;//+ forcePopulations_[CellDirection::rest]; // Adding force term

#if DEBUG
		auto t2 = Clock::now();
		loopTime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
		accumulatedLoopTime += loopTime;
		//accumulatedLoopTime += loopTime > 1000 * (averageLoopTime = accumulatedLoopTime / nodeCalcCounter) && nodeCalcCounter > 1? averageLoopTime : loopTime ;
		averageLoopTime = accumulatedLoopTime / nodeCalcCounter;
		
		// This is an attemt at counting cach misses by loop times deviating considerably from the average
		if (loopTime > 10 * averageLoopTime) {
			cacheMiss++;			
		}
#endif

	}


	/*void setReceived(std::array<field_t, nFieldDuplicates * nPopulations> &sourcePopulations, const CellDirection cellDirection) override {
		populations_.at(nextPopulationIndexOffset_ + cellDirection) = sourcePopulations.at(currentPopulationIndexOffset_ + cellDirection);
	}*/

	/*void setReceived(std::array<field_t, nFieldDuplicates * nPopulations> &sourcePopulations, const CellDirection &cellDirection) override {
		populations_[nextPopulationIndexOffset_ + cellDirection] = sourcePopulations[currentPopulationIndexOffset_ + cellDirection];
	}*/

	void setReceived(field_t *sourcePopulations, const CellDirection &cellDirection) override {
		populations_[nextPopulationIndexOffset_ + cellDirection] = sourcePopulations[currentPopulationIndexOffset_ + cellDirection];
	}

	char getCellTypeChar() const override {
		return 'B';
	}

	CellType getCellType() const override {
		return CellType::bulk;
	}
};
