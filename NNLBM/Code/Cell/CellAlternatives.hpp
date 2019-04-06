#pragma once
#include "../Globals.hpp"
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <array>
#include <memory>
//#include <intrin.h>

#if TESTRUN
#include <chrono>
using Clock = std::chrono::high_resolution_clock;
#endif

const enum CellType : uint_t
{
	solid = 0x00,
	cell = 0x01,
	bulk = 0x02,
	wall = 0x04
};

// The "rest" direction has been given the index 8, as opposed to the more usual 0.
// This was done for compability with the directions indexes starting at east = 0.
// This might not be necessary.
// TODO: Maybe find a way to use the more standard "rest = 0"? 
const enum  CellDirection : uint_t
{
	east = 0,
	northEast = 1,
	north = 2,
	northWest = 3,
	west = 4,
	southWest = 5,
	south = 6,
	southEast = 7,
	rest = 8
};

const enum SpatialDirection : uint_t
{
	x = 0,
	y = 1
};


	// TODO: Where should these static constants be? Here it will be accessible to all who include the hpp file.
	// If declareedd as a member, it has to be initialized in a cpp file...
	static const uint_t nCellDirections = 8;
	static const uint_t nPopulations = 9;		// Number of poulations for species
	//static const uint_t nFieldDuplicates = 2;	// Number of field "duplicats" (for temporary storage)(probably no more than 2)
	static const uint_t nDimensions = 2;
	static const field_t dt = 1.;

	// poiseuille flow
	static const field_t tau = F_TAU;
	static const field_t tauInverse = 1. / tau;
	static const field_t dtOverTau = dt / tau;
	static const field_t oneMinusDtOverTau = 1. - dtOverTau;
	static const std::array<field_t, 2> bodyForce = { F_BODY_FORCE_X, F_BODY_FORCE_Y };
	static const std::array<uint_t, nPopulations> reverseDirectionIndex = { 4,5,6,7,0,1,2,3 };



// TODO: See what can be done with templates instead of virtual functions
class Cell {

private:
	//Neighbours neighbours_;
	//std::array<std::shared_ptr<Cell>, nCellDirections> neighbours_;
	std::array<Cell*, nCellDirections> neighbours_;
	//std::array<field_t, nFieldDuplicates * nPopulations> populations_;
	field_t populations_[2 * nPopulations]; // An array holding two sets of populations
	std::array<field_t, nPopulations> populationsEq_;
	//field_t populationsEq_[nPopulations];
	//std::array<field_t, nPopulations> forcePopulations_;
	std::array<field_t, nDimensions> velocity_;
	field_t density_;
	CellType cellType = CellType::solid;
	//field_t tau;
	
	//std::array<field_t, nPopulations> movingPlateTerm;
	static bool runIndex_;
#if PRECOMPUTED
	static uint_t currentPopulationIndexOffset_;
	static uint_t nextPopulationIndexOffset_;
	static std::array<uint_t, nPopulations> currentPopulationIndexes_;
	static std::array<uint_t, nPopulations> nextPopulationIndexes_;
#elif EXPLICIT
	
#else
	static uint_t currentPopulationIndexOffset_;
	static uint_t nextPopulationIndexOffset_;
#endif
	

	


	field_t applyMovingSolid(field_t &population, const int cellDirection) const {
		// Add velocity from moving walls
		const field_t csSqrInvers = 3.;
		// ! Excludes rest population in the following definitions of c and w. !
		const std::array<int, nCellDirections * nDimensions> c = { 1,1,0,-1,-1,-1,0,1, 0,1,1,1,0,-1,-1,-1 };
		//const std::array<field_t, nCellDirections * nDimensions> c = { 1,1,0,-1,-1,-1,0,1, 0,1,1,1,0,-1,-1,-1 };
		const std::array<field_t, nCellDirections> w = { (1. / 9), (1. / 36), (1. / 9), (1. / 36), (1. / 9), (1. / 36), (1. / 9), (1. / 36) };

		//const std::array<field_t, nCellDirections> factors{ (6 / 9), (6 / 36), (6 / 9), (6 / 36), (6 / 9), (6 / 36), (6 / 9), (6 / 36) };
		return population - (2. * w.at(cellDirection) * density_ * csSqrInvers)
		//return population - (factors.at(cellDirection) * density_)
			* ((c.at((SpatialDirection::x * nCellDirections) + cellDirection) * velocity_[SpatialDirection::x])
				+ (c.at((SpatialDirection::y * nCellDirections) + cellDirection) * velocity_[SpatialDirection::y]));
	}


public:
	Cell() {	
		//tau = F_TAU;
	}

	~Cell() = default;
	
	
	//bool solidQ = true;

	// Returns the 1D array index which depend on the 2D (runIndex, direction).
	inline uint_t getArrayIndex(const bool runIndex, const uint_t direction) const {
		return (runIndex * nPopulations) + direction;
	}

	/* This function gives the opposite direction of what you put in.
	 This is also the reason why I chose to not follow the convension
	 having the rest direction be the 0 direction. 	*/
	/*const uint_t threeHalfDirection = (3 * nCellDirections) / 2;
	inline uint_t reverseDirectionIndex(const uint_t direction) const {
		return ((direction + threeHalfDirection) % nCellDirections);
	}*/

	//TODO: const
	void setNeighbour(Cell &cell, const CellDirection cellDirection) {
		neighbours_.at(cellDirection) = &cell;
	}

	/*std::shared_ptr<Cell> getNeighbour(const CellDirection cellDirection) const {
		return neighbours_.at(cellDirection);
	}*/

	void setSolid() {
		cellType = CellType::solid;
	}
	void setBulk() {
		cellType = CellType::bulk;
	}

	//*****************************************************************************************
	// Cell initialization

	//void initialize(const field_t density, const field_t xVelocity, const field_t yVelocity) {
	//	density_ = density;
	//	velocity_.at(SpatialDirection::x) = xVelocity;
	//	velocity_.at(SpatialDirection::y) = yVelocity;
	//	//std::cout << "density = " << density_ << "\txVelocity = " << velocity_.at(SpatialDirection::x) << "\tyVelocity = " << velocity_.at(SpatialDirection::y) << std::endl;
	//	computePopulationsEq();
	//	std::copy(populationsEq_.begin(), populationsEq_.end(), populations_.begin());
	//	std::copy(populationsEq_.begin(), populationsEq_.end(), populations_.begin() + nPopulations);
	//}

	void initialize(const field_t density, const field_t xVelocity, const field_t yVelocity) {
		density_ = density;
		velocity_.at(SpatialDirection::x) = xVelocity;
		velocity_.at(SpatialDirection::y) = yVelocity;
		//std::cout << "density = " << density_ << "\txVelocity = " << velocity_.at(SpatialDirection::x) << "\tyVelocity = " << velocity_.at(SpatialDirection::y) << std::endl;
		computePopulationsEq();
		for (uint_t populationIndex = 0; populationIndex < nPopulations; populationIndex++) {
			populations_[populationIndex] = populationsEq_[populationIndex];
			populations_[populationIndex + nPopulations] = populationsEq_[populationIndex];
		}	
	}

	//void initializeDensity(const field_t density) {
	//	density_ = density;
	//	computePopulationsEq();
	//	std::copy(populationsEq_.begin(), populationsEq_.end(), populations_.begin());
	//}

	//void initializeVelocity(const SpatialDirection direction, const field_t velocity) {
	//	velocity_.at(direction) = velocity; // -((bodyForce.at(direction) * dt) / (density_));
	//	computePopulationsEq();
	//	std::copy(populationsEq_.begin(), populationsEq_.end(), populations_.begin());
	//}

	//*****************************************************************************************
	// Do functions	

	static void setRunStep(const bool runIndex) {
		runIndex_ = runIndex;

#if PRECOMPUTED
		currentPopulationIndexOffset_ = nPopulations * runIndex;
		nextPopulationIndexOffset_ = nPopulations * !runIndex;
		currentPopulationIndexes_ = {0 + currentPopulationIndexOffset_,	1 + currentPopulationIndexOffset_, 2 + currentPopulationIndexOffset_,
			3 + currentPopulationIndexOffset_, 4 + currentPopulationIndexOffset_, 5 + currentPopulationIndexOffset_,
			6 + currentPopulationIndexOffset_,	7 + currentPopulationIndexOffset_,	8 + currentPopulationIndexOffset_ };
		nextPopulationIndexes_ = { 0 + nextPopulationIndexOffset_,	1 + nextPopulationIndexOffset_, 2 + nextPopulationIndexOffset_,
			3 + nextPopulationIndexOffset_, 4 + nextPopulationIndexOffset_, 5 + nextPopulationIndexOffset_,
			6 + nextPopulationIndexOffset_,	7 + nextPopulationIndexOffset_,	8 + nextPopulationIndexOffset_ };
#elif EXPLICIT
#else
		currentPopulationIndexOffset_ = nPopulations * runIndex;
		nextPopulationIndexOffset_ = nPopulations * !runIndex;
#endif	

	}		
	   	  
	void collideAndPropagate() {
		if (cellType == CellType::bulk) {			
			collideAndPropagateBulk();
		}
		else {
			collideAndPropagateSolid();
		}
	}

	void collideAndPropagateBulk() {				
		computeDensity();
		computeVelocity();		
		computePopulationsEq();		

		//computeForcePopulations();		
#if TESTRUN
		nodeCalcCounter++;
		auto t1 = Clock::now();
#endif

#if !PRECOMPUTED && !EXPLICIT
	uint_t currentPopulationIndexOffset = nPopulations * runIndex_;
#endif

		for (int cellDirection = 0; cellDirection < nCellDirections; cellDirection++) {
			//populations_[currentPopulationIndexOffset_ + cellDirection] = (1 - (dt / tau)) * populations_[currentPopulationIndexOffset_ + cellDirection] + (dt / tau) * populationsEq_[cellDirection] // Collision			
			//	;// +forcePopulations_[cellDirection]; // Adding force term
#if PRECOMPUTED
			populations_[currentPopulationIndexes_[cellDirection]] = oneMinusDtOverTau * populations_[currentPopulationIndexes_[cellDirection]] + dtOverTau * populationsEq_[cellDirection] // Collision			
				;// +forcePopulations_[cellDirection]; // Adding force term
#elif EXPLICIT
			populations_[getArrayIndex(runIndex_, cellDirection)] = oneMinusDtOverTau * populations_[getArrayIndex(runIndex_, cellDirection)] + dtOverTau * populationsEq_[cellDirection] // Collision			
				;// +forcePopulations_[cellDirection]; // Adding force term
#else			
			populations_[currentPopulationIndexOffset + cellDirection] = oneMinusDtOverTau * populations_[currentPopulationIndexOffset + cellDirection] + dtOverTau * populationsEq_[cellDirection] // Collision			
				;// +forcePopulations_[cellDirection]; // Adding force term
#endif
			if (neighbours_[cellDirection]->cellType == CellType::bulk) {
				neighbours_[cellDirection]->setReceivedBulk(populations_, cellDirection);
			}
			else {				
				neighbours_[cellDirection]->setReceivedSolid(populations_, cellDirection);
			}

			//neighbours_[cellDirection]->setReceived(populations_, cellDirection);
		}
		// Updating rest population (no propagation)
		//populations_[nextPopulationIndexOffset_ + CellDirection::rest] = (1 - (dt / tau)) * populations_[currentPopulationIndexOffset_ + CellDirection::rest] + (dt / tau) * populationsEq_[CellDirection::rest] // Collision			
		//	;//+ forcePopulations_[CellDirection::rest]; // Adding force term

#if PRECOMPUTED
		populations_[nextPopulationIndexes_[CellDirection::rest]] = oneMinusDtOverTau * populations_[currentPopulationIndexes_[CellDirection::rest]] + dtOverTau * populationsEq_[CellDirection::rest] // Collision			
			;//+ forcePopulations_[CellDirection::rest]; // Adding force term
#elif EXPLICIT
		populations_[(nPopulations * !runIndex_) + CellDirection::rest] = oneMinusDtOverTau * populations_[getArrayIndex(runIndex_, CellDirection::rest)] + dtOverTau * populationsEq_[CellDirection::rest] // Collision			
			;//+ forcePopulations_[CellDirection::rest]; // Adding force term
#else
		uint_t nextPopulationIndexOffset = nPopulations * !runIndex_;
		populations_[nextPopulationIndexOffset + CellDirection::rest] = oneMinusDtOverTau * populations_[currentPopulationIndexOffset + CellDirection::rest] + dtOverTau * populationsEq_[CellDirection::rest] // Collision			
			;//+ forcePopulations_[CellDirection::rest]; // Adding force term
#endif
#if TESTRUN
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

	void collideAndPropagateSolid() {
		/*std::cout << "\nSolidCell::collideAndProppagate()";
		system("pause");*/
	}

	void setReceivedBulk(const field_t *sourcePopulations, const uint_t cellDirection) {
#if PRECOMPUTED		
		populations_[nextPopulationIndexes_[cellDirection]] = sourcePopulations[currentPopulationIndexes_[cellDirection]];
#elif EXPLICIT
		populations_[(nPopulations * !runIndex_) + cellDirection] = sourcePopulations[getArrayIndex(runIndex_, cellDirection)];
#else
		uint_t currentPopulationIndexOffset = nPopulations * runIndex_;
		uint_t nextPopulationIndexOffset = nPopulations * !runIndex_;
		populations_[nextPopulationIndexOffset + cellDirection] = sourcePopulations[currentPopulationIndexOffset + cellDirection];
#endif
	}

	void setReceivedSolid(field_t *sourcePopulations, const uint_t cellDirection) {
#if PRECOMPUTED	
		sourcePopulations[nextPopulationIndexes_[reverseDirectionIndex[cellDirection]]]
			= applyMovingSolid(sourcePopulations[currentPopulationIndexes_[cellDirection]], cellDirection);
#elif EXPLICIT
		sourcePopulations[(nPopulations * !runIndex_) + reverseDirectionIndex[cellDirection]]
			= applyMovingSolid(sourcePopulations[getArrayIndex(runIndex_, cellDirection)], cellDirection);
#else
		uint_t currentPopulationIndexOffset = nPopulations * runIndex_;
		uint_t nextPopulationIndexOffset = nPopulations * !runIndex_;
		sourcePopulations[nextPopulationIndexOffset + reverseDirectionIndex[cellDirection]]
			= applyMovingSolid(sourcePopulations[currentPopulationIndexOffset + cellDirection], cellDirection);
#endif	
	
	}

	void setReceived(field_t *sourcePopulations, const uint_t cellDirection) {
		//if (cellType == CellType::bulk) {
		//	//populations_[nextPopulationIndexOffset_ + cellDirection] = sourcePopulations[currentPopulationIndexOffset_ + cellDirection];
		//	populations_[nextPopulationIndexes_[cellDirection]] = sourcePopulations[currentPopulationIndexes_[cellDirection]];
		//}
		//else {
		//	/*sourcePopulations[nextPopulationIndexOffset_ + reverseDirectionIndex[cellDirection]]
		//		= applyMovingSolid(sourcePopulations[currentPopulationIndexOffset_ + cellDirection], cellDirection);*/
		//	sourcePopulations[nextPopulationIndexes_[reverseDirectionIndex[cellDirection]]]
		//		= applyMovingSolid(sourcePopulations[currentPopulationIndexes_[cellDirection]], cellDirection);
		//}
	}
	
	//void setReceived(std::array<field_t, nFieldDuplicates * nPopulations> &sourcePopulations, const uint_t &cellDirection) = 0;
	//virtual void setReceived(std::shared_ptr<Cell> sourceCell, const CellDirection cellDirection) = 0;
	//virtual void setReceived(field_t *sourcePopulations, const CellDirection &cellDirection) = 0;
	

	inline void computeDensity() {
		
		//const uint_t nextPopulationIndexOffset_ = nPopulations * !runIndex_;
		/*density_ =
			populations_.at(currentPopulationIndexOffset_ + CellDirection::east)
			+ populations_.at(currentPopulationIndexOffset_ + CellDirection::northEast)
			+ populations_.at(currentPopulationIndexOffset_ + CellDirection::north)
			+ populations_.at(currentPopulationIndexOffset_ + CellDirection::northWest)
			+ populations_.at(currentPopulationIndexOffset_ + CellDirection::west)
			+ populations_.at(currentPopulationIndexOffset_ + CellDirection::southWest)
			+ populations_.at(currentPopulationIndexOffset_ + CellDirection::south)
			+ populations_.at(currentPopulationIndexOffset_ + CellDirection::southEast)
			+ populations_.at(currentPopulationIndexOffset_ + CellDirection::rest);*/
#if PRECOMPUTED
		density_ =
			populations_[currentPopulationIndexes_[CellDirection::east]]
			+ populations_[currentPopulationIndexes_[CellDirection::northEast]]
			+ populations_[currentPopulationIndexes_[CellDirection::north]]
			+ populations_[currentPopulationIndexes_[CellDirection::northWest]]
			+ populations_[currentPopulationIndexes_[CellDirection::west]]
			+ populations_[currentPopulationIndexes_[CellDirection::southWest]]
			+ populations_[currentPopulationIndexes_[CellDirection::south]]
			+ populations_[currentPopulationIndexes_[CellDirection::southEast]]
			+ populations_[currentPopulationIndexes_[CellDirection::rest]];
#elif EXPLICIT
		density_ =
			populations_[getArrayIndex(runIndex_, CellDirection::east)]
			+ populations_[getArrayIndex(runIndex_, CellDirection::northEast)]
			+ populations_[getArrayIndex(runIndex_, CellDirection::north)]
			+ populations_[getArrayIndex(runIndex_, CellDirection::northWest)]
			+ populations_[getArrayIndex(runIndex_, CellDirection::west)]
			+ populations_[getArrayIndex(runIndex_, CellDirection::southWest)]
			+ populations_[getArrayIndex(runIndex_, CellDirection::south)]
			+ populations_[getArrayIndex(runIndex_, CellDirection::southEast)]
			+ populations_[getArrayIndex(runIndex_, CellDirection::rest)];
#else
		const uint_t currentPopulationIndexOffset_ = nPopulations * runIndex_;
		density_ =
			populations_[currentPopulationIndexOffset_ + CellDirection::east]
			+ populations_[currentPopulationIndexOffset_ + CellDirection::northEast]
			+ populations_[currentPopulationIndexOffset_ + CellDirection::north]
			+ populations_[currentPopulationIndexOffset_ + CellDirection::northWest]
			+ populations_[currentPopulationIndexOffset_ + CellDirection::west]
			+ populations_[currentPopulationIndexOffset_ + CellDirection::southWest]
			+ populations_[currentPopulationIndexOffset_ + CellDirection::south]
			+ populations_[currentPopulationIndexOffset_ + CellDirection::southEast]
			+ populations_[currentPopulationIndexOffset_ + CellDirection::rest];
#endif
	}

	inline void computeVelocity() {
		
		//const uint_t nextPopulationIndexOffset_ = nPopulations * !runIndex_;

		/*	field_t densityInverse = 1 / density_;
			velocity_.at(SpatialDirection::x) =
				((populations_.at(currentPopulationIndexOffset_ + CellDirection::east)
					+ populations_.at(currentPopulationIndexOffset_ + CellDirection::northEast)
					+ populations_.at(currentPopulationIndexOffset_ + CellDirection::southEast)
					- populations_.at(currentPopulationIndexOffset_ + CellDirection::west)
					- populations_.at(currentPopulationIndexOffset_ + CellDirection::northWest)
					- populations_.at(currentPopulationIndexOffset_ + CellDirection::southWest)) * densityInverse) + ((1 * bodyForce.at(SpatialDirection::x) * dt) * 0.5 * densityInverse);

			velocity_.at(SpatialDirection::y) =
				((populations_.at(currentPopulationIndexOffset_ + CellDirection::north)
					+ populations_.at(currentPopulationIndexOffset_ + CellDirection::northEast)
					+ populations_.at(currentPopulationIndexOffset_ + CellDirection::northWest)
					- populations_.at(currentPopulationIndexOffset_ + CellDirection::south)
					- populations_.at(currentPopulationIndexOffset_ + CellDirection::southEast)
					- populations_.at(currentPopulationIndexOffset_ + CellDirection::southWest)) * densityInverse) + ((1 * bodyForce.at(SpatialDirection::y) * dt) * 0.5 * densityInverse);*/

		/*velocity_.at(SpatialDirection::x) =
			((populations_.at(currentPopulationIndexOffset_ + CellDirection::east)
				+ populations_.at(currentPopulationIndexOffset_ + CellDirection::northEast)
				+ populations_.at(currentPopulationIndexOffset_ + CellDirection::southEast)
				- populations_.at(currentPopulationIndexOffset_ + CellDirection::west)
				- populations_.at(currentPopulationIndexOffset_ + CellDirection::northWest)
				- populations_.at(currentPopulationIndexOffset_ + CellDirection::southWest)) / density_) + ((1 * bodyForce.at(SpatialDirection::x) * dt) / (2 * density_));

		velocity_.at(SpatialDirection::y) =
			((populations_.at(currentPopulationIndexOffset_ + CellDirection::north)
				+ populations_.at(currentPopulationIndexOffset_ + CellDirection::northEast)
				+ populations_.at(currentPopulationIndexOffset_ + CellDirection::northWest)
				- populations_.at(currentPopulationIndexOffset_ + CellDirection::south)
				- populations_.at(currentPopulationIndexOffset_ + CellDirection::southEast)
				- populations_.at(currentPopulationIndexOffset_ + CellDirection::southWest)) / density_) + ((1 * bodyForce.at(SpatialDirection::y) * dt) / (2 * density_));*/

		//TODO: 0.5*
		/*velocity_[SpatialDirection::x] =
			((populations_[currentPopulationIndexOffset_ + CellDirection::east]
				+ populations_[currentPopulationIndexOffset_ + CellDirection::northEast]
				+ populations_[currentPopulationIndexOffset_ + CellDirection::southEast]
				- populations_[currentPopulationIndexOffset_ + CellDirection::west]
				- populations_[currentPopulationIndexOffset_ + CellDirection::northWest]
				- populations_[currentPopulationIndexOffset_ + CellDirection::southWest]) / density_) + ((tau * bodyForce[SpatialDirection::x] * dt) / (1 * density_));

		velocity_[SpatialDirection::y] =
			((populations_[currentPopulationIndexOffset_ + CellDirection::north]
				+ populations_[currentPopulationIndexOffset_ + CellDirection::northEast]
				+ populations_[currentPopulationIndexOffset_ + CellDirection::northWest]
				- populations_[currentPopulationIndexOffset_ + CellDirection::south]
				- populations_[currentPopulationIndexOffset_ + CellDirection::southEast]
				- populations_[currentPopulationIndexOffset_ + CellDirection::southWest]) / density_) + ((tau * bodyForce[SpatialDirection::y] * dt) / (1 * density_));*/
	
		const field_t densityInverse = 1./density_;
		
		/*for (uint_t index = 0; index < nPopulations; index++) {
			currentPopulationIndexes_[index] = currentPopulationIndexOffset_ + index;
		}*/

#if PRECOMPUTED
		velocity_[SpatialDirection::x] =
			densityInverse * ((populations_[currentPopulationIndexes_[CellDirection::east]]
				+ populations_[currentPopulationIndexes_[CellDirection::northEast]]
				+ populations_[currentPopulationIndexes_[CellDirection::southEast]]
				- populations_[currentPopulationIndexes_[CellDirection::west]]
				- populations_[currentPopulationIndexes_[CellDirection::northWest]]
				- populations_[currentPopulationIndexes_[CellDirection::southWest]]) + (tau * bodyForce[SpatialDirection::x] * dt));

		velocity_[SpatialDirection::y] =
			densityInverse * ((populations_[currentPopulationIndexes_[CellDirection::north]]
				+ populations_[currentPopulationIndexes_[CellDirection::northEast]]
				+ populations_[currentPopulationIndexes_[CellDirection::northWest]]
				- populations_[currentPopulationIndexes_[CellDirection::south]]
				- populations_[currentPopulationIndexes_[CellDirection::southEast]]
				- populations_[currentPopulationIndexes_[CellDirection::southWest]]) + (tau * bodyForce[SpatialDirection::y] * dt));
#elif EXPLICIT
		velocity_[SpatialDirection::x] =
			densityInverse * ((populations_[getArrayIndex(runIndex_, CellDirection::east)]
				+ populations_[getArrayIndex(runIndex_, CellDirection::northEast)]
				+ populations_[getArrayIndex(runIndex_, CellDirection::southEast)]
				- populations_[getArrayIndex(runIndex_, CellDirection::west)]
				- populations_[getArrayIndex(runIndex_, CellDirection::northWest)]
				- populations_[getArrayIndex(runIndex_, CellDirection::southWest)]) + (tau * bodyForce[SpatialDirection::x] * dt));

		velocity_[SpatialDirection::y] =
			densityInverse * ((populations_[getArrayIndex(runIndex_, CellDirection::north)]
				+ populations_[getArrayIndex(runIndex_, CellDirection::northEast)]
				+ populations_[getArrayIndex(runIndex_, CellDirection::northWest)]
				- populations_[getArrayIndex(runIndex_, CellDirection::south)]
				- populations_[getArrayIndex(runIndex_, CellDirection::southEast)]
				- populations_[getArrayIndex(runIndex_, CellDirection::southWest)]) + (tau * bodyForce[SpatialDirection::y] * dt));
#else
		const uint_t currentPopulationIndexOffset_ = nPopulations * runIndex_;
		velocity_[SpatialDirection::x] =
			densityInverse * ((populations_[currentPopulationIndexOffset_ + CellDirection::east]
				+ populations_[currentPopulationIndexOffset_ + CellDirection::northEast]
				+ populations_[currentPopulationIndexOffset_ + CellDirection::southEast]
				- populations_[currentPopulationIndexOffset_ + CellDirection::west]
				- populations_[currentPopulationIndexOffset_ + CellDirection::northWest]
				- populations_[currentPopulationIndexOffset_ + CellDirection::southWest]) + (tau * bodyForce[SpatialDirection::x] * dt));

		velocity_[SpatialDirection::y] =
			densityInverse * ((populations_[currentPopulationIndexOffset_ + CellDirection::north]
				+ populations_[currentPopulationIndexOffset_ + CellDirection::northEast]
				+ populations_[currentPopulationIndexOffset_ + CellDirection::northWest]
				- populations_[currentPopulationIndexOffset_ + CellDirection::south]
				- populations_[currentPopulationIndexOffset_ + CellDirection::southEast]
				- populations_[currentPopulationIndexOffset_ + CellDirection::southWest]) + (tau * bodyForce[SpatialDirection::y] * dt));	
#endif
		/*uint_t northEast = currentPopulationIndexOffset_ + CellDirection::northEast;
		uint_t northWest = currentPopulationIndexOffset_ + CellDirection::northWest;
		uint_t southWest = currentPopulationIndexOffset_ + CellDirection::southWest;
		uint_t southEast = currentPopulationIndexOffset_ + CellDirection::southEast;

		velocity_[SpatialDirection::x] =
			((populations_[currentPopulationIndexes_[CellDirection::east]]
				+ populations_[northEast]
				+ populations_[southEast]
				- populations_[currentPopulationIndexOffset_ + CellDirection::west]
				- populations_[northWest]
				- populations_[southWest]) + (tau * bodyForce[SpatialDirection::x] * dt)) / (1 * density_);

		velocity_[SpatialDirection::y] =
			((populations_[currentPopulationIndexOffset_ + CellDirection::north]
				+ populations_[northEast]
				+ populations_[northWest]
				- populations_[currentPopulationIndexOffset_ + CellDirection::south]
				- populations_[southEast]
				- populations_[southWest]) + (tau * bodyForce[SpatialDirection::y] * dt)) / (1 * density_);*/

		}

	void computePopulationsEq() {
		/*field_t ux = velocity_.at(SpatialDirection::x);
		field_t uy = velocity_.at(SpatialDirection::y);*/

		field_t ux = velocity_[SpatialDirection::x];
		field_t uy = velocity_[SpatialDirection::y];

		field_t uxSqr = ux * ux;
		field_t uySqr = uy * uy;
		field_t uxuy = ux * uy;
		field_t uSqr = uxSqr + uySqr;

		//// Weight fctors
		//const field_t weightFactorR = 2. / 9;	// Rest
		//const field_t weightFactorHV = 1. / 18;	// Horizontal/Vertical
		//const field_t weightFactorD = 1. / 36;		// Diagonal

		//// Weights
		//field_t weightR = density_ * weightFactorR;	// Rest
		//field_t weightHV = density_ * weightFactorHV;	// Horizontal/Vertical
		//field_t weightD = density_ * weightFactorD;		// Diagonal

		// Weights
		field_t weightR = (2 * density_) / 9;	// Rest
		field_t weightHV = density_ / 18;	// Horizontal/Vertical
		field_t weightD = density_ / 36;		// Diagonal

		//// Calculate the rest equlibrium field component
		//populationsEq_.at(CellDirection::rest) = weightR * (2 - (3 * uSqr));
		//// Calculate horizontal and vertical equlibrium field components
		//populationsEq_.at(CellDirection::east) = weightHV * (2 + (6 * ux) + (6 * uxSqr) - (3 * uySqr));
		//populationsEq_.at(CellDirection::north) = weightHV * (2 + (6 * uy) + (6 * uySqr) - (3 * uxSqr));
		//populationsEq_.at(CellDirection::west) = weightHV * (2 - (6 * ux) + (6 * uxSqr) - (3 * uySqr));
		//populationsEq_.at(CellDirection::south) = weightHV * (2 - (6 * uy) + (6 * uySqr) - (3 * uxSqr));
		//// Calculate diagonal equlibrium field components
		//populationsEq_.at(CellDirection::northEast) = weightD * (1 + (3 * (ux + uy)) + (9 * uxuy) + (3 * uSqr));
		//populationsEq_.at(CellDirection::northWest) = weightD * (1 - (3 * (ux - uy)) - (9 * uxuy) + (3 * uSqr));
		//populationsEq_.at(CellDirection::southWest) = weightD * (1 - (3 * (ux + uy)) + (9 * uxuy) + (3 * uSqr));
		//populationsEq_.at(CellDirection::southEast) = weightD * (1 + (3 * (ux - uy)) - (9 * uxuy) + (3 * uSqr));

		// Calculate the rest equlibrium field component
		populationsEq_[CellDirection::rest] = weightR * (2 - (3 * uSqr));
		//populationsEq_[CellDirection::rest] =

		// Calculate horizontal and vertical equlibrium field components
		populationsEq_[CellDirection::east] = weightHV * (2 + (6 * ux) + (6 * uxSqr) - (3 * uySqr));
		populationsEq_[CellDirection::north] = weightHV * (2 + (6 * uy) + (6 * uySqr) - (3 * uxSqr));
		populationsEq_[CellDirection::west] = weightHV * (2 - (6 * ux) + (6 * uxSqr) - (3 * uySqr));
		populationsEq_[CellDirection::south] = weightHV * (2 - (6 * uy) + (6 * uySqr) - (3 * uxSqr));

		// Calculate diagonal equlibrium field components
		populationsEq_[CellDirection::northEast] = weightD * (1 + (3 * (ux + uy)) + (9 * uxuy) + (3 * uSqr));
		populationsEq_[CellDirection::northWest] = weightD * (1 - (3 * (ux - uy)) - (9 * uxuy) + (3 * uSqr));
		populationsEq_[CellDirection::southWest] = weightD * (1 - (3 * (ux + uy)) + (9 * uxuy) + (3 * uSqr));
		populationsEq_[CellDirection::southEast] = weightD * (1 + (3 * (ux - uy)) - (9 * uxuy) + (3 * uSqr));

	}

	//void computeForcePopulations() {
	
	//	// TODO: Why do I have to divide by 2 in "tauFactor"???
	//	const field_t tauFactor = (1. - (dt / (2 * tau))) / 2;
	//	const field_t weightHV = tauFactor / 3;
	//	const field_t weightD = tauFactor / 12;
	//	const field_t weightR = 4 * tauFactor / 9;
	//	field_t Fx = bodyForce.at(SpatialDirection::x);
	//	field_t Fy = bodyForce.at(SpatialDirection::y);
	//	field_t ux = velocity_.at(SpatialDirection::x);
	//	field_t uy = velocity_.at(SpatialDirection::y);
	//	field_t FxUx = Fx * ux;
	//	field_t FxUy = Fx * uy;
	//	field_t FyUy = Fy * uy;
	//	field_t FyUx = Fy * ux;

	//	// Calculate horizontal and vertical force components
	//	forcePopulations_.at(CellDirection::east) = weightHV * (2 * Fx + 2 * FxUx - FyUy);
	//	forcePopulations_.at(CellDirection::north) = weightHV * (2 * Fy + 2 * FyUy - FxUx);
	//	forcePopulations_.at(CellDirection::west) = weightHV * (-2 * Fx + 2 * FxUx - FyUy);
	//	forcePopulations_.at(CellDirection::south) = weightHV * (-2 * Fy + 2 * FyUy - FxUx);

	//	// Calculate diagonal force components
	//	forcePopulations_.at(CellDirection::northEast) = weightD * (2 * Fx + 2 * Fy + 2 * FxUx + 3 * FyUx + 3 * FxUy + 2 * FyUy);
	//	forcePopulations_.at(CellDirection::northWest) = -weightD * (2 * Fx - 2 * Fy - 2 * FxUx + 3 * FyUx + 3 * FxUy - 2 * FyUy);
	//	forcePopulations_.at(CellDirection::southWest) = weightD * (-2 * Fx - 2 * Fy + 2 * FxUx + 3 * FyUx + 3 * FxUy + 2 * FyUy);
	//	forcePopulations_.at(CellDirection::southEast) = -weightD * (-2 * Fx + 2 * Fy - 2 * FxUx + 3 * FyUx + 3 * FxUy - 2 * FyUy);

	//	// Calculate rest force component
	//	forcePopulations_.at(CellDirection::rest) = -weightR * (3 * FxUx + 3 * FyUy);

	//	/*for (uint_t populationIndex = 0; populationIndex < nPopulations; populationIndex++) {
	//		std::cout << "\n forcePopulations_[" << populationIndex << "] = " << forcePopulations_.at(populationIndex);
	//	}
	//	std::cout << "\n\n";
	//	system("pause");		*/
	//}

	/*field_t dotProduct(const std::array<field_t, nDimensions> leftVector, const std::array<field_t, nDimensions> rightVector) const{
		field_t result = 0.;
		for (uint_t spatialDirection = 0; spatialDirection < nDimensions; spatialDirection++) {
			result += leftVector.at(spatialDirection) * rightVector.at(spatialDirection);
		}
		return result;*/
	
	//}

	//void computeForcePopulations() {
	
	//	const field_t csSqr = 1. / 3;
	//	///*const std::array<int, nPopulations * nDimensions> c = { 1,1,0,-1,-1,-1,0,1,0, 0,1,1,1,0,-1,-1,-1,0 };*/
	//	const std::array<field_t, nDimensions * nPopulations> c = { 1,0, 1,1, 0,1, -1,1, -1,0, -1,-1, 0,-1, 1,-1, 0,0 };
	//	const std::array<field_t, nPopulations> w = { (1. / 9), (1. / 36), (1. / 9), (1. / 36), (1. / 9), (1. / 36), (1. / 9), (1. / 36), (4. / 9) };
	//	
	//	field_t cDotF = 0;
	//	field_t cDotu = 0;
	//	field_t uDotF = 0;

	//	for (uint_t populationIndex = 0; populationIndex < nPopulations; populationIndex++) {
	//		std::array<field_t, nDimensions> cTemp = { c.at((2 * populationIndex) + 0) , c.at((2 * populationIndex) + 1) };
	//		cDotF = dotProduct(cTemp, bodyForce);
	//		cDotu = dotProduct(cTemp, velocity_);
	//		uDotF = dotProduct(velocity_, bodyForce);		

	//		//std::cout << "\ncDotF = " << cDotF;
	//		/*std::cout << "\ncDotu: ("
	//			<< c.at((2 * populationIndex) + 0) << " * " << velocity_.at(0) << ") + ("
	//			<< c.at((2 * populationIndex) + 1) << " * " << velocity_.at(1)
	//			<< ") = " << cDotu;*/
	//		//std::cout << "\tuDotF = " << uDotF;
	//		/*std::cout << "\nuDotF: ("
	//			<< velocity_.at(0) << " * " << bodyForce.at(0) << ") + ("
	//			<< velocity_.at(1) << " * " << bodyForce.at(1)
	//			<< ") = " << uDotF;*/
	//		//std::cout << "\n popultionIndex = " << populationIndex;

	//		forcePopulations_.at(populationIndex) = w.at(populationIndex) * (1 - (dt/(2 * tau))) * ((cDotF / csSqr) + (((cDotF * cDotu) - (csSqr * uDotF)) / (csSqr * csSqr)));
	//		std::cout << "\n forcePopulations_[" << populationIndex << "] = " << forcePopulations_.at(populationIndex);
	//	}
	//	//std::cout << "\n forcePopulations_[" << 8 << "] = " << forcePopulations_.at(8);
	//	std::cout << "\n\n";
	//	system("pause");		
	
	//}


	//*****************************************************************************************
	// Set functions
	/*void setPopulation(const uint_t populationIndex, const field_t population) {
		populations_.at(currentPopulationIndexOffset_ + populationIndex) = population;
	}*/
	void setPopulation(const uint_t populationIndex, const field_t population) {
		populations_[(nPopulations * runIndex_) + populationIndex] = population;
	}

	/*void setNextPopulation(const uint_t populationIndex, const field_t population) {
		populations_.at(nextPopulationIndexOffset_ + populationIndex) = population;
	}*/

	/*void setNextPopulation(const uint_t populationIndex, const field_t population) {
		populations_[nextPopulationIndexOffset_ + populationIndex] = population;
	}*/

	void setDensity(const field_t density) {
		density_ = density;
	}

	void setVelocity(const SpatialDirection direction, const field_t velocity) {
		velocity_.at(direction) = velocity;
	}


	//*****************************************************************************************
	// Get functions

	//Neighbours &getCellNeighbours() {
	//	return neighbours_;
	//}

	/*const std::string getPopulationsList(const bool runIndex) const {
		std::ostringstream populationsListStream;

		populationsListStream << "{" << std::setprecision(9) << populations_.at(getArrayIndex(runIndex, 0));
		for (uint_t i = 1; i < nPopulations; i++) {
			populationsListStream << ", " << populations_.at(getArrayIndex(runIndex, i));
		}
		populationsListStream << "}";
		return populationsListStream.str();
	}*/

	const std::string getPopulationsList(const bool runIndex) const {
		std::ostringstream populationsListStream;

		populationsListStream << "{" << std::setprecision(9) << populations_
			[(nPopulations * runIndex_) + 0];
		for (uint_t i = 1; i < nPopulations; i++) {
			populationsListStream << ", " << populations_[(nPopulations * runIndex_) + i];
		}
		populationsListStream << "}";
		return populationsListStream.str();
	}

	const std::string getNonEqPopulationsList(const bool runIndex) {
		std::ostringstream populationsListStream;
		std::array<field_t, nPopulations> nonEquilibriumPoulation;
		computePopulationsEq();
		for (uint_t populationIndex = 0; populationIndex < nPopulations; populationIndex++) {
			nonEquilibriumPoulation.at(populationIndex) = populations_[(nPopulations * runIndex_) + populationIndex] - populationsEq_.at(populationIndex);
			//nonEquilibriumPoulation.at(populationIndex) = populations_[getArrayIndex(runIndex, populationIndex)] - populationsEq_[populationIndex];
		}


		populationsListStream << "{" << std::setprecision(9) << nonEquilibriumPoulation.at(0);
		for (uint_t i = 1; i < nPopulations; i++) {
			populationsListStream << ", " << nonEquilibriumPoulation.at(i);
		}
		populationsListStream << "}";
		return populationsListStream.str();
	}


	/*const field_t getPopulation(const uint_t populationIndex) const {
		return populations_.at(currentPopulationIndexOffset_ + populationIndex);
	}*/

	const field_t getPopulation(const uint_t populationIndex) const {
		return populations_[(nPopulations * runIndex_) + populationIndex];
	}

	const field_t getDensity() const {
		return density_;
	}


	const field_t getVelocity(const SpatialDirection direction) const {
		return velocity_.at(direction);
	}

	const std::string getVelocityList() const {
		std::ostringstream velocityListStream;
		velocityListStream << "{" << std::setprecision(5) << velocity_.at(SpatialDirection::x);
		velocityListStream << ", " << velocity_.at(SpatialDirection::y) << std::setprecision(3) << "}";
		return velocityListStream.str();
	}

	char getCellTypeChar() const {
		return 'C';
	}


	CellType getCellType() const {
		return CellType::cell;
	}


	//const uint_t getNumberOfPopulations() const {
	//	return nPopulations;
	//}

	//const uint_t getNumberOfFieldDuplicates() const {
	//	return nFieldDuplicates;
	//}


	//*****************************************************************************************
	// Print functions

	/*void printPopulations(const bool runIndex) {
		std::cout << "populations at runIndex = " << runIndex << ": \n\n";
		for (uint_t i = 0; i < nPopulations; i++) {
			std::cout << populations_.at(getArrayIndex(runIndex, i)) << "\t";
		}
		std::cout << std::endl;
	}*/

	void printPopulationsEq(const bool runIndex) const {
		std::cout << "populationsEq at runIndex = " << runIndex << ": \n\n";
		for (uint_t i = 0; i < nPopulations; i++) {
			std::cout << populationsEq_.at(i) << "\t";
			//std::cout << populationsEq_[i] << "\t";
		}
		std::cout << std::endl;
	}
};

//**************************!!!!************************
// To be put in Cell.cpp

bool Cell::runIndex_{ false };
#if PRECOMPUTED
	uint_t Cell::currentPopulationIndexOffset_{ 0 };
	uint_t Cell::nextPopulationIndexOffset_{ 0 };
	std::array<uint_t, nPopulations> Cell::currentPopulationIndexes_{ 0,1,2,3,4,5,6,7,8 };
	std::array<uint_t, nPopulations> Cell::nextPopulationIndexes_{ 9,10,11,12,13,14,15,16,17 };
#elif EXPLICIT
#else
	uint_t Cell::currentPopulationIndexOffset_{ 0 };
	uint_t Cell::nextPopulationIndexOffset_{ 0 };
#endif


