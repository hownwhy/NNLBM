#pragma once
#include "Globals.hpp"
#include "Neighbours.hpp"
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <array>
#include <memory>


const enum CellType : int
{
	cell = 0,
	bulkCell = 1,
	solidCell = 2,
	periodicBulk = 3,
	periodicSolid = 4
};

const enum SpatialDirection : int
{
	x = 0,
	y = 1
};

const enum SimType : int {
	streaming,
	poiseuille
};

// TODO: Where should these static constants be? Here it will be accessible to all who include the hpp file.
// If declareedd as a member, it has to be initialized in a cpp file...
static const int nPopulations = 9;		// Number of poulations for species
static const int nFieldDuplicates = 2;	// Number of field "duplicats" (for temporary storage)(probably no more than 2)
static const int nDimensions = 2;
static const field_t dt = 1;




//// Streaming example
//static field_t tau = 1000*dt;
//static std::array<field_t, 2> force = { 0, 0 };

//// Collision example
//static field_t tau = 1*dt;
//static std::array<field_t, 2> force = { 0, 0 };

// poiseuille flow
static const field_t tau = F_TAU;
static const std::array<field_t, 2> force = { F_BODY_FORCE_X, F_BODY_FORCE_Y };



// TODO: See what can be done with templates instead of virtual functions
class Cell {

protected:
	Neighbours neighbours_;
	std::array<field_t, nFieldDuplicates * nPopulations> populations_;
	std::array<field_t, nPopulations> populationsEq_;
	std::array<field_t, nPopulations> forcePopulations_;
	std::array<field_t, nDimensions> velocity_;
	field_t density_;
	
	


public:

	virtual ~Cell() = default;
	Cell() = default;

	// Returns the 1D array index which depend on the 2D (runIndex, direction).
	static inline int getArrayIndex(bool runIndex, int direction) {
		return (runIndex * nPopulations) + direction;
	}
	//*****************************************************************************************
	// Cell initialization

	void initialize(const field_t density, const field_t xVelocity, const field_t yVelocity) {
		const bool runIndex = 0;
		density_ = density;
		velocity_[SpatialDirection::x] = xVelocity;
		velocity_[SpatialDirection::y] = yVelocity;
		//std::cout << "density = " << density_ << "\txVelocity = " << velocity_[SpatialDirection::x] << "\tyVelocity = " << velocity_[SpatialDirection::y] << std::endl;
		computePopulationsEq();
		//printPopulationsEq(0);
		//printPopulations(0);
		std::copy(populationsEq_.begin(), populationsEq_.end(), populations_.begin());
		//std::copy(populationsEq_.begin() + (runIndex * nPopulations), populationsEq_.end() - (!runIndex * nPopulations), populations_.begin() + (runIndex * nPopulations));
		//printPopulations(0);
		//system("pause");
	}	

	void initializeDensity(const field_t density) {
		const bool runIndex = 0;
		density_ = density;
		computePopulationsEq();
		//std::copy(populationsEq_.at(runIndex * nPopulations), populationsEq_.at((runIndex * nPopulations) + nPopulations), populations_.at(runIndex * nPopulations));
		std::copy(populationsEq_.begin(), populationsEq_.end(), populations_.begin());
		//std::copy(populationsEq_.begin() + (runIndex * nPopulations), populationsEq_.end() - (!runIndex * nPopulations), populations_.begin() + (runIndex * nPopulations));
	}

	void initializeVelocity(const SpatialDirection direction, const field_t velocity) {
		const bool runIndex = 0;
		velocity_[direction] = velocity;
		computePopulationsEq();		
		std::copy(populationsEq_.begin(), populationsEq_.end(), populations_.begin());
		//std::copy(populationsEq_.begin() + (runIndex * nPopulations), populationsEq_.end() - (!runIndex * nPopulations), populations_.begin() + (runIndex * nPopulations));
	}

	//*****************************************************************************************
	// Do functions

	// Propagate populations to neighbouring cells (notice the: !runIndex) used
	// so that it will not overwrite populations that will be used in the current run (loop).	
	// TODO: The local variables used for readability may (or may not) decrease performance. Find out. 
	void propageteTo(const bool runIndex) const {
		field_t currentPopulation;
		std::shared_ptr<Cell> targetCell;

		for (int direction = 0; direction < nDirections; direction++) {
			currentPopulation = populations_[getArrayIndex(runIndex, direction)];
			targetCell = neighbours_.getNeighbour(direction);
			targetCell->setReceived(!runIndex, direction, currentPopulation);
		}
	}

	virtual void setReceived(const bool runIndex, const int populationIndex, const field_t fieldValue) = 0;

	virtual void collide(const bool runIndex) = 0;

	virtual void collideAndPropagate(const bool runIndex) = 0;

	void computeForcePopulations() {
		static const field_t weightHV = 1. / 3;
		static const field_t weightD = 1. / 12;
		static const field_t Fx = force[SpatialDirection::x];
		static const field_t Fy = force[SpatialDirection::y];
		static field_t ux = velocity_[SpatialDirection::x];
		static field_t uy = velocity_[SpatialDirection::y];
		static field_t FxUx = Fx * ux;
		static field_t FxUy = Fx * uy;
		static field_t FyUy = Fy * uy;
		static field_t FyUx = Fy * ux;


		// Calculate horizontal and vertical source field components
		forcePopulations_[CellDirection::east] = weightHV * (2 * Fx + 3 * FxUx);
		forcePopulations_[CellDirection::north] = weightHV * (2 * Fy + 3 * FyUy);
		forcePopulations_[CellDirection::west] = weightHV * (-2 * Fx + 3 * FxUx);
		forcePopulations_[CellDirection::south] = weightHV * (-2 * Fy + 3 * FyUy);

		// Calculate diagonal source field components
		forcePopulations_[CellDirection::northEast] = weightD * (2 * Fx + 2 * Fy + 3 * FxUx + 3 * FyUx + 3 * FxUy + 3 * FyUy);
		forcePopulations_[CellDirection::northWest] = -weightD * (2 * Fx - 2 * Fy - 3 * FxUx + 3 * FyUx + 3 * FxUy - 3 * FyUy);
		forcePopulations_[CellDirection::southWest] = weightD * (-2 * Fx - 2 * Fy + 3 * FxUx + 3 * FyUx + 3 * FxUy + 3 * FyUy);
		forcePopulations_[CellDirection::southEast] = -weightD * (-2 * Fx + 2 * Fy - 3 * FxUx + 3 * FyUx + 3 * FxUy - 3 * FyUy);

	}

	void computeDensity(const bool runIndex) {
		density_ =
			populations_[getArrayIndex(runIndex, CellDirection::east)] 
			+ populations_[getArrayIndex(runIndex, CellDirection::northEast)] 
			+ populations_[getArrayIndex(runIndex, CellDirection::north)] 
			+ populations_[getArrayIndex(runIndex, CellDirection::northWest)] 
			+ populations_[getArrayIndex(runIndex, CellDirection::west)] 
			+ populations_[getArrayIndex(runIndex, CellDirection::southWest)] 
			+ populations_[getArrayIndex(runIndex, CellDirection::south)] 
			+ populations_[getArrayIndex(runIndex, CellDirection::southEast)] 
			+ populations_[getArrayIndex(runIndex, CellDirection::rest)];
	}

	void computeVelocity(const bool runIndex, const std::array<field_t, 2> bodyForce) {
		velocity_[SpatialDirection::x] =
			((populations_[getArrayIndex(runIndex, CellDirection::east)]
				+ populations_[getArrayIndex(runIndex, CellDirection::northEast)]
				+ populations_[getArrayIndex(runIndex, CellDirection::southEast)]
				- populations_[getArrayIndex(runIndex, CellDirection::west)]
				- populations_[getArrayIndex(runIndex, CellDirection::northWest)]
				- populations_[getArrayIndex(runIndex, CellDirection::southWest)]) / density_) + bodyForce[SpatialDirection::x] / (2 * density_);

		velocity_[SpatialDirection::y] =
			((populations_[getArrayIndex(runIndex, CellDirection::north)]
				+ populations_[getArrayIndex(runIndex, CellDirection::northEast)]
				+ populations_[getArrayIndex(runIndex, CellDirection::northWest)]
				- populations_[getArrayIndex(runIndex, CellDirection::south)]
				- populations_[getArrayIndex(runIndex, CellDirection::southEast)]
				- populations_[getArrayIndex(runIndex, CellDirection::southWest)]) / density_) + bodyForce[SpatialDirection::y] / (2 * density_);
	}

	void computePopulationsEq() {
		field_t ux = velocity_[ SpatialDirection::x];
		field_t uy = velocity_[SpatialDirection::y];
		field_t uxSqr = ux * ux;
		field_t uySqr = uy * uy;
		field_t uxuy = ux * uy;
		field_t uSqr = uxSqr + uySqr;

		// Weights
		field_t weightR = (2 * density_) / 9;	// Rest
		field_t weightHV = density_ / 18;	// Horizontal/Vertical
		field_t weightD = density_ / 36;		// Diagonal

		// Calculate the rest equlibrium field component
		populationsEq_[CellDirection::rest] = weightR * (2 - (3 * uSqr));

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


	//*****************************************************************************************
	// Set functions
	void setPopulation(const bool runIndex, const int populationIndex, const field_t fieldValue) {
		int arrayIndex = getArrayIndex(runIndex, populationIndex);
		populations_[arrayIndex] = fieldValue;
	}

	void setDensity(const field_t density) {
		density_ = density;
	}

	void setVelocity(const SpatialDirection direction, const field_t velocity) {
		velocity_[direction] = velocity;
	}


	//*****************************************************************************************
	// Get functions

	Neighbours &getCellNeighbours() {
		return neighbours_;
	}

	const std::string getPopulationsList(const bool runIndex) const {
		std::ostringstream populationsListStream;
		populationsListStream << "{" << std::setprecision(9) << populations_[getArrayIndex(runIndex, 0)];
		for (int i = 1; i < nPopulations; i++) {
			populationsListStream << ", " << populations_[getArrayIndex(runIndex, i)];
		}
		populationsListStream << "}";
		return populationsListStream.str();
	}


	const field_t getPolulation(const bool runIndex, const int populationIndex) const {
		int arrayIndex = getArrayIndex(runIndex, populationIndex);
		//std::cout << "getPopulation ARRAY_INDEX : " << arrayIndex << std::endl;;
		assert(("getPopulations: arrayIndex is negative", arrayIndex >= 0));
		assert(("getPopulations: arrayIndex to high", arrayIndex < nFieldDuplicates * nPopulations));
		return populations_[arrayIndex];
	}

	const field_t getDensity() const {
		return density_;
	}


	const field_t getVelocity(const SpatialDirection direction) const {
		return velocity_[direction];
	}

	const std::string getVelocityList() {
		std::ostringstream velocityListStream;
		velocityListStream << "{" << std::setprecision(3) << velocity_[SpatialDirection::x];
		velocityListStream << ", " << velocity_[SpatialDirection::y] << std::setprecision(3) << "}";
		return velocityListStream.str();
	}

	virtual char getCellTypeChar() const {
		return 'C';
	}


	virtual CellType getCellType() const {
		return CellType::cell;
	}


	const int getNumberOfPopulations() {
		return nPopulations;
	}

	const int getNumberOfFieldDuplicates() {
		return nFieldDuplicates;
	}


	//*****************************************************************************************
	// Print functions
	
	void printPopulations(const bool runIndex) {
		std::cout << "populations at runIndex = " << runIndex << ": \n\n";
		for (int i = 0; i < nPopulations; i++) {
			std::cout << populations_[getArrayIndex(runIndex, i)] << "\t";
		}
		std::cout << std::endl;
	}

	void printPopulationsEq(const bool runIndex) {
		std::cout << "populationsEq at runIndex = " << runIndex << ": \n\n";
		for(int i = 0; i < nPopulations; i++) {
			std::cout << populationsEq_[i] << "\t";
		}
		std::cout << std::endl;
	}
};