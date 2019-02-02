#pragma once
#include "../Globals.hpp"
#include "../Neighbours.hpp"
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
//#include <array>
#include <memory>



const enum CellType : int
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
const enum  CellDirection : int
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

const enum SpatialDirection : int
{
	x = 0,
	y = 1
};


// TODO: Where should these static constants be? Here it will be accessible to all who include the hpp file.
// If declareedd as a member, it has to be initialized in a cpp file...
static const int nPopulations = 9;		// Number of poulations for species
static const int nFieldDuplicates = 2;	// Number of field "duplicats" (for temporary storage)(probably no more than 2)
static const int nDimensions = 2;
static const field_t dt = 1;

// poiseuille flow
static const field_t tau = F_TAU;
static const std::array<field_t, 2> bodyForce = { F_BODY_FORCE_X, F_BODY_FORCE_Y };

//Couette flow
//TODO: Make this a property of wall cells with an array of east, west, north and south wall velocities.
static field_t wallVelocityX = F_TOP_PLATE_VELOCITY;
static field_t wallVelocityY = 0;

// TODO: See what can be done with templates instead of virtual functions
class Cell {

protected:
	Neighbours neighbours_;
	std::array<field_t, nFieldDuplicates * nPopulations> populations_;
	std::array<field_t, nPopulations> populationsEq_;	
	std::array<field_t, nPopulations> forcePopulations_;
	std::array<field_t, nDimensions> velocity_;
	field_t density_;
	std::array<field_t, nPopulations> movingPlateTerm;

public:

	~Cell() = default;
	Cell() = default;

	// Returns the 1D array index which depend on the 2D (runIndex, direction).
	inline int getArrayIndex(const bool runIndex, const int direction) const {
		return (runIndex * nPopulations) + direction;
	}

	// This function gives the opposite direction of what you put in.
	// This is also the reason why I chose to not follow the convension
	// having the rest direction be the 0 direction. 	
	static const int threeHalfDirection = (3 * nDirections) / 2;
	inline int reverseDirectionIndex(const int direction) const {
		return ((direction + threeHalfDirection) % nDirections);
	}

	//*****************************************************************************************
	// Cell initialization

	void initialize(const field_t density, const field_t xVelocity, const field_t yVelocity) {
		density_ = density;
		velocity_[SpatialDirection::x] = xVelocity;
		velocity_[SpatialDirection::y] = yVelocity;
		//std::cout << "density = " << density_ << "\txVelocity = " << velocity_[SpatialDirection::x] << "\tyVelocity = " << velocity_[SpatialDirection::y] << std::endl;
		computePopulationsEq();
		std::copy(populationsEq_.begin(), populationsEq_.end(), populations_.begin());
		std::copy(populationsEq_.begin(), populationsEq_.end(), populations_.begin() + nPopulations);
	}

	void initializeDensity(const field_t density) {
		density_ = density;
		computePopulationsEq();
		std::copy(populationsEq_.begin(), populationsEq_.end(), populations_.begin());
	}

	void initializeVelocity(const SpatialDirection direction, const field_t velocity) {
		velocity_[direction] = velocity; // -((bodyForce[direction] * dt) / (density_));
		computePopulationsEq();
		std::copy(populationsEq_.begin(), populationsEq_.end(), populations_.begin());
	}

	//*****************************************************************************************
	// Do functions	

	//virtual void collide(const bool runIndex) = 0;
	void collide(const bool runIndex) {
		//std::cout << "\n\nDensity before computeDensity = " << density_;
		computeDensity(runIndex);
		//std::cout << "\n\nDensity after computeDensity = " << density_;
		computeVelocity(runIndex);
		computePopulationsEq();
		computeForcePopulations();
		int populationIndex = 0;
		for (int cellDirection = 0; cellDirection < nDirections; cellDirection++) {
			populationIndex = getArrayIndex(runIndex, cellDirection);
			populations_[populationIndex]
				= populations_[populationIndex] - (dt * (populations_[populationIndex] - populationsEq_[cellDirection]) / tau);
			// Add force term
			populations_[populationIndex] += forcePopulations_[cellDirection];
		}
	}

	virtual void propageteTo(const bool runIndex) = 0;

	virtual void collideAndPropagate(const bool runIndex) = 0;

	void setReceived(const bool runIndex, const int populationIndex, const field_t fieldValue) {
		int arrayIndex = getArrayIndex(runIndex, populationIndex);
		populations_[arrayIndex] = fieldValue;
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

	void computeVelocity(const bool runIndex) {
		velocity_[SpatialDirection::x] =
			((populations_[getArrayIndex(runIndex, CellDirection::east)]
				+ populations_[getArrayIndex(runIndex, CellDirection::northEast)]
				+ populations_[getArrayIndex(runIndex, CellDirection::southEast)]
				- populations_[getArrayIndex(runIndex, CellDirection::west)]
				- populations_[getArrayIndex(runIndex, CellDirection::northWest)]
				- populations_[getArrayIndex(runIndex, CellDirection::southWest)]) / density_) + ((1 * bodyForce[SpatialDirection::x] * dt) / (2 * density_));

		velocity_[SpatialDirection::y] =
			((populations_[getArrayIndex(runIndex, CellDirection::north)]
				+ populations_[getArrayIndex(runIndex, CellDirection::northEast)]
				+ populations_[getArrayIndex(runIndex, CellDirection::northWest)]
				- populations_[getArrayIndex(runIndex, CellDirection::south)]
				- populations_[getArrayIndex(runIndex, CellDirection::southEast)]
				- populations_[getArrayIndex(runIndex, CellDirection::southWest)]) / density_) + ((1 * bodyForce[SpatialDirection::y] * dt) / (2 * density_));
	}
	
	void computePopulationsEq() {
		field_t ux = velocity_[SpatialDirection::x];
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

	void computeForcePopulations() {
		// TODO: Why do I have to divide by 2 in "tauFactor"???
		const field_t tauFactor = (1. - (dt / (2 * tau))) / 2;
		const field_t weightHV = tauFactor / 3;
		const field_t weightD = tauFactor / 12;
		const field_t weightR = 4 * tauFactor / 9;
		field_t Fx = bodyForce[SpatialDirection::x];
		field_t Fy = bodyForce[SpatialDirection::y];
		field_t ux = velocity_[SpatialDirection::x];
		field_t uy = velocity_[SpatialDirection::y];
		field_t FxUx = Fx * ux;
		field_t FxUy = Fx * uy;
		field_t FyUy = Fy * uy;
		field_t FyUx = Fy * ux;

		// Calculate horizontal and vertical force components
		forcePopulations_[CellDirection::east] = weightHV * (2 * Fx + 2 * FxUx - FyUy);
		forcePopulations_[CellDirection::north] = weightHV * (2 * Fy + 2 * FyUy - FxUx);
		forcePopulations_[CellDirection::west] = weightHV * (-2 * Fx + 2 * FxUx - FyUy);
		forcePopulations_[CellDirection::south] = weightHV * (-2 * Fy + 2 * FyUy - FxUx);

		// Calculate diagonal force components
		forcePopulations_[CellDirection::northEast] = weightD * (2 * Fx + 2 * Fy + 2 * FxUx + 3 * FyUx + 3 * FxUy + 2 * FyUy);
		forcePopulations_[CellDirection::northWest] = -weightD * (2 * Fx - 2 * Fy - 2 * FxUx + 3 * FyUx + 3 * FxUy - 2 * FyUy);
		forcePopulations_[CellDirection::southWest] = weightD * (-2 * Fx - 2 * Fy + 2 * FxUx + 3 * FyUx + 3 * FxUy + 2 * FyUy);
		forcePopulations_[CellDirection::southEast] = -weightD * (-2 * Fx + 2 * Fy - 2 * FxUx + 3 * FyUx + 3 * FxUy - 2 * FyUy);

		// Calculate rest force component
		forcePopulations_[CellDirection::rest] = -weightR * (3 * FxUx + 3 * FyUy);

		/*for (int populationIndex = 0; populationIndex < nPopulations; populationIndex++) {
			std::cout << "\n forcePopulations_[" << populationIndex << "] = " << forcePopulations_[populationIndex];
		}
		std::cout << "\n\n";
		system("pause");		*/
	}

	/*field_t dotProduct(const std::array<field_t, nDimensions> leftVector, const std::array<field_t, nDimensions> rightVector) const{
		field_t result = 0.;
		for (int spatialDirection = 0; spatialDirection < nDimensions; spatialDirection++) {
			result += leftVector[spatialDirection] * rightVector[spatialDirection];
		}
		return result;
	}*/

	//void computeForcePopulations() {
	//	const field_t csSqr = 1. / 3;
	//	///*const std::array<int, nPopulations * nDimensions> c = { 1,1,0,-1,-1,-1,0,1,0, 0,1,1,1,0,-1,-1,-1,0 };*/
	//	const std::array<field_t, nDimensions * nPopulations> c = { 1,0, 1,1, 0,1, -1,1, -1,0, -1,-1, 0,-1, 1,-1, 0,0 };
	//	const std::array<field_t, nPopulations> w = { (1. / 9), (1. / 36), (1. / 9), (1. / 36), (1. / 9), (1. / 36), (1. / 9), (1. / 36), (4. / 9) };
	//	
	//	field_t cDotF = 0;
	//	field_t cDotu = 0;
	//	field_t uDotF = 0;

	//	for (int populationIndex = 0; populationIndex < nPopulations; populationIndex++) {
	//		std::array<field_t, nDimensions> cTemp = { c[(2 * populationIndex) + 0] , c[(2 * populationIndex) + 1] };
	//		cDotF = dotProduct(cTemp, bodyForce);
	//		cDotu = dotProduct(cTemp, velocity_);
	//		uDotF = dotProduct(velocity_, bodyForce);		

	//		//std::cout << "\ncDotF = " << cDotF;
	//		/*std::cout << "\ncDotu: ("
	//			<< c[(2 * populationIndex) + 0] << " * " << velocity_[0] << ") + ("
	//			<< c[(2 * populationIndex) + 1] << " * " << velocity_[1]
	//			<< ") = " << cDotu;*/
	//		//std::cout << "\tuDotF = " << uDotF;
	//		/*std::cout << "\nuDotF: ("
	//			<< velocity_[0] << " * " << bodyForce[0] << ") + ("
	//			<< velocity_[1] << " * " << bodyForce[1]
	//			<< ") = " << uDotF;*/
	//		//std::cout << "\n popultionIndex = " << populationIndex;

	//		forcePopulations_[populationIndex] = w[populationIndex] * (1 - (dt/(2 * tau))) * ((cDotF / csSqr) + (((cDotF * cDotu) - (csSqr * uDotF)) / (csSqr * csSqr)));
	//		std::cout << "\n forcePopulations_[" << populationIndex << "] = " << forcePopulations_[populationIndex];
	//	}
	//	//std::cout << "\n forcePopulations_[" << 8 << "] = " << forcePopulations_[8];

	//	std::cout << "\n\n";
	//	system("pause");		
	//}


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

	const std::string getNonEqPopulationsList(const bool runIndex) {
		std::ostringstream populationsListStream;
		std::array<field_t, nPopulations> nonEquilibriumPoulation;
		computePopulationsEq();
		for (int populationIndex = 0; populationIndex < nPopulations; populationIndex++) {
			nonEquilibriumPoulation[populationIndex] = populations_[getArrayIndex(runIndex, populationIndex)] - populationsEq_[populationIndex];
		}


		populationsListStream << "{" << std::setprecision(9) << nonEquilibriumPoulation[0];
		for (int i = 1; i < nPopulations; i++) {
			populationsListStream << ", " << nonEquilibriumPoulation[i];
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
		velocityListStream << "{" << std::setprecision(5) << velocity_[SpatialDirection::x];
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
		for (int i = 0; i < nPopulations; i++) {
			std::cout << populationsEq_[i] << "\t";
		}
		std::cout << std::endl;
	}
};