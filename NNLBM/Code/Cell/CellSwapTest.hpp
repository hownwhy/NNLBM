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


static const uint_t nStreamDirections = 8;	// Number of directions populations can propagate
static const uint_t nCellDirections = 9;		// Number of poulations for species
//static const uint_t nFieldDuplicates = 2;	// Number of field "duplicats" (for temporary storage)(probably no more than 2)
static const uint_t nDimensions = 2;
static const field_t dt = 1.;


static const field_t tau = F_TAU;
static const field_t tauInverse = 1. / tau;
static const field_t dtOverTau = dt / tau;
static const field_t oneMinusDtOverTau = 1. - dtOverTau;
static const std::array<field_t, 2> bodyForce = { F_BODY_FORCE_X, F_BODY_FORCE_Y };
static const std::array<uint_t, nCellDirections> reverseDirectionIndex = { 4,5,6,7,0,1,2,3 };


// TODO: See what can be done with templates instead of virtual functions
class Cell {

private:
	std::array<Cell*, nStreamDirections> neighbours_;
	//std::array<field_t, nFieldDuplicates * nCellDirections> populations_;
	//field_t populations_[2 * nCellDirections]; // An array holding two sets of populations
	field_t p1[nCellDirections];
	field_t p2[nCellDirections];
	field_t *populationsNow_ = p1;
	field_t *populationsNext_ = p2;
	//std::array<field_t, nCellDirections> forcePopulations_;
	std::array<field_t, nDimensions> velocity_;
	field_t density_;
	CellType cellType_;
	//field_t tau;

	/*static bool runIndex_;
	static std::array<uint_t, nCellDirections> currentPopulationIndexes_;
	static std::array<uint_t, nCellDirections> nextPopulationIndexes_;*/

	// Add velocity from moving walls
	field_t solidToBulkVelocityTransfer(const field_t &population, const uint_t cellDirection) const {
		/*field_t solidToBulkVelocityTransfer(const field_t &population, const field_t density_,
			const field_t ux, const field_t uy, const uint_t cellDirection) const {*/

		const field_t csSqrInvers = 3.;
		// ! Rest population is EXCLUDED in the following definitions of c and w. !
		const std::array<int, nStreamDirections * nDimensions> c = { 1,1,0,-1,-1,-1,0,1, 0,1,1,1,0,-1,-1,-1 };
		//const std::array<field_t, nStreamDirections * nDimensions> c = { 1,1,0,-1,-1,-1,0,1, 0,1,1,1,0,-1,-1,-1 };
		const std::array<field_t, nStreamDirections> w = { (1. / 9), (1. / 36), (1. / 9), (1. / 36), (1. / 9), (1. / 36), (1. / 9), (1. / 36) };
		//const std::array<field_t, nStreamDirections> factors{ (6 / 9), (6 / 36), (6 / 9), (6 / 36), (6 / 9), (6 / 36), (6 / 9), (6 / 36) };
		return population - (2. * w.at(cellDirection) * density_ * csSqrInvers)
			//return population - (factors.at(cellDirection) * density_)
			/** ((c.at((SpatialDirection::x * nStreamDirections) + cellDirection) * ux)
				+ (c.at((SpatialDirection::y * nStreamDirections) + cellDirection) * uy));*/
			* ((c.at((SpatialDirection::x * nStreamDirections) + cellDirection) * velocity_[SpatialDirection::x])
				+ (c.at((SpatialDirection::y * nStreamDirections) + cellDirection) * velocity_[SpatialDirection::y]));
	}

public:
	Cell() {}
	Cell(const CellType cellType) : cellType_(cellType) {}

	~Cell() = default;



	//*****************************************************************************************
	// Cell initialization

	void initialize(const field_t density, const field_t ux, const field_t uy) {
		density_ = density;
		velocity_.at(SpatialDirection::x) = ux;
		velocity_.at(SpatialDirection::y) = uy;

		field_t populationsEq_[nCellDirections];
		computePopulationsEq(populationsEq_, density, ux, uy);
		for (uint_t cellDirection = 0; cellDirection < nCellDirections; cellDirection++) {
			populationsNow_[cellDirection] = populationsEq_[cellDirection];
			populationsNext_[cellDirection] = populationsEq_[cellDirection];
		}
	}


	//*****************************************************************************************
	// Do functions	

	static void setRunStep(const bool runIndex) {
		/*runIndex_ = runIndex;
		const uint_t currentPopulationIndexOffset_ = nCellDirections * runIndex;
		const uint_t nextPopulationIndexOffset_ = nCellDirections * !runIndex;

		currentPopulationIndexes_ = { 0 + currentPopulationIndexOffset_,	1 + currentPopulationIndexOffset_, 2 + currentPopulationIndexOffset_,
			3 + currentPopulationIndexOffset_, 4 + currentPopulationIndexOffset_, 5 + currentPopulationIndexOffset_,
			6 + currentPopulationIndexOffset_,	7 + currentPopulationIndexOffset_,	8 + currentPopulationIndexOffset_ };

		nextPopulationIndexes_ = { 0 + nextPopulationIndexOffset_,	1 + nextPopulationIndexOffset_, 2 + nextPopulationIndexOffset_,
			3 + nextPopulationIndexOffset_, 4 + nextPopulationIndexOffset_, 5 + nextPopulationIndexOffset_,
			6 + nextPopulationIndexOffset_,	7 + nextPopulationIndexOffset_,	8 + nextPopulationIndexOffset_ };*/
	}

	void swapPop() {
		std::swap(populationsNow_, populationsNext_);
	}

	void collideAndPropagate() {
		if (cellType_ == CellType::bulk) {
			collideAndPropagateBulk();
		}
		else {
			collideAndPropagateSolid();
		}
	}

	void collideAndPropagateBulk() {

		field_t density;
		computeDensity(density);
		field_t ux;
		field_t uy;
		computeVelocity(density, ux, uy);

		field_t populationsEq[nCellDirections];
		computePopulationsEq(populationsEq, density, ux, uy);

		// TODO: remains after half done optimization
		// Should probably not update these for every loop since
		// it causes unnecessary writes to RAM
		density_ = density;
		velocity_[SpatialDirection::x] = ux;
		velocity_[SpatialDirection::y] = uy;

		//computeForcePopulations();	

		for (uint_t cellDirection = 0; cellDirection < nStreamDirections; cellDirection++) {
			//populations_[currentPopulationIndexOffset_ + cellDirection] = (1 - (dt / tau)) * populations_[currentPopulationIndexOffset_ + cellDirection] + (dt / tau) * populationsEq_[cellDirection] // Collision			
			//	;// +forcePopulations_[cellDirection]; // Adding force term

			populationsNow_[cellDirection] = oneMinusDtOverTau * populationsNow_[cellDirection] + dtOverTau * populationsEq[cellDirection] // Collision			
				;// +forcePopulations_[cellDirection]; // Adding force term

			if (neighbours_[cellDirection]->getCellType() == CellType::bulk) {
				neighbours_[cellDirection]->setReceivedBulk(populationsNow_, cellDirection);
			}
			else {
				neighbours_[cellDirection]->setReceivedSolid(populationsNow_, populationsNext_, cellDirection);
			}

			//neighbours_[cellDirection]->setReceived(populations_, cellDirection);
		}

		//// Updating rest population (NO PROPAGATION!)
		//populations_[nextPopulationIndexOffset_ + CellDirection::rest] = (1 - (dt / tau)) * populations_[currentPopulationIndexOffset_ + CellDirection::rest] + (dt / tau) * populationsEq_[CellDirection::rest] // Collision			
		//	;//+ forcePopulations_[CellDirection::rest]; // Adding force term

		// Updating rest population (NO PROPAGATION!)
		populationsNext_[CellDirection::rest] = oneMinusDtOverTau * populationsNow_[CellDirection::rest] + dtOverTau * populationsEq[CellDirection::rest] // Collision			
			;//+ forcePopulations_[CellDirection::rest]; // Adding force term	

		
	}

	void collideAndPropagateSolid() {
		std::cout << "\nSolidCell::collideAndProppagate()";
		system("pause");
	}

	void setReceivedBulk(const field_t *sourcePopulationsNow, const uint_t cellDirection) {
		populationsNext_[cellDirection] = sourcePopulationsNow[cellDirection];
	}

	/*void setReceivedBulk(Cell *sourceCell, const uint_t cellDirection) {
		populationsNext_[cellDirection] = sourceCell->getPopulation(cellDirection);
	}
*/
	void setReceivedSolid(const field_t *sourcePopulationsNow, field_t *sourcePopulationsNext, const uint_t cellDirection) {
		sourcePopulationsNext[reverseDirectionIndex[cellDirection]]
			= solidToBulkVelocityTransfer(sourcePopulationsNow[cellDirection], cellDirection);
	}

	/*void setReceivedSolid(Cell *sourceCell, const uint_t cellDirection) {
		sourceCell->setPopulationNext(solidToBulkVelocityTransfer(sourceCell->getPopulation(cellDirection), cellDirection), reverseDirectionIndex[cellDirection]);
	}*/


	void computeDensity(field_t& density_) const {
		density_ =
			populationsNow_[CellDirection::east]
			+ populationsNow_[CellDirection::northEast]
			+ populationsNow_[CellDirection::north]
			+ populationsNow_[CellDirection::northWest]
			+ populationsNow_[CellDirection::west]
			+ populationsNow_[CellDirection::southWest]
			+ populationsNow_[CellDirection::south]
			+ populationsNow_[CellDirection::southEast]
			+ populationsNow_[CellDirection::rest];
	}

	void computeVelocity(const field_t density_, field_t& ux, field_t& uy) const {
		// Use this code if source population term is used
		/*const field_t halfDensityInverse = 0.5 / density_;
		velocity_[SpatialDirection::x] =
			halfDensityInverse * (2 * (populationsNow_[CellDirection::east]]
				+ populationsNow_[CellDirection::northEast]]
				+ populationsNow_[CellDirection::southEast]]
				- populationsNow_[CellDirection::west]]
				- populationsNow_[CellDirection::northWest]]
				- populationsNow_[CellDirection::southWest]]) + (bodyForce[SpatialDirection::x] * dt));

		velocity_[SpatialDirection::y] =
			halfDensityInverse * (2 * (populationsNow_[CellDirection::north]]
				+ populationsNow_[CellDirection::northEast]]
				+ populationsNow_[CellDirection::northWest]]
				- populationsNow_[CellDirection::south]]
				- populationsNow_[CellDirection::southEast]]
				- populationsNow_[CellDirection::southWest]]) + (bodyForce[SpatialDirection::y] * dt));*/

				// Use this code if the source population term is NOT used
		const field_t densityInverse = 1. / density_;
		//velocity_[SpatialDirection::x] =
		ux = densityInverse * ((populationsNow_[CellDirection::east]
			+ populationsNow_[CellDirection::northEast]
			+ populationsNow_[CellDirection::southEast]
			- populationsNow_[CellDirection::west]
			- populationsNow_[CellDirection::northWest]
			- populationsNow_[CellDirection::southWest]) + (tau * bodyForce[SpatialDirection::x] * dt));

		//velocity_[SpatialDirection::y] =
		uy = densityInverse * ((populationsNow_[CellDirection::north]
			+ populationsNow_[CellDirection::northEast]
			+ populationsNow_[CellDirection::northWest]
			- populationsNow_[CellDirection::south]
			- populationsNow_[CellDirection::southEast]
			- populationsNow_[CellDirection::southWest]) + (tau * bodyForce[SpatialDirection::y] * dt));
	}

	// !!!!!
	// About 20% speed increase by passing the local array "populationsEq[]" as a function parameter 
	// compared to using a permanently stored member array "populationsEq_[]"!
	void computePopulationsEq(field_t populationsEq_[], const field_t density, const field_t ux, const field_t uy) const {

		/*const field_t ux = velocity_[SpatialDirection::x];
		const field_t uy = velocity_[SpatialDirection::y];*/

		const field_t uxSqr = ux * ux;
		const field_t uySqr = uy * uy;
		const field_t uxuy = ux * uy;
		const field_t uSqr = uxSqr + uySqr;

		// Weight factors
		static const field_t weightFactorR = 2. / 9;	// Rest
		static const field_t weightFactorHV = 1. / 18;	// Horizontal/Vertical
		static const field_t weightFactorD = 1. / 36;		// Diagonal

		// Weights
		const field_t weightR = density * weightFactorR;	// Rest
		const field_t weightHV = density * weightFactorHV;	// Horizontal/Vertical
		const field_t weightD = density * weightFactorD;		// Diagonal

		const field_t eastWest = (2 + (6 * uxSqr) - (3 * uySqr));
		const field_t northSouth = (2 + (6 * uySqr) - (3 * uxSqr));

		const field_t nEsW = (1 + (9 * uxuy) + (3 * uSqr));
		const field_t nWsE = (1 - (9 * uxuy) + (3 * uSqr));

		//// Calculate the rest equlibrium field component
		populationsEq_[CellDirection::rest] = weightR * (2 - (3 * uSqr));

		// Calculate horizontal and vertical equlibrium field components
		populationsEq_[CellDirection::east] = weightHV * (eastWest + (6 * ux));
		populationsEq_[CellDirection::north] = weightHV * (northSouth + (6 * uy));
		populationsEq_[CellDirection::west] = weightHV * (eastWest - (6 * ux));
		populationsEq_[CellDirection::south] = weightHV * (northSouth - (6 * uy));

		// Calculate diagonal equlibrium field components
		populationsEq_[CellDirection::northEast] = weightD * (nEsW + (3 * (ux + uy)));
		populationsEq_[CellDirection::northWest] = weightD * (nWsE - (3 * (ux - uy)));
		populationsEq_[CellDirection::southWest] = weightD * (nEsW - (3 * (ux + uy)));
		populationsEq_[CellDirection::southEast] = weightD * (nWsE + (3 * (ux - uy)));


		// More explicit code (about 5% slower)
		/*const field_t ux = velocity_[SpatialDirection::x];
		const field_t uy = velocity_[SpatialDirection::y];

		const field_t uxSqr = ux * ux;
		const field_t uySqr = uy * uy;
		const field_t uxuy = ux * uy;
		const field_t uSqr = uxSqr + uySqr;

		// Weights
		const field_t weightR = (2 * density_) / 9;	// Rest
		const field_t weightHV = density_ / 18;	// Horizontal/Vertical
		const field_t weightD = density_ / 36;		// Diagonal

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
		populationsEq_[CellDirection::southEast] = weightD * (1 + (3 * (ux - uy)) - (9 * uxuy) + (3 * uSqr));*/
	}





	//*****************************************************************************************
	// Set functions

	void setNeighbour(Cell& cell, const CellDirection cellDirection) {
		neighbours_.at(cellDirection) = &cell;
	}

	void setPopulationNow(const field_t population, const uint_t cellDirection) {
		populationsNow_[cellDirection] = population;
	}

	void setPopulationNext(const field_t population, const uint_t cellDirection) {
		populationsNext_[cellDirection] = population;
	}

	void setDensity(const field_t density) {
		density_ = density;
	}

	void setVelocity(const SpatialDirection spatialDirection, const field_t velocity) {
		velocity_.at(spatialDirection) = velocity;
	}


	//*****************************************************************************************
	// Get functions

	const std::string getPopulationsList(const bool runIndex) const {
		std::ostringstream populationsListStream;

		populationsListStream << "{" << std::setprecision(9) << populationsNow_[0];
		for (uint_t cellDirection = 1; cellDirection < nCellDirections; cellDirection++) {
			populationsListStream << ", " << populationsNow_[cellDirection];
		}
		populationsListStream << "}";
		return populationsListStream.str();
	}

	const std::string getNonEqPopulationsList(const bool runIndex) const {
		field_t density;
		computeDensity(density);
		field_t ux;
		field_t uy;
		computeVelocity(density, ux, uy);
		field_t populationsEq[nCellDirections];
		computePopulationsEq(populationsEq, density, ux, uy);

		std::ostringstream populationsListStream;
		std::array<field_t, nCellDirections> nonEquilibriumPopulations;
		for (uint_t cellDirection = 0; cellDirection < nCellDirections; cellDirection++) {
			nonEquilibriumPopulations.at(cellDirection) = populationsNow_[cellDirection] - populationsEq[cellDirection];
		}

		populationsListStream << "{" << std::setprecision(9) << nonEquilibriumPopulations.at(0);
		for (uint_t i = 1; i < nCellDirections; i++) {
			populationsListStream << ", " << nonEquilibriumPopulations.at(i);
		}
		populationsListStream << "}";
		return populationsListStream.str();
	}

	const field_t getPopulation(const uint_t cellDirection) const {
		return populationsNow_[cellDirection];
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

	void setCellType(const CellType cellType) {
		cellType_ = cellType;
	}

	CellType getCellType() const {
		return cellType_;
	}

	char getCellTypeChar() const {
		char cellChar;
		if (cellType_ == CellType::bulk) {
			cellChar = 'B';
		}
		else if (cellType_ == CellType::solid) {
			cellChar = 'S';
		}
		else {
			cellChar = 'N';
		}
		return cellChar;
	}


	//*****************************************************************************************
	// Print functions

	void printPopulations(const bool runIndex) {
		std::cout << "populations at runIndex = " << runIndex << ": \n\n";
		for (uint_t i = 0; i < nCellDirections; i++) {
			std::cout << populationsNow_[i] << "\t";
		}
		std::cout << std::endl;
	}

	//void printPopulationsEq(const bool runIndex) const {
	//	std::cout << "populationsEq at runIndex = " << runIndex << ": \n\n";
	//	for (uint_t i = 0; i < nCellDirections; i++) {
	//		std::cout << populationsEq_.at(i) << "\t";
	//		//std::cout << populationsEq_[i] << "\t";
	//	}
	//	std::cout << std::endl;
	//}
};

//**************************!!!!************************
// To be put in Cell.cpp

//bool Cell::runIndex_{ false };
//std::array<uint_t, nCellDirections> Cell::currentPopulationIndexes_{ 0,1,2,3,4,5,6,7,8 };
//std::array<uint_t, nCellDirections> Cell::nextPopulationIndexes_{ 9,10,11,12,13,14,15,16,17 };






