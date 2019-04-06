#pragma once
#include "Globals.hpp"
#include "Cell/Cell.hpp"
//#include "./Cell/BulkCell.hpp"
//#include "./Cell/SolidCell.hpp"

#include <iostream>
#include <string>
#include <sstream>
//#include <iomanip>
//#include <fstream>
#include <memory>
#include <vector>
#include <thread>
#include <execution>



//template<size_t X_DIM_FLUID, size_t Y_DIM_FLUID, GeometryType GEOMETRY_TYPE>
class Grid {

private:
#if TEST_TYPE & STREAM_TEST_PIPE + COUETTE_TEST + POISEUILLE_TEST
	static const int N_X_SOLID_CELLS = 0;
	static const int N_Y_SOLID_CELLS = 1;
#elif TEST_TYPE & STREAM_TEST_BOX + CAVITY_TEST
	static const int N_X_SOLID_CELLS = 1;
	static const int N_Y_SOLID_CELLS = 1;
#endif
	static const uint_t xDimFluid = N_GRID_X_DIM_FLUID;
	static const uint_t yDimFluid = N_GRID_Y_DIM_FLUID;
	static const int xDimTotal = xDimFluid + 2 * N_X_SOLID_CELLS;
	static const int yDimTotal = yDimFluid + 2 * N_Y_SOLID_CELLS;

	//std::array<uint_t, xDimTotal * yDimTotal> geometry;
	std::vector<uint_t> geometry;
	__declspec(align(64)) std::vector<Cell> grid_;
	//std::vector<Cell> grid_;
	__declspec(align(64)) std::vector<Cell*> fluidGridPtrs;

public:
	Grid() {
		geometry.resize(xDimTotal * yDimTotal);
		grid_.resize(xDimTotal * yDimTotal);
	}
	~Grid() = default;



public:

	//******************************************************************************************************************
	//******************************************************************************************************************
	// Grid possition 
	//******************************************************************************************************************

	// Transforms 2D possition to 1D. A 2D grid_ can then be represented by a 1D array or vector.
	inline uint_t gridPosition(const uint_t x_, const uint_t y_) const {
		return (y_ * xDimTotal) + x_;
	}

	// Returns the grid_ possition for a neighbour in a given direction
	uint_t gridNeigbourPossition(const int x, const int y, const uint_t direction) const {

		int dx = 0;
		int dy = 0;

#if TEST_TYPE & STREAM_TEST_PIPE + COUETTE_TEST + POISEUILLE_TEST
		bool periodic = (x == 0 + N_X_SOLID_CELLS || x == xDimTotal - 1 - N_X_SOLID_CELLS);
#elif TEST_TYPE & STREAM_TEST_BOX + CAVITY_TEST
		bool periodic = 0;
#endif

		// non-periodic cell
		if (!periodic) {

			switch (direction) {
			case CellDirection::east:		dx = +1;
				break;
			case CellDirection::northEast:	dx = +1;
				dy = -1;
				break;
			case CellDirection::north:		dy = -1;
				break;
			case CellDirection::northWest:	dx = -1;
				dy = -1;
				break;
			case CellDirection::west:		dx = -1;
				break;
			case CellDirection::southWest:	dx = -1;
				dy = +1;
				break;
			case CellDirection::south:		dy = +1;
				break;
			case CellDirection::southEast:	dx = +1;
				dy = +1;
				break;
			case CellDirection::rest:
			default:						dx = 0;
				dy = 0;
			}
		}
		// periodic cell
		else {

			switch (direction) {
			case CellDirection::east:		dx = (x == xDimTotal - 1 - N_X_SOLID_CELLS) ? -(xDimTotal - 1 - (2 * N_X_SOLID_CELLS)) : +1;
				break;
			case CellDirection::northEast:	dx = (x == xDimTotal - 1 - N_X_SOLID_CELLS) ? -(xDimTotal - 1 - (2 * N_X_SOLID_CELLS)) : +1;
				dy = -1;
				break;
			case CellDirection::north:		dy = -1;
				break;
			case CellDirection::northWest:	dx = (x == 0 + N_X_SOLID_CELLS) ? (xDimTotal - 1 - (2 * N_X_SOLID_CELLS)) : -1;
				dy = -1;
				break;
			case CellDirection::west:		dx = (x == 0 + N_X_SOLID_CELLS) ? (xDimTotal - 1 - (2 * N_X_SOLID_CELLS)) : -1;
				break;
			case CellDirection::southWest:	dx = (x == 0 + N_X_SOLID_CELLS) ? (xDimTotal - 1 - (2 * N_X_SOLID_CELLS)) : -1;
				dy = +1;
				break;
			case CellDirection::south:		dy = +1;
				break;
			case CellDirection::southEast:	dx = (x == xDimTotal - 1 - N_X_SOLID_CELLS) ? -(xDimTotal - 1 - (2 * N_X_SOLID_CELLS)) : +1;
				dy = +1;
				break;
			case CellDirection::rest:
			default:						dx = 0;
				dy = 0;
			}
		}

		return gridPosition(x + dx, y + dy);
	}




	//******************************************************************************************************************
	//******************************************************************************************************************
	// Prepare and initialize
	//******************************************************************************************************************

	void makePipeGeometry() {
		uint_t cellType = 0x00;
		for (uint_t y = 0; y < yDimTotal; y++) {
			for (uint_t x = 0; x < xDimTotal; x++) {
				if (y < N_Y_SOLID_CELLS || y > yDimTotal - (1 + N_Y_SOLID_CELLS) || x < N_X_SOLID_CELLS || x > xDimTotal - (1 + N_X_SOLID_CELLS)) {
					cellType = CellType::solid;
				}
				else if (y == (0 + N_Y_SOLID_CELLS) || y == yDimTotal - (1 + N_Y_SOLID_CELLS)) {

					cellType = CellType::bulk;

				}
				else {
					cellType = CellType::bulk;
				}
				geometry.at((y * xDimTotal) + x) = cellType;
			}
		}
	}

	void makeBoxGeometry() {
		uint_t cellType = 0x00;
		for (uint_t y = 0; y < yDimTotal; y++) {
			for (uint_t x = 0; x < xDimTotal; x++) {

				if (y < (0 + N_Y_SOLID_CELLS) || y > yDimTotal - (1 + N_Y_SOLID_CELLS) || (x < 0 + N_X_SOLID_CELLS) || x > xDimTotal - (1 + N_X_SOLID_CELLS)) {
					cellType = CellType::solid;
				}
				else if (y == (0 + N_Y_SOLID_CELLS) || y == yDimTotal - (1 + N_Y_SOLID_CELLS) || (x == 0 + N_X_SOLID_CELLS) || x == xDimTotal - (1 + N_X_SOLID_CELLS)) {
					//cellType = CellType::wall;
					cellType = CellType::bulk;
				}
				else {
					cellType = CellType::bulk;
				}
				geometry.at((y * xDimTotal) + x) = cellType;
			}
		}
	}

	// An array of pointers to Cell type objects is filled with BulkCells and SolidCells 
	// according to the specified geometry.
	// TODO: To my understanding, this will not ensure that Cell objects are stored in a consecutive manner in memory.
	// grid_[] is an array of addresses, not an array of object data, as far as I understand.
	// Therefore, I should look into the possibility to allocate space for the objects themselves.
	// Will it work to allocate space for base class objects, while storing inherited class objects?
	// TODO: Also include: Neighbour linking through constructor as the cell objects are created.
	void makeGrid() {

		for (uint_t y = 0; y < yDimTotal; y++) {
			for (uint_t x = 0; x < xDimTotal; x++) {
				switch (geometry.at(gridPosition(x, y))) {
				case CellType::solid:
					grid_.at(gridPosition(x, y)).setCellType(CellType::solid);
					break;

				case CellType::bulk:
					grid_.at(gridPosition(x, y)).setCellType(CellType::bulk);
					fluidGridPtrs.push_back(&(grid_.at(gridPosition(x, y))));
					break;

				default:
					std::cout << "\nIn makeGrid(): No cell type match";
				}
			}
		}
	}

	// Assign neighbour cells for each cell object
	// TODO: Generalize: Giving directions explicitly would probably be quite impractical for 3D implementations. 
	void linkNeighbours() {
		// TODO: The margins were addeed as a quick fix to avoid dealing with wall cell neighbours.
		// Find a better solution to this problem.
		const uint_t xMargin = N_X_SOLID_CELLS;
		const uint_t yMargin = N_Y_SOLID_CELLS;
		for (uint_t y = yMargin; y < yDimTotal - yMargin; y++) {
			for (uint_t x = xMargin; x < xDimTotal - xMargin; x++) {
				grid_.at(gridPosition(x, y)).setNeighbour(grid_.at(gridNeigbourPossition(x, y, CellDirection::east)), CellDirection::east);
				grid_.at(gridPosition(x, y)).setNeighbour(grid_.at(gridNeigbourPossition(x, y, CellDirection::northEast)), CellDirection::northEast);
				grid_.at(gridPosition(x, y)).setNeighbour(grid_.at(gridNeigbourPossition(x, y, CellDirection::north)), CellDirection::north);
				grid_.at(gridPosition(x, y)).setNeighbour(grid_.at(gridNeigbourPossition(x, y, CellDirection::northWest)), CellDirection::northWest);
				grid_.at(gridPosition(x, y)).setNeighbour(grid_.at(gridNeigbourPossition(x, y, CellDirection::west)), CellDirection::west);
				grid_.at(gridPosition(x, y)).setNeighbour(grid_.at(gridNeigbourPossition(x, y, CellDirection::southWest)), CellDirection::southWest);
				grid_.at(gridPosition(x, y)).setNeighbour(grid_.at(gridNeigbourPossition(x, y, CellDirection::south)), CellDirection::south);
				grid_.at(gridPosition(x, y)).setNeighbour(grid_.at(gridNeigbourPossition(x, y, CellDirection::southEast)), CellDirection::southEast);
			}
		}
	}


	void gridInitialize() {
		const uint_t xMargin = 0;// N_X_SOLID_CELLS;
		const uint_t yMargin = 0;// N_Y_SOLID_CELLS;
		field_t rho = 1.0;
		field_t xVelocity = 0.0;
		field_t yVelocity = 0.0;
		field_t topPlateVelocity = F_LID_VELOCITY;

		for (uint_t y = yMargin; y < yDimTotal - yMargin; y++) {
			for (uint_t x = xMargin; x < xDimTotal - xMargin; x++) {
				if (y == 0 + yMargin) {
					grid_.at(gridPosition(x, y)).initialize(rho, topPlateVelocity, yVelocity);
				}
				else {
					grid_.at(gridPosition(x, y)).initialize(rho, xVelocity, yVelocity);
				}
			}
		}
	}

	//******************************************************************************************************************
	//******************************************************************************************************************
	// Do routines
	//******************************************************************************************************************


	void collideAndPropagate(const bool runIndex) {
		/*const uint_t xMargin = N_X_SOLID_CELLS;
		const uint_t yMargin = N_Y_SOLID_CELLS;*/
		Cell::setRunStep(runIndex);

		//fluidGridPtrs.at(2)->collideAndPropagate();

		// C++17 parallel for loop
		std::for_each(
			std::execution::seq,
			//std::execution::par,
			//std::execution::par_unseq,
			fluidGridPtrs.begin(),
			fluidGridPtrs.end(),
			[](auto cellPtr) {cellPtr->collideAndPropagate(); }
		);

	}

	//******************************************************************************************************************
	//******************************************************************************************************************
	// Print routines
	//******************************************************************************************************************

	void printGeometry() const {
		std::cout << std::endl;
		for (uint_t y = 0; y < yDimTotal; y++) {
			for (uint_t x = 0; x < xDimTotal; x++) {
				std::cout << geometry.at((y * xDimTotal) + x);
			}
			std::cout << std::endl;
		}
	}

	//******************************************************************************************************************
	// Cell type - print routine
	void printCellType() const {
		std::cout << std::endl;
		for (uint_t y = 0; y < yDimTotal; y++) {
			for (uint_t x = 0; x < xDimTotal; x++) {
				std::cout << grid_.at((y * xDimTotal) + x).getCellTypeChar();
				std::cout << " ";
			}
			std::cout << std::endl;
		}
	}


	std::string appendGridPoplulationsList(const bool runIndex, std::string& populationLists) const {
		//std::string populationLists;
		populationLists += ((populationLists == "") ? "{" : ",\n\n{");
		for (uint_t y = 0; y < yDimTotal; y++) {
			populationLists += "{";
			for (uint_t x = 0; x < xDimTotal; x++) {
				populationLists += grid_.at(gridPosition(x, y)).getPopulationsList(runIndex) + ((x < xDimTotal - 1) ? ",\n" : "");
			}
			populationLists += ((y < yDimTotal - 1) ? "},\n\n" : "}");
		}
		populationLists += "}";

		return populationLists;
	}


	void appendGridVelocityList(const bool runIndex, std::string& velocityLists) const {
		const uint_t xMargin = N_X_SOLID_CELLS;
		const uint_t yMargin = N_Y_SOLID_CELLS;
		std::ostringstream velocityStringStream;
		velocityStringStream << std::setprecision(4);
		velocityStringStream << velocityLists;
		if (velocityStringStream.str() == "") {
			velocityStringStream << "{" << xDimTotal - (2 * N_X_SOLID_CELLS) << "," << yDimTotal - (2 * N_Y_SOLID_CELLS) << "," << N_RUN << "," << F_TAU << "," << RE << "," << F_LID_VELOCITY << "," << F_BODY_FORCE_X << "," << F_BODY_FORCE_Y
				<< std::setprecision(4) << "},\n\n" << "{";
		}
		else {
			velocityStringStream << ",\n\n{";
		}
		/*velocityStringStream << velocityLists << ((velocityLists == "") ? "{" + std::to_string(xDimTotal) + "," + std::to_string(yDimTotal) + "},\n\n" + "{" : ",\n\n{");*/
		for (uint_t y = 0 + yMargin; y < yDimTotal - yMargin; y++) {
			velocityStringStream << "{";
			for (uint_t x = 0 + xMargin; x < xDimTotal - xMargin; x++) {
				//if (grid_.at(gridPosition(x, y)) != nullptr) {
				velocityStringStream << grid_.at(gridPosition(x, y)).getVelocityList() << ((x < xDimTotal - xMargin - 1) ? ",\n" : "");
				//}
			}

			velocityStringStream << ((y < yDimTotal - yMargin - 1) ? "},\n\n" : "}");
		}
		velocityStringStream << "}";

		velocityLists = velocityStringStream.str();

		//return velocityLists;

	}

	void appendGridDensityList(const bool runIndex, std::string& densityLists) const {
		std::ostringstream densityStringStream;
		densityStringStream << std::setprecision(3);
		densityStringStream << densityLists;
		if (densityStringStream.str() == "") {
			densityStringStream << "{" << xDimTotal << "," << yDimTotal << "," << F_TAU << "," << std::setprecision(12) << std::fixed << F_BODY_FORCE_X << "," << F_BODY_FORCE_Y
				<< std::setprecision(3) << std::defaultfloat << "},\n\n" << "{";
		}
		else {
			densityStringStream << ",\n\n{";
		}
		/*densityStringStream << densityLists << ((densityLists == "") ? "{" + std::to_string(xDimTotal) + "," + std::to_string(yDimTotal) + "},\n\n" + "{" : ",\n\n{");*/
		for (uint_t y = 0; y < yDimTotal; y++) {
			densityStringStream << "{";
			for (uint_t x = 0; x < xDimTotal; x++) {
				//if (grid_.at(gridPosition(x, y)) != nullptr) {
				densityStringStream << std::setprecision(12) << std::fixed << grid_.at(gridPosition(x, y)).getDensity() << std::setprecision(3) << std::defaultfloat << ((x < xDimTotal - 1) ? ",\n" : "");
				//}
			}

			densityStringStream << ((y < yDimTotal - 1) ? "},\n\n" : "}");
		}
		densityStringStream << "}";

		densityLists = densityStringStream.str();
	}
};
