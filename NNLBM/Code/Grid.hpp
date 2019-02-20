#pragma once
#include "Globals.hpp"
//#include "Geometry.hpp"
#include "Cell/Cell.hpp"
#include "./Cell/BulkCell.hpp"
#include "./Cell/WallCell.hpp"
#include "./Cell/SolidCell.hpp"

#include <iostream>
#include <string>
#include <sstream>
//#include <iomanip>
//#include <fstream>
#include <memory>

//template<size_t X_DIM_FLUID, size_t Y_DIM_FLUID, GeometryType GEOMETRY_TYPE>
class Grid {
	
private:
	static const int xDimFluid = N_GRID_X_DIM_FLUID;
	static const int yDimFluid = N_GRID_Y_DIM_FLUID;
	static const int N_X_SOLID_CELLS = 0;
	static const int N_Y_SOLID_CELLS = 1;
	static const int xDimTotal = xDimFluid + 2 * N_X_SOLID_CELLS;
	static const int yDimTotal = yDimFluid + 2 * N_Y_SOLID_CELLS;

	//Geometry< xDimTotal, yDimTotal,GeometryType::pipe> geometry_;
	std::array<int, xDimTotal * yDimTotal> geometry;
	std::array<std::shared_ptr<Cell>, xDimTotal * yDimTotal> grid;

public:
	Grid() {

	}
	~Grid() = default;



public:

	//******************************************************************************************************************
	//******************************************************************************************************************
	// Grid possition 
	//******************************************************************************************************************

	// Transforms 2D possition to 1D. A 2D grid can then be represented by a 1D array or vector.
	inline int gridPosition(const int x_, const int y_) const {
		return (y_ * xDimTotal) + x_;
	}

	// Returns the grid possition for a neighbour in a given direction
	int gridNeigbourPossition(const int x, const int y, const int direction) const {

		int dx = 0;
		int dy = 0;
		bool periodic = (x == 0 + N_X_SOLID_CELLS || x == xDimTotal - 1 - N_X_SOLID_CELLS);
		
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
			case CellDirection::northEast:	dx = (x == xDimTotal - 1- N_X_SOLID_CELLS) ? -(xDimTotal - 1 - (2 * N_X_SOLID_CELLS)) : +1;
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
		int cellType = 0x00;
		for (int y = 0; y < yDimTotal; y++) {
			for (int x = 0; x < xDimTotal; x++) {
				if (y < N_Y_SOLID_CELLS || y > yDimTotal - (1 + N_Y_SOLID_CELLS) || x < N_X_SOLID_CELLS || x > xDimTotal - (1 + N_X_SOLID_CELLS)) {
					cellType = CellType::solid;
				}				
				else if (y == (0 + N_Y_SOLID_CELLS) || y == yDimTotal - (1 + N_Y_SOLID_CELLS)) {
#if BOOK_BOUNCE_BACK
					cellType = CellType::wall;
#else
					cellType = CellType::solid;
#endif
				}
				else {
					cellType = CellType::bulk;
				}
				geometry.at((y * xDimTotal) + x) = cellType;
			}
		}		
	}

	void makeBoxGeometry() {
		int cellType = 0x00;
		for (int y = 0; y < yDimTotal; y++) {
			for (int x = 0; x < xDimTotal; x++) {

				if (y < (0 + N_Y_SOLID_CELLS) || y > yDimTotal - (1 + N_Y_SOLID_CELLS) || (x < 0 + N_X_SOLID_CELLS) || x > xDimTotal - (1 + N_X_SOLID_CELLS)) {
					cellType = CellType::solid;
				}
				else if (y == (0 + N_Y_SOLID_CELLS) || y == yDimTotal - (1 + N_Y_SOLID_CELLS) || (x == 0 + N_X_SOLID_CELLS) || x == xDimTotal - (1 + N_X_SOLID_CELLS)) {
					cellType = CellType::wall;
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
	// grid[] is an array of addresses, not an array of object data, as far as I understand.
	// Therefore, I should look into the possibility to allocate space for the objects themselves.
	// Will it work to allocate space for base class objects, while storing inherited class objects?
	// TODO: Also include: Neighbour linking through constructor as the cell objects are created.
	void makeGrid() {

		for (int y = 0; y < yDimTotal; y++) {
			for (int x = 0; x < xDimTotal; x++) {
				switch (geometry.at(gridPosition(x, y))) {
				case CellType::solid:
					grid.at(gridPosition(x, y)) = std::make_shared<SolidCell>();
					break;

				case CellType::bulk:
					grid.at(gridPosition(x, y)) = std::make_shared<BulkCell>();
					break;

				case CellType::wall:
					grid.at(gridPosition(x, y)) = std::make_shared<WallCell>();
					break;

				default:
					std::cout << "\nIn makeGrid(): No cell type match";
				}				
				
				for (int cellDirection = 0; cellDirection < nPopulations; cellDirection++) {
					grid.at(gridPosition(x, y))->setPopulation(0, cellDirection, 0);
					grid.at(gridPosition(x, y))->setPopulation(1, cellDirection, 0);
				}
			}
		}
	}

	// Assign neighbour cells for each cell object
	// TODO: Generalize: Giving directions explicitly would probably be quite impractical for 3D implementations. 
	void linkNeighbours() const {
		// TODO: The margins were addeed as a quick fix to avoid dealing with wall cell neighbours.
		// Find a better solution to this problem.
		const int xMargin = N_X_SOLID_CELLS;
		const int yMargin = N_Y_SOLID_CELLS;
		for (int y = yMargin; y < yDimTotal - yMargin; y++) {
			for (int x = xMargin; x < xDimTotal - xMargin; x++) {
				grid.at(gridPosition(x, y))->getCellNeighbours().setNeighbour(CellDirection::east, grid.at(gridNeigbourPossition(x, y, CellDirection::east)));
				grid.at(gridPosition(x, y))->getCellNeighbours().setNeighbour(CellDirection::northEast, grid.at(gridNeigbourPossition(x, y, CellDirection::northEast)));
				grid.at(gridPosition(x, y))->getCellNeighbours().setNeighbour(CellDirection::north, grid.at(gridNeigbourPossition(x, y, CellDirection::north)));
				grid.at(gridPosition(x, y))->getCellNeighbours().setNeighbour(CellDirection::northWest, grid.at(gridNeigbourPossition(x, y, CellDirection::northWest)));
				grid.at(gridPosition(x, y))->getCellNeighbours().setNeighbour(CellDirection::west, grid.at(gridNeigbourPossition(x, y, CellDirection::west)));
				grid.at(gridPosition(x, y))->getCellNeighbours().setNeighbour(CellDirection::southWest, grid.at(gridNeigbourPossition(x, y, CellDirection::southWest)));
				grid.at(gridPosition(x, y))->getCellNeighbours().setNeighbour(CellDirection::south, grid.at(gridNeigbourPossition(x, y, CellDirection::south)));
				grid.at(gridPosition(x, y))->getCellNeighbours().setNeighbour(CellDirection::southEast, grid.at(gridNeigbourPossition(x, y, CellDirection::southEast)));
			}
		}
	}

	// Set rho and velocity for all cells, exxept for the "ghost" cells.
	//void gridInitialize(const bool runIndex) const {
	//	const int xMargin = N_X_SOLID_CELLS;
	//	const int yMargin = N_Y_SOLID_CELLS;
	//	const field_t rho = 1.0;
	//	const field_t xVelocity = 0.0;
	//	const field_t yVelocity = 0.0;
	//	//const 
	//	
	//	for (int runIndex = 0; runIndex < 1; runIndex++) {
	//		for (int y = yMargin; y < yDimTotal - yMargin; y++) {
	//			for (int x = xMargin; x < xDimTotal - xMargin; x++) {
	//				grid.at(gridPosition(x, y))->initialize(runIndex, rho, xVelocity, yVelocity);
	//			}
	//		}
	//	}
	//}

	void gridInitialize() const{
		const int xMargin = N_X_SOLID_CELLS;
		const int yMargin = N_Y_SOLID_CELLS;
		field_t rho = 1.0;
		field_t xVelocity = 0.0;
		field_t yVelocity = 0.0;
		field_t topPlateVelocity = F_TOP_PLATE_VELOCITY;
		//field_t topPlateVelocity = 0.0;

		for (int y = yMargin; y < yDimTotal - yMargin; y++) {
			for (int x = xMargin; x < xDimTotal - xMargin; x++) {
				if (y == 0 + yMargin) {
					grid.at(gridPosition(x, y))->initialize(rho, topPlateVelocity, yVelocity);
				}
				else {
					grid.at(gridPosition(x, y))->initialize(rho, xVelocity, yVelocity);
				}
			}
		}
	}

	//******************************************************************************************************************
	//******************************************************************************************************************
	// Do routines
	//******************************************************************************************************************
	void propagate(const bool runIndex) const {
#if BOOK_BOUNCE_BACK
		const int xMargin = N_X_SOLID_CELLS;
		const int yMargin = N_Y_SOLID_CELLS;
#else
		const int xMargin = N_X_SOLID_CELLS;
		const int yMargin = N_Y_SOLID_CELLS;
#endif
		for (int y = yMargin; y < yDimTotal - yMargin; y++) {
			for (int x = xMargin; x < xDimTotal - xMargin; x++) {
				grid.at(gridPosition(x, y))->propageteTo(runIndex);
			}
		}
	}
	   
	void collide(const bool runIndex) const {
#if BOOK_BOUNCE_BACK
		const int xMargin = N_X_SOLID_CELLS;
		const int yMargin = N_Y_SOLID_CELLS;
#else
		const int xMargin = N_X_SOLID_CELLS;
		const int yMargin = 2*N_Y_SOLID_CELLS;
#endif
		for (int y = yMargin; y < yDimTotal - yMargin; y++) {
			for (int x = xMargin; x < xDimTotal - xMargin; x++) {
				grid.at(gridPosition(x, y))->collide(runIndex);
			}
		}
	}

	void collideAndPropagate(const bool runIndex) const {
		const int xMargin = N_X_SOLID_CELLS;
		const int yMargin = N_Y_SOLID_CELLS;
		for (int y = yMargin; y < yDimTotal - yMargin; y++) {
			for (int x = xMargin; x < xDimTotal - xMargin; x++) {
				grid.at(gridPosition(x, y))->collideAndPropagate(runIndex);
			}
		}
	}

	//void moveBoundary(const bool runIndex) const {
	//	const int xMargin = N_X_SOLID_CELLS;
	//	const int yMargin = N_Y_SOLID_CELLS;
	//	for (int y = yMargin; y < yDimTotal - yMargin; y++) {
	//		for (int x = xMargin; x < xDimTotal - xMargin; x++) {
	//			if (y == 0 + yMargin) {
	//				grid.at(gridPosition(x, y))->addMovingBoundaryTerm(runIndex);
	//				/*std::cout << "addMoving....";
	//				system("pause");*/
	//			}				
	//		}
	//	}
	//}

	void computeAllRho(const bool runIndex) const {
		const int xMargin = N_X_SOLID_CELLS;
		const int yMargin = N_Y_SOLID_CELLS;
		for (int y = yMargin; y < yDimTotal - yMargin; y++) {
			for (int x = xMargin; x < xDimTotal - xMargin; x++) {
				grid.at(gridPosition(x, y))->computeDensity(runIndex);
			}
		}
	}


	//******************************************************************************************************************
	//******************************************************************************************************************
	// Print routines
	//******************************************************************************************************************

	void printGeometry() const {
		std::cout << std::endl;
		for (int y = 0; y < yDimTotal; y++) {
			for (int x = 0; x < xDimTotal; x++) {
				std::cout << geometry.at((y * xDimTotal) + x);
			}
			std::cout << std::endl;
		}
	}

	//******************************************************************************************************************
	// Cell type - print routine
	void printCellType() const {
		std::cout << std::endl;
		for (int y = 0; y < yDimTotal; y++) {
			for (int x = 0; x < xDimTotal; x++) {
				if (grid.at(gridPosition(x, y)) != nullptr) {
					std::cout << grid.at((y * xDimTotal) + x)->getCellTypeChar();
				}
				else {
					std::cout << "N";
				}
				std::cout << " ";
			}
			std::cout << std::endl;
		}
	}

	//******************************************************************************************************************
	// Neighbours - print routines
	void printNeighboursCellType(const int x, const int y) const {
		std::shared_ptr<Cell> tempCell;
		char E, NE, N, NW, W, SW, S, SE, R;


		if ((tempCell = grid.at(gridPosition(x, y))->getCellNeighbours().getNeighbour(CellDirection::east)) != nullptr) {
			E = tempCell->getCellTypeChar();
		}
		else {
			E = 'N';
		}
		if ((tempCell = grid.at(gridPosition(x, y))->getCellNeighbours().getNeighbour(CellDirection::northEast)) != nullptr) {
			NE = tempCell->getCellTypeChar();
		}
		else {
			NE = 'N';
		}
		if ((tempCell = grid.at(gridPosition(x, y))->getCellNeighbours().getNeighbour(CellDirection::north)) != nullptr) {
			N = tempCell->getCellTypeChar();
		}
		else {
			N = 'N';
		}
		if ((tempCell = grid.at(gridPosition(x, y))->getCellNeighbours().getNeighbour(CellDirection::northWest)) != nullptr) {
			NW = tempCell->getCellTypeChar();
		}
		else {
			NW = 'N';
		}
		if ((tempCell = grid.at(gridPosition(x, y))->getCellNeighbours().getNeighbour(CellDirection::west)) != nullptr) {
			W = tempCell->getCellTypeChar();
		}
		else {
			W = 'N';
		}
		if ((tempCell = grid.at(gridPosition(x, y))->getCellNeighbours().getNeighbour(CellDirection::southWest)) != nullptr) {
			SW = tempCell->getCellTypeChar();
		}
		else {
			SW = 'N';
		}
		if ((tempCell = grid.at(gridPosition(x, y))->getCellNeighbours().getNeighbour(CellDirection::south)) != nullptr) {
			S = tempCell->getCellTypeChar();
		}
		else {
			S = 'N';
		}
		if ((tempCell = grid.at(gridPosition(x, y))->getCellNeighbours().getNeighbour(CellDirection::southEast)) != nullptr) {
			SE = tempCell->getCellTypeChar();
		}
		else {
			SE = 'N';
		}

		// Rest shows cell at (x, y)
		if ((tempCell = grid.at(gridPosition(x, y))) != nullptr) {
			R = tempCell->getCellTypeChar();
		}
		else {
			R = 'N';
		}

		std::cout << std::endl;
		std::cout << NW
			<< " " << N
			<< " " << NE << std::endl;

		std::cout << " \\|/ " << std::endl;

		std::cout << W
			<< "-" << R
			<< "-" << E << std::endl;

		std::cout << " /|\\ " << std::endl;

		std::cout << SW
			<< " " << S
			<< " " << SE << std::endl;
	}

	void printNeighboursCellType() const {
		const int xMargin = N_X_SOLID_CELLS;
		const int yMargin = N_Y_SOLID_CELLS;
		for (int y = yMargin; y < yDimTotal - yMargin; y++) {
			for (int x = xMargin; x < xDimTotal - xMargin; x++) {
				printNeighboursCellType(x, y);
			}
		}
	}


	//******************************************************************************************************************
	// Cell Population - print routines

	void printCellPopulation(const bool runIndex, const int x, const int y) const {
		std::cout << std::endl;
		std::cout << grid.at(gridPosition(x, y))->getPolulation(runIndex, CellDirection::northWest)
			<< " " << grid.at(gridPosition(x, y))->getPolulation(runIndex, CellDirection::north)
			<< " " << grid.at(gridPosition(x, y))->getPolulation(runIndex, CellDirection::northEast) << std::endl;

		std::cout << " \\|/ " << std::endl;

		std::cout << grid.at(gridPosition(x, y))->getPolulation(runIndex, CellDirection::west)
			<< "-" << grid.at(gridPosition(x, y))->getPolulation(runIndex, CellDirection::rest)
			<< "-" << grid.at(gridPosition(x, y))->getPolulation(runIndex, CellDirection::east) << std::endl;

		std::cout << " /|\\ " << std::endl;

		std::cout << grid.at(gridPosition(x, y))->getPolulation(runIndex, CellDirection::southWest)
			<< " " << grid.at(gridPosition(x, y))->getPolulation(runIndex, CellDirection::south)
			<< " " << grid.at(gridPosition(x, y))->getPolulation(runIndex, CellDirection::southEast) << std::endl;
	}

	void  printCellPopulation(const bool runIndex) const {
		const int xMargin = N_X_SOLID_CELLS;
		const int yMargin = N_Y_SOLID_CELLS;
		for (int y = yMargin; y < yDimTotal - yMargin; y++) {
			for (int x = xMargin; x < xDimTotal - xMargin; x++) {
				printCellPopulation(runIndex, x, y);
			}
		}
	}


	//******************************************************************************************************************
	// Cell rho - print routine
	void printCellRho(const bool runIndex) const {
		/*std::cout << std::endl;*/
		for (int y = 0; y < yDimTotal; y++) {
			for (int x = 0; x < xDimTotal; x++) {
				if (grid.at(gridPosition(x, y)) != nullptr) {
					std::cout << grid.at(gridPosition(x, y))->getDensity();
					std::cout << " - ";
				}
			}
			std::cout << std::endl;
		}
	}


	//******************************************************************************************************************
	// Cell velocity - print routine
	void printCellVelocity(const bool runIndex) const {
		std::cout << std::endl;
		for (int y = 0; y < yDimTotal; y++) {
			for (int x = 0; x < xDimTotal; x++) {
				if (grid.at(gridPosition(x, y)) != nullptr) {
					std::cout << "(" << grid.at(gridPosition(x, y))->getVelocity(SpatialDirection::x)
						<< ","
						<< grid.at(gridPosition(x, y))->getVelocity(SpatialDirection::y) << ")";
					std::cout << " - ";
				}
			}
			std::cout << std::endl;
		}
	}


	std::string appendGridPoplulationsList(const bool runIndex, std::string& populationLists) const {
		//std::string populationLists;
		populationLists += ((populationLists == "") ? "{" : ",\n\n{");
		for (int y = 0; y < yDimTotal; y++) {
			populationLists += "{";
			for (int x = 0; x < xDimTotal; x++) {
				populationLists += grid.at(gridPosition(x, y))->getPopulationsList(runIndex) + ((x < xDimTotal - 1) ? ",\n" : "");				
			}
			populationLists += ((y < yDimTotal - 1) ? "},\n\n" : "}");
		}
		populationLists += "}";

		return populationLists;
	}

	std::string appendGridNonEqPoplulationsList(const bool runIndex, std::string& populationLists) const {
		//std::string populationLists;
		populationLists += ((populationLists == "") ? "{" : ",\n\n{");
		for (int y = 0; y < yDimTotal; y++) {
			populationLists += "{";
			for (int x = 0; x < xDimTotal; x++) {
				populationLists += grid.at(gridPosition(x, y))->getNonEqPopulationsList(runIndex) + ((x < xDimTotal - 1) ? ",\n" : "");			
			}
			populationLists += ((y < yDimTotal - 1) ? "},\n\n" : "}");
		}
		populationLists += "}";

		return populationLists;
	}


	void appendGridVelocityList(const bool runIndex, std::string& velocityLists) const {
		const int xMargin = N_X_SOLID_CELLS;
		const int yMargin = N_Y_SOLID_CELLS;
		std::ostringstream velocityStringStream;
		velocityStringStream << std::setprecision(3);
		velocityStringStream << velocityLists;
		if (velocityStringStream.str() == "") {
			velocityStringStream << "{" << xDimTotal - (2 * N_X_SOLID_CELLS) << "," << yDimTotal - (2 * N_Y_SOLID_CELLS) << "," << F_TAU << "," << std::setprecision(12) << std::fixed << F_BODY_FORCE_X << "," << F_BODY_FORCE_Y
				<< std::setprecision(3) << std::defaultfloat << "},\n\n" << "{";
		}
		else {
			velocityStringStream << ",\n\n{";
		}
		/*velocityStringStream << velocityLists << ((velocityLists == "") ? "{" + std::to_string(xDimTotal) + "," + std::to_string(yDimTotal) + "},\n\n" + "{" : ",\n\n{");*/
		for (int y = 0 + yMargin; y < yDimTotal - yMargin; y++) {
			velocityStringStream << "{";
			for (int x = 0 + xMargin; x < xDimTotal - xMargin; x++) {
				//if (grid.at(gridPosition(x, y)) != nullptr) {
				velocityStringStream << grid.at(gridPosition(x, y))->getVelocityList() << ((x < xDimTotal - xMargin - 1) ? ",\n" : "");
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
		for (int y = 0; y < yDimTotal; y++) {
			densityStringStream << "{";
			for (int x = 0; x < xDimTotal; x++) {
				//if (grid.at(gridPosition(x, y)) != nullptr) {
				densityStringStream << std::setprecision(12) << std::fixed << grid.at(gridPosition(x, y))->getDensity() << std::setprecision(3) << std::defaultfloat << ((x < xDimTotal - 1) ? ",\n" : "");
				//}
			}

			densityStringStream << ((y < yDimTotal - 1) ? "},\n\n" : "}");
		}
		densityStringStream << "}";

		densityLists = densityStringStream.str();

		//return densityLists;

	}

	std::shared_ptr<Cell> getCell(const int x, const int y) const {
		return grid.at(gridPosition(x, y));
	}
};