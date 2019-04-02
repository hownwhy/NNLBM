#pragma once
#include "Globals.hpp"
#include "Cell/Cell.hpp"
#include "./Cell/BulkCell.hpp"
#include "./Cell/SolidCell.hpp"

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
	static const uint_t N_X_SOLID_CELLS = 0;
	static const uint_t N_Y_SOLID_CELLS = 1;
#elif TEST_TYPE & STREAM_TEST_BOX + CAVITY_TEST
	static const int N_X_SOLID_CELLS = 1;
	static const int N_Y_SOLID_CELLS = 1;
#endif
	static const uint_t xDimFluid = N_GRID_X_DIM_FLUID;
	static const uint_t yDimFluid = N_GRID_Y_DIM_FLUID;
	static const int xDimTotal = xDimFluid + 2 * N_X_SOLID_CELLS;
	static const int yDimTotal = yDimFluid + 2 * N_Y_SOLID_CELLS;

	//Geometry< xDimTotal, yDimTotal,GeometryType::pipe> geometry_;
	//std::array<uint_t, xDimTotal * yDimTotal> geometry;
	std::vector<uint_t> geometry;
	//std::array<std::shared_ptr<Cell>, xDimTotal * yDimTotal> grid;
	std::vector < std::shared_ptr<Cell>> grid;
	
public:
	Grid() {
		geometry.resize(xDimTotal * yDimTotal);
		grid.resize(xDimTotal * yDimTotal);
	}
	~Grid() = default;



public:

	//******************************************************************************************************************
	//******************************************************************************************************************
	// Grid possition 
	//******************************************************************************************************************

	// Transforms 2D possition to 1D. A 2D grid can then be represented by a 1D array or vector.
	inline uint_t gridPosition(const uint_t x_, const uint_t y_) const {
		return (y_ * xDimTotal) + x_;
	}

	// Returns the grid possition for a neighbour in a given direction
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
		uint_t cellType = 0x00;
		for (uint_t y = 0; y < yDimTotal; y++) {
			for (uint_t x = 0; x < xDimTotal; x++) {
				if (y < N_Y_SOLID_CELLS || y > yDimTotal - (1 + N_Y_SOLID_CELLS) || x < N_X_SOLID_CELLS || x > xDimTotal - (1 + N_X_SOLID_CELLS)) {
					cellType = CellType::solid;
				}				
				else if (y == (0 + N_Y_SOLID_CELLS) || y == yDimTotal - (1 + N_Y_SOLID_CELLS)) {

					cellType = CellType::wall;

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
					cellType = CellType::wall;
					//cellType = CellType::bulk;
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
		
		for (uint_t y = 0; y < yDimTotal; y++) {
			for (uint_t x = 0; x < xDimTotal; x++) {
				switch (geometry.at(gridPosition(x, y))) {
				case CellType::solid:
					grid.at(gridPosition(x, y)) = std::make_shared<SolidCell>();
					break;

				case CellType::bulk:
					grid.at(gridPosition(x, y)) = std::make_shared<BulkCell>();
					break;

				case CellType::wall:
					//grid.at(gridPosition(x, y)) = std::make_shared<WallCell>();
					grid.at(gridPosition(x, y)) = std::make_shared<BulkCell>();
					break;

				default:
					std::cout << "\nIn makeGrid(): No cell type match";
				}					
			}
		}
	}

	// Assign neighbour cells for each cell object
	// TODO: Generalize: Giving directions explicitly would probably be quite impractical for 3D implementations. 
	void linkNeighbours() const {
		// TODO: The margins were addeed as a quick fix to avoid dealing with wall cell neighbours.
		// Find a better solution to this problem.
		const uint_t xMargin = N_X_SOLID_CELLS;
		const uint_t yMargin = N_Y_SOLID_CELLS;
		for (uint_t y = yMargin; y < yDimTotal - yMargin; y++) {
			for (uint_t x = xMargin; x < xDimTotal - xMargin; x++) {
				grid.at(gridPosition(x, y))->setNeighbour(grid.at(gridNeigbourPossition(x, y, CellDirection::east)), CellDirection::east);
				grid.at(gridPosition(x, y))->setNeighbour(grid.at(gridNeigbourPossition(x, y, CellDirection::northEast)), CellDirection::northEast);
				grid.at(gridPosition(x, y))->setNeighbour(grid.at(gridNeigbourPossition(x, y, CellDirection::north)), CellDirection::north);
				grid.at(gridPosition(x, y))->setNeighbour(grid.at(gridNeigbourPossition(x, y, CellDirection::northWest)), CellDirection::northWest);
				grid.at(gridPosition(x, y))->setNeighbour(grid.at(gridNeigbourPossition(x, y, CellDirection::west)), CellDirection::west);
				grid.at(gridPosition(x, y))->setNeighbour(grid.at(gridNeigbourPossition(x, y, CellDirection::southWest)), CellDirection::southWest);
				grid.at(gridPosition(x, y))->setNeighbour(grid.at(gridNeigbourPossition(x, y, CellDirection::south)), CellDirection::south);
				grid.at(gridPosition(x, y))->setNeighbour(grid.at(gridNeigbourPossition(x, y, CellDirection::southEast)), CellDirection::southEast);

				/*grid.at(gridPosition(x, y))->getCellNeighbours().setNeighbour(CellDirection::east, grid.at(gridNeigbourPossition(x, y, CellDirection::east)));
				grid.at(gridPosition(x, y))->getCellNeighbours().setNeighbour(CellDirection::northEast, grid.at(gridNeigbourPossition(x, y, CellDirection::northEast)));
				grid.at(gridPosition(x, y))->getCellNeighbours().setNeighbour(CellDirection::north, grid.at(gridNeigbourPossition(x, y, CellDirection::north)));
				grid.at(gridPosition(x, y))->getCellNeighbours().setNeighbour(CellDirection::northWest, grid.at(gridNeigbourPossition(x, y, CellDirection::northWest)));
				grid.at(gridPosition(x, y))->getCellNeighbours().setNeighbour(CellDirection::west, grid.at(gridNeigbourPossition(x, y, CellDirection::west)));
				grid.at(gridPosition(x, y))->getCellNeighbours().setNeighbour(CellDirection::southWest, grid.at(gridNeigbourPossition(x, y, CellDirection::southWest)));
				grid.at(gridPosition(x, y))->getCellNeighbours().setNeighbour(CellDirection::south, grid.at(gridNeigbourPossition(x, y, CellDirection::south)));
				grid.at(gridPosition(x, y))->getCellNeighbours().setNeighbour(CellDirection::southEast, grid.at(gridNeigbourPossition(x, y, CellDirection::southEast)));*/
			}
		}
	}

	// Set rho and velocity for all cells, exxept for the "ghost" cells.
	//void gridInitialize(const bool runIndex) const {
	//	const uint_t xMargin = N_X_SOLID_CELLS;
	//	const uint_t yMargin = N_Y_SOLID_CELLS;
	//	const field_t rho = 1.0;
	//	const field_t xVelocity = 0.0;
	//	const field_t yVelocity = 0.0;
	//	//const 
	//	
	//	for (uint_t runIndex = 0; runIndex < 1; runIndex++) {
	//		for (uint_t y = yMargin; y < yDimTotal - yMargin; y++) {
	//			for (uint_t x = xMargin; x < xDimTotal - xMargin; x++) {
	//				grid.at(gridPosition(x, y))->initialize(runIndex, rho, xVelocity, yVelocity);
	//			}
	//		}
	//	}
	//}

	void gridInitialize(){
		const uint_t xMargin = 0;// N_X_SOLID_CELLS;
		const uint_t yMargin = 0;// N_Y_SOLID_CELLS;
		field_t rho = 1.0;
		field_t xVelocity = 0.0;
		field_t yVelocity = 0.0;
		field_t topPlateVelocity = F_LID_VELOCITY;
		
		for (uint_t y = yMargin; y < yDimTotal - yMargin; y++) {
			for (uint_t x = xMargin; x < xDimTotal - xMargin; x++) {
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
		
	//template<std::size_t SIZE>
	//static void collideAndPropagate2(std::array<std::shared_ptr<Cell>, SIZE> grid3, const uint_t beginIndex, const uint_t endIndex) {
	//	for (uint_t i = beginIndex; i < endIndex; ++i)
	//	{
	//		grid3.at(i)->collideAndPropagate();
	//		
	//	}
	//	//std::cout << std::endl << "Begin index = " << beginIndex;
	//}

	void collideAndPropagate(const bool runIndex) const {
		const uint_t xMargin = N_X_SOLID_CELLS;
		const uint_t yMargin = N_Y_SOLID_CELLS;
		Cell::setRunStep(runIndex);

		// C++17 parallel for loop
		std::for_each(			
			//std::execution::seq,
			//std::execution::par,
			std::execution::par_unseq,
			grid.begin(),
			grid.end(), 
			[](auto cell) {cell->collideAndPropagate(); }
		);

		//for (uint_t y = yMargin; y < yDimTotal - yMargin; y++) {
		//	//std::cout << std::endl;
		//	for (uint_t x = xMargin; x < xDimTotal - xMargin; x++) {
		//		//std::cout << "(" << x << "," << y << ")";
		//		grid.at(gridPosition(x, y))->collideAndPropagate();
		//	}
		//}

	}	

	//void moveBoundary(const bool runIndex) const {
	//	const uint_t xMargin = N_X_SOLID_CELLS;
	//	const uint_t yMargin = N_Y_SOLID_CELLS;
	//	for (uint_t y = yMargin; y < yDimTotal - yMargin; y++) {
	//		for (uint_t x = xMargin; x < xDimTotal - xMargin; x++) {
	//			if (y == 0 + yMargin) {
	//				grid.at(gridPosition(x, y))->addMovingBoundaryTerm(runIndex);
	//				/*std::cout << "addMoving....";
	//				system("pause");*/
	//			}				
	//		}
	//	}
	//}

	//void computeAllRho(const bool runIndex) const {
	//	const uint_t xMargin = N_X_SOLID_CELLS;
	//	const uint_t yMargin = N_Y_SOLID_CELLS;
	//	for (uint_t y = yMargin; y < yDimTotal - yMargin; y++) {
	//		for (uint_t x = xMargin; x < xDimTotal - xMargin; x++) {
	//			grid.at(gridPosition(x, y))->computeDensity(runIndex);
	//		}
	//	}
	//}


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
	void printNeighboursCellType(const uint_t x, const uint_t y) const {
		std::shared_ptr<Cell> tempCell;
		char E, NE, N, NW, W, SW, S, SE, R;


		if ((tempCell = grid.at(gridPosition(x, y))->getNeighbour(CellDirection::east)) != nullptr) {
			E = tempCell->getCellTypeChar();
		}
		else {
			E = 'N';
		}
		if ((tempCell = grid.at(gridPosition(x, y))->getNeighbour(CellDirection::northEast)) != nullptr) {
			NE = tempCell->getCellTypeChar();
		}
		else {
			NE = 'N';
		}
		if ((tempCell = grid.at(gridPosition(x, y))->getNeighbour(CellDirection::north)) != nullptr) {
			N = tempCell->getCellTypeChar();
		}
		else {
			N = 'N';
		}
		if ((tempCell = grid.at(gridPosition(x, y))->getNeighbour(CellDirection::northWest)) != nullptr) {
			NW = tempCell->getCellTypeChar();
		}
		else {
			NW = 'N';
		}
		if ((tempCell = grid.at(gridPosition(x, y))->getNeighbour(CellDirection::west)) != nullptr) {
			W = tempCell->getCellTypeChar();
		}
		else {
			W = 'N';
		}
		if ((tempCell = grid.at(gridPosition(x, y))->getNeighbour(CellDirection::southWest)) != nullptr) {
			SW = tempCell->getCellTypeChar();
		}
		else {
			SW = 'N';
		}
		if ((tempCell = grid.at(gridPosition(x, y))->getNeighbour(CellDirection::south)) != nullptr) {
			S = tempCell->getCellTypeChar();
		}
		else {
			S = 'N';
		}
		if ((tempCell = grid.at(gridPosition(x, y))->getNeighbour(CellDirection::southEast)) != nullptr) {
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
		const uint_t xMargin = N_X_SOLID_CELLS;
		const uint_t yMargin = N_Y_SOLID_CELLS;
		for (uint_t y = yMargin; y < yDimTotal - yMargin; y++) {
			for (uint_t x = xMargin; x < xDimTotal - xMargin; x++) {
				printNeighboursCellType(x, y);
			}
		}
	}


	//******************************************************************************************************************
	// Cell Population - print routines

	void printCellPopulation(const bool runIndex, const uint_t x, const uint_t y) const {
		std::cout << std::endl;
		std::cout << grid.at(gridPosition(x, y))->getPopulation(CellDirection::northWest)
			<< " " << grid.at(gridPosition(x, y))->getPopulation(CellDirection::north)
			<< " " << grid.at(gridPosition(x, y))->getPopulation(CellDirection::northEast) << std::endl;

		std::cout << " \\|/ " << std::endl;

		std::cout << grid.at(gridPosition(x, y))->getPopulation(CellDirection::west)
			<< "-" << grid.at(gridPosition(x, y))->getPopulation(CellDirection::rest)
			<< "-" << grid.at(gridPosition(x, y))->getPopulation(CellDirection::east) << std::endl;

		std::cout << " /|\\ " << std::endl;

		std::cout << grid.at(gridPosition(x, y))->getPopulation(CellDirection::southWest)
			<< " " << grid.at(gridPosition(x, y))->getPopulation(CellDirection::south)
			<< " " << grid.at(gridPosition(x, y))->getPopulation(CellDirection::southEast) << std::endl;
	}

	void  printCellPopulation(const bool runIndex) const {
		const uint_t xMargin = N_X_SOLID_CELLS;
		const uint_t yMargin = N_Y_SOLID_CELLS;
		for (uint_t y = yMargin; y < yDimTotal - yMargin; y++) {
			for (uint_t x = xMargin; x < xDimTotal - xMargin; x++) {
				printCellPopulation(runIndex, x, y);
			}
		}
	}


	//******************************************************************************************************************
	// Cell rho - print routine
	void printCellRho(const bool runIndex) const {
		/*std::cout << std::endl;*/
		for (uint_t y = 0; y < yDimTotal; y++) {
			for (uint_t x = 0; x < xDimTotal; x++) {
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
		for (uint_t y = 0; y < yDimTotal; y++) {
			for (uint_t x = 0; x < xDimTotal; x++) {
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
		for (uint_t y = 0; y < yDimTotal; y++) {
			populationLists += "{";
			for (uint_t x = 0; x < xDimTotal; x++) {
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
		for (uint_t y = 0; y < yDimTotal; y++) {
			populationLists += "{";
			for (uint_t x = 0; x < xDimTotal; x++) {
				populationLists += grid.at(gridPosition(x, y))->getNonEqPopulationsList(runIndex) + ((x < xDimTotal - 1) ? ",\n" : "");			
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
		for (uint_t y = 0; y < yDimTotal; y++) {
			densityStringStream << "{";
			for (uint_t x = 0; x < xDimTotal; x++) {
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

	std::shared_ptr<Cell> getCell(const uint_t x, const uint_t y) const {
		return grid.at(gridPosition(x, y));
	}
};