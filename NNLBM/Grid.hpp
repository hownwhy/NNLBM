#pragma once
#include "Globals.hpp"
#include "Cell.hpp"
#include "BulkCell.hpp"
#include "SolidCell.hpp"
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <memory>



class Grid {

public:
	Grid() = default;
	~Grid() = default;

private:
	static const int xDim = N_GRID_X_DIM;
	static const int yDim = N_GRID_Y_DIM;

	std::array<int, xDim * yDim> geometry;
	std::array<std::shared_ptr<Cell>, xDim * yDim> grid;

	
	//******************************************************************************************************************
	//******************************************************************************************************************
	// Grid possition 
	//******************************************************************************************************************

	// Transforms 2D possition to 1D. A 2D grid can then be represented by a 1D array or vector.
	inline int gridPosition(const int x_, const int y_) const{
		return (y_ * xDim) + x_;
	}

	// Returns the grid possition for a neighbour in a given direction
	int gridNeigbourPossition(const int x, const int y, const int direction) const {

		int dx = 0;
		int dy = 0;
		int gridPos = 0;

		switch (geometry[(y * xDim) + x]) {
		case CellType::bulkCell:
		case CellType::solidCell:
			
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
			break;

		case CellType::periodicBulk:
		case CellType::periodicSolid:
			switch (direction) {
			case CellDirection::east:		(x == xDim - 2) ? dx = -(xDim - 3) : dx = +1;
				break;
			case CellDirection::northEast:	(x == xDim - 2) ? dx = -(xDim - 3) : dx = +1;
											dy = -1;
				break;
			case CellDirection::north:		dy = -1;
				break;
			case CellDirection::northWest:	(x == 1) ? dx = xDim - 3 : dx = -1;
				dy = -1;
				break;
			case CellDirection::west:		(x == 1) ? dx = xDim - 3 : dx = -1;
				break;
			case CellDirection::southWest:	(x == 1) ? dx = xDim - 3 : dx = -1;
				dy = +1;
				break;
			case CellDirection::south:		dy = +1;
				break;
			case CellDirection::southEast:	(x == xDim - 2) ? dx = -(xDim - 3) : dx = +1;
				dy = +1;
				break;
			case CellDirection::rest:
			default:						dx = 0;
				dy = 0;
			}
			
			break;
			
		}
			gridPos = gridPosition(x + dx, y + dy);
			/*std::cout << "pos (" << x << "," << y << ") : " << "(" << x + dx << ", " << y + dy << ")" << std::endl;
			system("pause");*/
		return gridPos;
		
	}


public:

	//******************************************************************************************************************
	//******************************************************************************************************************
	// Prepare and initialize
	//******************************************************************************************************************

	// Makes a description of a rectangular box filled with fluid (1 = solid, and 0 = fluid/bulk).
	// The reulting array "geometry" for a 6 x 5 box will look like this :
	// 111111
	// 100001
	// 100001
	// 100001
	// 111111	
	// -----------------------
	void makeGeometry() {
		for (int y = 0; y < yDim; y++) {
			for (int x = 0; x < xDim; x++) {
				if (y <= 1 || y >= yDim - 2 || x <= 0 || x >= xDim - 1) {
					if (y == 1 || y == yDim - 2 || x == 0 || x == xDim - 1) {
						geometry[(y * xDim) + x] = CellType::periodicSolid;
					}
					else {
						geometry[(y * xDim) + x] = CellType::solidCell;
					}
				}
				else if (x == 1 || x == xDim - 2) {
					geometry[(y * xDim) + x] = CellType::periodicBulk;
				}
				else {
					geometry[(y * xDim) + x] = CellType::bulkCell;
				}
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

		for (int y = 0; y < yDim; y++) {
			for (int x = 0; x < xDim; x++) {
				if (geometry[gridPosition(x, y)] == CellType::bulkCell) {
					grid[gridPosition(x, y)] = std::make_shared<BulkCell>();
					for (int i = 0; i < nFieldDuplicates; i++) {
						for (int cellDirection = 0; cellDirection < nPopulations; cellDirection++) {
							grid[gridPosition(x, y)]->setPopulation(i, cellDirection, 0);
						}
					}
				}
				else if (geometry[gridPosition(x, y)] == CellType::periodicBulk) {
					grid[gridPosition(x, y)] = std::make_shared<BulkCell>();
					for (int i = 0; i < nFieldDuplicates; i++) {
						for (int cellDirection = 0; cellDirection < nPopulations; cellDirection++) {
							grid[gridPosition(x, y)]->setPopulation(i, cellDirection, 0);
						}
					}
				}
				else if (geometry[gridPosition(x, y)] == CellType::periodicSolid) {
					grid[gridPosition(x, y)] = std::make_shared<SolidCell>();
					for (int i = 0; i < nFieldDuplicates; i++) {
						for (int cellDirection = 0; cellDirection < nPopulations; cellDirection++) {
							grid[gridPosition(x, y)]->setPopulation(i, cellDirection, 0);
						}
					}
				}
				else if (geometry[gridPosition(x, y)] == CellType::solidCell) {
					grid[gridPosition(x, y)] = std::make_shared<SolidCell>();
					for (int i = 0; i < nFieldDuplicates; i++) {
						for (int cellDirection = 0; cellDirection < nPopulations; cellDirection++) {
							grid[gridPosition(x, y)]->setPopulation(i, cellDirection, 0);
						}
					}
				}
			}
		}
	}

	// Assign neighbour cells for each cell object
	// TODO: Generalize: Giving directions explicitly would probably be quite impractical for 3D implementations. 
	void linkNeighbours() const {
		// TODO: The margins were addeed as a quick fix to avoid dealing with wall cell neighbours.
		// Find a better solution to this problem.
		const int xMargin = 1;
		const int yMargin = 1;
		for (int y = yMargin; y < yDim - yMargin; y++) {
			for (int x = xMargin; x < xDim- xMargin; x++) {				
				grid[gridPosition(x, y)]->getCellNeighbours().setNeighbour(CellDirection::east,			grid[gridNeigbourPossition(x, y, CellDirection::east)]);
				grid[gridPosition(x, y)]->getCellNeighbours().setNeighbour(CellDirection::northEast,	grid[gridNeigbourPossition(x, y, CellDirection::northEast)]);
				grid[gridPosition(x, y)]->getCellNeighbours().setNeighbour(CellDirection::north,		grid[gridNeigbourPossition(x, y, CellDirection::north)]);
				grid[gridPosition(x, y)]->getCellNeighbours().setNeighbour(CellDirection::northWest,	grid[gridNeigbourPossition(x, y, CellDirection::northWest)]);
				grid[gridPosition(x, y)]->getCellNeighbours().setNeighbour(CellDirection::west,			grid[gridNeigbourPossition(x, y, CellDirection::west)]);
				grid[gridPosition(x, y)]->getCellNeighbours().setNeighbour(CellDirection::southWest,	grid[gridNeigbourPossition(x, y, CellDirection::southWest)]);
				grid[gridPosition(x, y)]->getCellNeighbours().setNeighbour(CellDirection::south,		grid[gridNeigbourPossition(x, y, CellDirection::south)]);
				grid[gridPosition(x, y)]->getCellNeighbours().setNeighbour(CellDirection::southEast,	grid[gridNeigbourPossition(x, y, CellDirection::southEast)]);				
			}
		}
	}	   	

	// Set rho and velocity for all cells, exxept for the "ghost" cells.
	//void gridInitialize(const bool runIndex) const {
	//	const int xMargin = 1;
	//	const int yMargin = 1;
	//	const field_t rho = 1.0;
	//	const field_t xVelocity = 0.0;
	//	const field_t yVelocity = 0.0;
	//	//const 
	//	
	//	for (int runIndex = 0; runIndex < 1; runIndex++) {
	//		for (int y = yMargin; y < yDim - yMargin; y++) {
	//			for (int x = xMargin; x < xDim - xMargin; x++) {
	//				grid[gridPosition(x, y)]->initialize(runIndex, rho, xVelocity, yVelocity);
	//			}
	//		}
	//	}
	//}

	void gridInitialize(const bool runIndex) const {
		const int xMargin = 1;
		const int yMargin = 1;
		field_t rho = 1.0;
		field_t xVelocity = 0.0;
		field_t yVelocity = 0.0;	
		
		for (int y = yMargin; y < yDim - yMargin; y++) {
			for (int x = xMargin; x < xDim - xMargin; x++) {
				grid[gridPosition(x, y)]->initialize(rho, xVelocity, yVelocity);
			}
		}		
	}
	
	//******************************************************************************************************************
	//******************************************************************************************************************
	// Do routines
	//******************************************************************************************************************
	void propagate(const bool runIndex) const{
		const int xMargin = 1;
		const int yMargin = 1;
		for (int y = yMargin; y < yDim - yMargin; y++) {
			for (int x = xMargin; x < xDim - xMargin; x++) {
				grid[gridPosition(x, y)]->propageteTo(runIndex);
			}
		}
	}

	void collide(const bool runIndex) const{
		const int xMargin = 1;
		const int yMargin = 1;
		for (int y = yMargin; y < yDim - yMargin; y++) {
			for (int x = xMargin; x < xDim - xMargin; x++) {
				grid[gridPosition(x, y)]->collide(runIndex);
			}
		}
	}

	void collideAndPropagate(const bool runIndex) const {
		const int xMargin = 1;
		const int yMargin = 1;
		for (int y = yMargin; y < yDim - yMargin; y++) {
			for (int x = xMargin; x < xDim - xMargin; x++) {
				grid[gridPosition(x, y)]->collideAndPropagate(runIndex);
			}
		}
	}

	void computeAllRho(const bool runIndex) const{
		const int xMargin = 1;
		const int yMargin = 1;
		for (int y = yMargin; y < yDim - yMargin; y++) {
			for (int x = xMargin; x < xDim - xMargin; x++) {
				grid[gridPosition(x, y)]->computeDensity(runIndex);
			}
		}
	}

	
	//******************************************************************************************************************
	//******************************************************************************************************************
	// Print routines
	//******************************************************************************************************************

	void printGeometry() const{
		std::cout << std::endl;
		for (int y = 0; y < yDim; y++) {
			for (int x = 0; x < xDim; x++) {
				std::cout << geometry[(y * xDim) + x];
			}
			std::cout << std::endl;
		}
	}

	//******************************************************************************************************************
	// Cell type - print routine
	void printCellType() const{
		std::cout << std::endl;
		for (int y = 0; y < yDim; y++) {
			for (int x = 0; x < xDim; x++) {
				if (grid[gridPosition(x, y)] != nullptr) {
					std::cout << grid[(y * xDim) + x]->getCellTypeChar();
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
	void printNeighboursCellType(const int x, const int y) const{
		std::shared_ptr<Cell> tempCell;
		char E, NE, N, NW, W, SW, S, SE, R;


		if ((tempCell = grid[gridPosition(x, y)]->getCellNeighbours().getNeighbour(CellDirection::east)) != nullptr) {
			E = tempCell->getCellTypeChar();
		}
		else {
			E = 'N';
		}
		if ((tempCell = grid[gridPosition(x, y)]->getCellNeighbours().getNeighbour(CellDirection::northEast)) != nullptr) {
			NE = tempCell->getCellTypeChar();
		}
		else {
			NE = 'N';
		}
		if ((tempCell = grid[gridPosition(x, y)]->getCellNeighbours().getNeighbour(CellDirection::north)) != nullptr) {
			N = tempCell->getCellTypeChar();
		}
		else {
			N = 'N';
		}
		if ((tempCell = grid[gridPosition(x, y)]->getCellNeighbours().getNeighbour(CellDirection::northWest)) != nullptr) {
			NW = tempCell->getCellTypeChar();
		}
		else {
			NW = 'N';
		}
		if ((tempCell = grid[gridPosition(x, y)]->getCellNeighbours().getNeighbour(CellDirection::west)) != nullptr) {
			W = tempCell->getCellTypeChar();
		}
		else {
			W = 'N';
		}
		if ((tempCell = grid[gridPosition(x, y)]->getCellNeighbours().getNeighbour(CellDirection::southWest)) != nullptr) {
			SW = tempCell->getCellTypeChar();
		}
		else {
			SW = 'N';
		}
		if ((tempCell = grid[gridPosition(x, y)]->getCellNeighbours().getNeighbour(CellDirection::south)) != nullptr) {
			S = tempCell->getCellTypeChar();
		}
		else {
			S = 'N';
		}
		if ((tempCell = grid[gridPosition(x, y)]->getCellNeighbours().getNeighbour(CellDirection::southEast)) != nullptr) {
			SE = tempCell->getCellTypeChar();
		}
		else {
			SE = 'N';
		}

		// Rest shows cell at (x, y)
		if ((tempCell = grid[gridPosition(x, y)]) != nullptr) {
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
		const int xMargin = 1;
		const int yMargin = 1;
		for (int y = yMargin; y < yDim - yMargin; y++) {
			for (int x = xMargin; x < xDim - xMargin; x++) {
				printNeighboursCellType(x, y);
			}
		}
	}


	//******************************************************************************************************************
	// Cell Population - print routines

	void printCellPopulation(const bool runIndex, const int x, const int y) const{
		std::cout << std::endl;
		std::cout << grid[gridPosition(x, y)]->getPolulation(runIndex, CellDirection::northWest)
			<< " " << grid[gridPosition(x, y)]->getPolulation(runIndex, CellDirection::north)
			<< " " << grid[gridPosition(x, y)]->getPolulation(runIndex, CellDirection::northEast) << std::endl;

		std::cout << " \\|/ " << std::endl;

		std::cout << grid[gridPosition(x, y)]->getPolulation(runIndex, CellDirection::west)
			<< "-" << grid[gridPosition(x, y)]->getPolulation(runIndex, CellDirection::rest)
			<< "-" << grid[gridPosition(x, y)]->getPolulation(runIndex, CellDirection::east) << std::endl;

		std::cout << " /|\\ " << std::endl;

		std::cout << grid[gridPosition(x, y)]->getPolulation(runIndex, CellDirection::southWest)
			<< " " << grid[gridPosition(x, y)]->getPolulation(runIndex, CellDirection::south)
			<< " " << grid[gridPosition(x, y)]->getPolulation(runIndex, CellDirection::southEast) << std::endl;
	}

	void  printCellPopulation(const bool runIndex) const{
		const int xMargin = 1;
		const int yMargin = 1;
		for (int y = yMargin; y < yDim - yMargin; y++) {
			for (int x = xMargin; x < xDim - xMargin; x++) {
				printCellPopulation(runIndex, x, y);
			}
		}
	}


	//******************************************************************************************************************
	// Cell rho - print routine
	void printCellRho(const bool runIndex) const{
		/*std::cout << std::endl;*/
		for (int y = 0; y < yDim; y++) {
			for (int x = 0; x < xDim; x++) {
				if (grid[gridPosition(x, y)] != nullptr) {
					std::cout << grid[gridPosition(x, y)]->getDensity();
					std::cout << " - ";
				}
			}
			std::cout << std::endl;
		}
	}


	//******************************************************************************************************************
	// Cell velocity - print routine
	void printCellVelocity(const bool runIndex) const{
		std::cout << std::endl;
		for (int y = 0; y < yDim; y++) {
			for (int x = 0; x < xDim; x++) {
				if (grid[gridPosition(x, y)] != nullptr) {
					std::cout << "(" << grid[gridPosition(x, y)]->getVelocity(SpatialDirection::x) 
						<< "," 
						<< grid[gridPosition(x, y)]->getVelocity(SpatialDirection::y) << ")";
					std::cout << " - ";
				}
			}
			std::cout << std::endl;
		}
	}
	   

	std::string appendGridPolulationsList(const bool runIndex, std::string& populationLists) const {
		//std::string populationLists;
		populationLists += ((populationLists == "") ? "{" : ",\n\n{");
		for (int y = 0; y < yDim; y++) {
			populationLists += "{";
			for (int x = 0; x < xDim; x++) {
				//if (grid[gridPosition(x, y)] != nullptr) {
				populationLists += grid[gridPosition(x, y)]->getPopulationsList(runIndex) + ((x < xDim - 1) ? ",\n" : "");
				//}
			}

			populationLists += ((y < yDim - 1) ? "},\n\n" : "}");
		}
		populationLists += "}";	

		return populationLists;
		
	}
		
	
	void appendGridVelocityList(const bool runIndex, std::string& velocityLists) const {
		std::ostringstream velocityStringStream;
		velocityStringStream << std::setprecision(3);
		velocityStringStream << velocityLists;
		if (velocityStringStream.str() == "") {
			velocityStringStream << "{" << xDim << "," << yDim << "," << F_TAU << "," << std::setprecision(12) << std::fixed << F_BODY_FORCE_X << "," << F_BODY_FORCE_Y
				<< std::setprecision(3) << std::defaultfloat << "},\n\n" << "{";
		}
		else {
			velocityStringStream << ",\n\n{";
		}
		/*velocityStringStream << velocityLists << ((velocityLists == "") ? "{" + std::to_string(xDim) + "," + std::to_string(yDim) + "},\n\n" + "{" : ",\n\n{");*/
		for (int y = 0; y < yDim; y++) {
			velocityStringStream << "{";
			for (int x = 0; x < xDim; x++) {
				//if (grid[gridPosition(x, y)] != nullptr) {
				velocityStringStream << grid[gridPosition(x, y)]->getVelocityList() << ((x < xDim - 1) ? ",\n" : "");
				//}
			}

			velocityStringStream << ((y < yDim - 1) ? "},\n\n" : "}");
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
			densityStringStream << "{" << xDim << "," << yDim << "," << F_TAU << "," << std::setprecision(12) << std::fixed << F_BODY_FORCE_X << "," << F_BODY_FORCE_Y
				<< std::setprecision(3) << std::defaultfloat << "},\n\n" << "{";
		}
		else {
			densityStringStream << ",\n\n{";
		}
		/*densityStringStream << densityLists << ((densityLists == "") ? "{" + std::to_string(xDim) + "," + std::to_string(yDim) + "},\n\n" + "{" : ",\n\n{");*/
		for (int y = 0; y < yDim; y++) {
			densityStringStream << "{";
			for (int x = 0; x < xDim; x++) {
				//if (grid[gridPosition(x, y)] != nullptr) {
				densityStringStream << grid[gridPosition(x, y)]->getDensity() << ((x < xDim - 1) ? ",\n" : "");
				//}
			}

			densityStringStream << ((y < yDim - 1) ? "},\n\n" : "}");
		}
		densityStringStream << "}";

		densityLists = densityStringStream.str();

		//return densityLists;

	}

	std::shared_ptr<Cell> getCell(const int x, const int y) const{
		return grid[gridPosition(x, y)];
	}	
};