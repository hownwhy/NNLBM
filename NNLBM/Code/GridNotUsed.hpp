#pragma once


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
//				grid_.at(gridPosition(x, y))->initialize(runIndex, rho, xVelocity, yVelocity);
//			}
//		}
//	}
//}

//// C++17 parallel for loop
		//std::for_each(
		//	//std::execution::seq,
		//	//std::execution::par,
		//	std::execution::par_unseq,
		//	grid_.begin(),
		//	grid_.end(),
		//	[](auto &cell) {cell.collideAndPropagate(); }
		//);

		//for (uint_t y = yMargin; y < yDimTotal - yMargin; y++) {
		//	//std::cout << std::endl;
		//	for (uint_t x = xMargin; x < xDimTotal - xMargin; x++) {
		//		//std::cout << "(" << x << "," << y << ")";
		//		grid_.at(gridPosition(x, y)).collideAndPropagate();
		//	}
		//}
//	}

	//void moveBoundary(const bool runIndex) const {
	//	const uint_t xMargin = N_X_SOLID_CELLS;
	//	const uint_t yMargin = N_Y_SOLID_CELLS;
	//	for (uint_t y = yMargin; y < yDimTotal - yMargin; y++) {
	//		for (uint_t x = xMargin; x < xDimTotal - xMargin; x++) {
	//			if (y == 0 + yMargin) {
	//				grid_.at(gridPosition(x, y))->addMovingBoundaryTerm(runIndex);
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
	//			grid_.at(gridPosition(x, y))->computeDensity(runIndex);
	//		}
	//	}
	//}


//******************************************************************************************************************
	// Neighbours - print routines
	//void printNeighboursCellType(const uint_t x, const uint_t y) const {
	//	std::shared_ptr<Cell> tempCell;
	//	char E, NE, N, NW, W, SW, S, SE, R;


	//	if ((tempCell = grid_.at(gridPosition(x, y))->getNeighbour(CellDirection::east)) != nullptr) {
	//		E = tempCell->getCellTypeChar();
	//	}
	//	else {
	//		E = 'N';
	//	}
	//	if ((tempCell = grid_.at(gridPosition(x, y))->getNeighbour(CellDirection::northEast)) != nullptr) {
	//		NE = tempCell->getCellTypeChar();
	//	}
	//	else {
	//		NE = 'N';
	//	}
	//	if ((tempCell = grid_.at(gridPosition(x, y))->getNeighbour(CellDirection::north)) != nullptr) {
	//		N = tempCell->getCellTypeChar();
	//	}
	//	else {
	//		N = 'N';
	//	}
	//	if ((tempCell = grid_.at(gridPosition(x, y))->getNeighbour(CellDirection::northWest)) != nullptr) {
	//		NW = tempCell->getCellTypeChar();
	//	}
	//	else {
	//		NW = 'N';
	//	}
	//	if ((tempCell = grid_.at(gridPosition(x, y))->getNeighbour(CellDirection::west)) != nullptr) {
	//		W = tempCell->getCellTypeChar();
	//	}
	//	else {
	//		W = 'N';
	//	}
	//	if ((tempCell = grid_.at(gridPosition(x, y))->getNeighbour(CellDirection::southWest)) != nullptr) {
	//		SW = tempCell->getCellTypeChar();
	//	}
	//	else {
	//		SW = 'N';
	//	}
	//	if ((tempCell = grid_.at(gridPosition(x, y))->getNeighbour(CellDirection::south)) != nullptr) {
	//		S = tempCell->getCellTypeChar();
	//	}
	//	else {
	//		S = 'N';
	//	}
	//	if ((tempCell = grid_.at(gridPosition(x, y))->getNeighbour(CellDirection::southEast)) != nullptr) {
	//		SE = tempCell->getCellTypeChar();
	//	}
	//	else {
	//		SE = 'N';
	//	}

	//	// Rest shows cell at (x, y)
	//	if ((tempCell = grid_.at(gridPosition(x, y))) != nullptr) {
	//		R = tempCell->getCellTypeChar();
	//	}
	//	else {
	//		R = 'N';
	//	}

	//	std::cout << std::endl;
	//	std::cout << NW
	//		<< " " << N
	//		<< " " << NE << std::endl;

	//	std::cout << " \\|/ " << std::endl;

	//	std::cout << W
	//		<< "-" << R
	//		<< "-" << E << std::endl;

	//	std::cout << " /|\\ " << std::endl;

	//	std::cout << SW
	//		<< " " << S
	//		<< " " << SE << std::endl;
	//}

	/*void printNeighboursCellType() const {
		const uint_t xMargin = N_X_SOLID_CELLS;
		const uint_t yMargin = N_Y_SOLID_CELLS;
		for (uint_t y = yMargin; y < yDimTotal - yMargin; y++) {
			for (uint_t x = xMargin; x < xDimTotal - xMargin; x++) {
				printNeighboursCellType(x, y);
			}
		}
	}*/


	//******************************************************************************************************************
	// Cell Population - print routines

	/*void printCellPopulation(const bool runIndex, const uint_t x, const uint_t y) const {
		std::cout << std::endl;
		std::cout << grid_.at(gridPosition(x, y))->getPopulation(CellDirection::northWest)
			<< " " << grid_.at(gridPosition(x, y))->getPopulation(CellDirection::north)
			<< " " << grid_.at(gridPosition(x, y))->getPopulation(CellDirection::northEast) << std::endl;

		std::cout << " \\|/ " << std::endl;

		std::cout << grid_.at(gridPosition(x, y))->getPopulation(CellDirection::west)
			<< "-" << grid_.at(gridPosition(x, y))->getPopulation(CellDirection::rest)
			<< "-" << grid_.at(gridPosition(x, y))->getPopulation(CellDirection::east) << std::endl;

		std::cout << " /|\\ " << std::endl;

		std::cout << grid_.at(gridPosition(x, y))->getPopulation(CellDirection::southWest)
			<< " " << grid_.at(gridPosition(x, y))->getPopulation(CellDirection::south)
			<< " " << grid_.at(gridPosition(x, y))->getPopulation(CellDirection::southEast) << std::endl;
	}*/

	/*void  printCellPopulation(const bool runIndex) const {
		const uint_t xMargin = N_X_SOLID_CELLS;
		const uint_t yMargin = N_Y_SOLID_CELLS;
		for (uint_t y = yMargin; y < yDimTotal - yMargin; y++) {
			for (uint_t x = xMargin; x < xDimTotal - xMargin; x++) {
				printCellPopulation(runIndex, x, y);
			}
		}
	}*/


	//******************************************************************************************************************
	// Cell rho - print routine
	//void printCellRho(const bool runIndex) const {
	//	/*std::cout << std::endl;*/
	//	for (uint_t y = 0; y < yDimTotal; y++) {
	//		for (uint_t x = 0; x < xDimTotal; x++) {
	//			if (grid_.at(gridPosition(x, y)) != nullptr) {
	//				std::cout << grid_.at(gridPosition(x, y))->getDensity();
	//				std::cout << " - ";
	//			}
	//		}
	//		std::cout << std::endl;
	//	}
	//}


	//******************************************************************************************************************
	// Cell velocity - print routine
	/*void printCellVelocity(const bool runIndex) const {
		std::cout << std::endl;
		for (uint_t y = 0; y < yDimTotal; y++) {
			for (uint_t x = 0; x < xDimTotal; x++) {
				if (grid_.at(gridPosition(x, y)) != nullptr) {
					std::cout << "(" << grid_.at(gridPosition(x, y))->getVelocity(SpatialDirection::x)
						<< ","
						<< grid_.at(gridPosition(x, y))->getVelocity(SpatialDirection::y) << ")";
					std::cout << " - ";
				}
			}
			std::cout << std::endl;
		}
	}*/



	//std::string appendGridNonEqPoplulationsList(const bool runIndex, std::string& populationLists)  {
	//	//std::string populationLists;
	//	populationLists += ((populationLists == "") ? "{" : ",\n\n{");
	//	for (uint_t y = 0; y < yDimTotal; y++) {
	//		populationLists += "{";
	//		for (uint_t x = 0; x < xDimTotal; x++) {
	//			populationLists += grid_.at(gridPosition(x, y)).getNonEqPopulationsList(runIndex) + ((x < xDimTotal - 1) ? ",\n" : "");
	//		}
	//		populationLists += ((y < yDimTotal - 1) ? "},\n\n" : "}");
	//	}
	//	populationLists += "}";

	//	return populationLists;
	//}