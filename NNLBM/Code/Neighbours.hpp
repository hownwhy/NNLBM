#pragma once
#include "Globals.hpp"
#include "Cell/Cell.hpp"
#include <array>
#include <memory>
#include <iostream>

class Cell;

class Neighbours {

public:
	~Neighbours() = default;
	Neighbours() = default;

private:
	static const int nNeighbours = nDirections;
	std::array<std::shared_ptr<Cell>, nNeighbours> neighbours;

public:
	void setNeighbour(const int direction, const std::shared_ptr<Cell> cell) {		
		//std::cout << "setNeighbour ARRAY_INDEX : " << direction << std::endl;
		assert(("setNeighbour: arrayIndex is negative", direction >= 0));
		assert(("setNeighbour: arrayIndex to high", direction < nNeighbours));
		neighbours.at(direction) = cell;
	}	

	std::shared_ptr<Cell> getNeighbour(const int direction) const {
		//std::cout << "getNeighbour ARRAY_INDEX : " << direction << std::endl;
		assert(("getNeighbour: arrayIndex is negative", direction >= 0));
		assert(("getNeighbour: arrayIndex to high", direction < nNeighbours));
		return neighbours.at(direction);
	}

	/*template<typename CellType>
	void initializeNeighbour(const int direction) {
		neighbours.at(direction) = std::make_shared<CellType>();
	}*/
};