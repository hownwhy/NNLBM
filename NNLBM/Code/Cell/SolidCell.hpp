#pragma once
#include "../Globals.hpp"
#include "Cell.hpp"
#include <iostream>
#include <algorithm>


class SolidCell : public Cell
{
private:

	field_t applyMovingSolid(field_t &population, const CellDirection cellDirection) {
		// Add velocity from moving walls
		const field_t csSqr = 1. / 3;
		// ! Excludes rest population in the following definitions of c and w. !
		const std::array<int, nCellDirections * nDimensions> c = { 1,1,0,-1,-1,-1,0,1, 0,1,1,1,0,-1,-1,-1 };
		const std::array<field_t, nCellDirections> w = { (1. / 9), (1. / 36), (1. / 9), (1. / 36), (1. / 9), (1. / 36), (1. / 9), (1. / 36) };

		auto getLatticeVelocity = [](const int spatialDirection, const int cellDirection) { return (spatialDirection * nCellDirections) + cellDirection; };
		return population - (2 * w.at(cellDirection) * density_ / csSqr)
			* ((c.at(getLatticeVelocity(SpatialDirection::x, cellDirection)) * velocity_[SpatialDirection::x])
				+ (c.at(getLatticeVelocity(SpatialDirection::y, cellDirection)) * velocity_[SpatialDirection::y]));
	}

public:

	~SolidCell() = default;
	SolidCell() = default;


	//*****************************************************************************************
	// Do functions

	void collideAndPropagate() override {
		std::cout << "\nSolidCell::collideAndProppagate()";
		system("pause");
	}

	/*void setReceived(std::array<field_t, nFieldDuplicates * nPopulations> &sourcePopulations, const CellDirection &cellDirection) override {		
		sourcePopulations[nextPopulationIndexOffset_ + reverseDirectionIndex(cellDirection)] = applyMovingSolid(sourcePopulations[currentPopulationIndexOffset_ + cellDirection], cellDirection);
	}*/

	void setReceived(field_t *sourcePopulations, const CellDirection &cellDirection) override {
		sourcePopulations[nextPopulationIndexOffset_ + reverseDirectionIndex(cellDirection)] 
			= applyMovingSolid(sourcePopulations[currentPopulationIndexOffset_ + cellDirection], cellDirection);
	}

	char getCellTypeChar() const override {
		return 'S';
	}

	CellType getCellType() const override {
		return CellType::solid;
	}
};
