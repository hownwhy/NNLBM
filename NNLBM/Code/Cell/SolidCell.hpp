#pragma once
#include "../Globals.hpp"
#include "Cell.hpp"
#include <iostream>
#include <algorithm>


class SolidCell : public Cell
{
private:
	

public:
	
	~SolidCell() = default;
	SolidCell() = default;

	
	//*****************************************************************************************
	// Do functions

	//void addMovingBoundaryTerm(const bool runIndex) override {
	//	std::cout << "\nSolidCell::addMovingBoundaryTerm()";
	//	system("pause");
	//}
		
	void collide(const bool runIndex) {
		std::cout << "\nSolidCell::collide()";
		system("pause");
	}

	

	void collideAndPropagate(const bool runIndex) override{
		std::cout << "\nSolidCell::collideAndProppagate()";
		system("pause");
	}
	
	////Half way bounce back
	//void setReceived(const bool runIndex, const int populationIndex, const field_t fieldValue) {
	//	std::cout << "\nSolidCell::setReceived()";
	//	system("pause");
	//}

#if BOOK_BOUNCE_BACK
	void propageteTo(const bool runIndex) override{
		std::cout << "\nSolidCell::propagateTo()";
		system("pause");
	}
#else
	//Half way bounce back
	void setReceived(const bool runIndex, const int populationIndex, const field_t fieldValue) override {
		int arrayIndex = getArrayIndex(runIndex, reverseDirectionIndex(populationIndex));
		////std::cout << "setPopulation ARRAY_INDEX : " << arrayIndex << std::endl;;
		//assert(("setPopulations: arrayIndex is negative", arrayIndex >= 0));
		//assert(("setPopulations: arrayIndex to high", arrayIndex < nFieldDuplicates * nPopulations));
		populations_.at(arrayIndex) = fieldValue;
	}
#endif
	char getCellTypeChar() const override {
		return 'S';
	}
	
	CellType getCellType() const override{
		return CellType::solid;
	}
};