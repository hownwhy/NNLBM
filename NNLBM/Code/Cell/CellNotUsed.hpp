//#pragma once
//
////Returns the 1D array index which depend on the 2D(runIndex, direction).
//inline uint_t getArrayIndex(const bool runIndex, const uint_t direction) const {
//	return (runIndex * nCellDirections) + direction;
//}
//
///* This function gives the opposite direction of what you put in.
// This is also the reason why I chose to not follow the convension
// having the rest direction be the 0 direction. 	*/
//const uint_t threeHalfDirection = (3 * nStreamDirections) / 2;
//inline uint_t reverseDirectionIndex(const uint_t direction) const {
//	return ((direction + threeHalfDirection) % nStreamDirections);
//}
//
//
//void initializeDensity(const field_t density) {
//	density_ = density;
//	computePopulationsEq();
//	std::copy(populationsEq_.begin(), populationsEq_.end(), populations_.begin());
//}
//
//void initializeVelocity(const SpatialDirection direction, const field_t velocity) {
//	velocity_.at(direction) = velocity; // -((bodyForce.at(direction) * dt) / (density_));
//	computePopulationsEq();
//	std::copy(populationsEq_.begin(), populationsEq_.end(), populations_.begin());
//}
//
//
//void setReceived(field_t *sourcePopulations, const uint_t cellDirection) {
//	if (cellType_ == CellType::bulk) {
//		populations_[nextPopulationIndexes_[cellDirection]] = sourcePopulations[currentPopulationIndexes_[cellDirection]];
//	}
//	else {
//		sourcePopulations[nextPopulationIndexes_[reverseDirectionIndex[cellDirection]]]
//			= solidToBulkVelocityTransfer(sourcePopulations[currentPopulationIndexes_[cellDirection]], cellDirection);
//	}
//}

// The code below is to be used ogether with anoter method to include body forces
	//void computeForcePopulations() {	
		/*// TODO: Why do I have to divide by 2 in "tauFactor"???
		const field_t tauFactor = (1. - (dt / (2 * tau))) / 2;
		const field_t weightHV = tauFactor / 3;
		const field_t weightD = tauFactor / 12;
		const field_t weightR = 4 * tauFactor / 9;
		field_t Fx = bodyForce.at(SpatialDirection::x);
		field_t Fy = bodyForce.at(SpatialDirection::y);
		field_t ux = velocity_.at(SpatialDirection::x);
		field_t uy = velocity_.at(SpatialDirection::y);
		field_t FxUx = Fx * ux;
		field_t FxUy = Fx * uy;
		field_t FyUy = Fy * uy;
		field_t FyUx = Fy * ux;

		// Calculate horizontal and vertical force components
		forcePopulations_.at(CellDirection::east) = weightHV * (2 * Fx + 2 * FxUx - FyUy);
		forcePopulations_.at(CellDirection::north) = weightHV * (2 * Fy + 2 * FyUy - FxUx);
		forcePopulations_.at(CellDirection::west) = weightHV * (-2 * Fx + 2 * FxUx - FyUy);
		forcePopulations_.at(CellDirection::south) = weightHV * (-2 * Fy + 2 * FyUy - FxUx);

		// Calculate diagonal force components
		forcePopulations_.at(CellDirection::northEast) = weightD * (2 * Fx + 2 * Fy + 2 * FxUx + 3 * FyUx + 3 * FxUy + 2 * FyUy);
		forcePopulations_.at(CellDirection::northWest) = -weightD * (2 * Fx - 2 * Fy - 2 * FxUx + 3 * FyUx + 3 * FxUy - 2 * FyUy);
		forcePopulations_.at(CellDirection::southWest) = weightD * (-2 * Fx - 2 * Fy + 2 * FxUx + 3 * FyUx + 3 * FxUy + 2 * FyUy);
		forcePopulations_.at(CellDirection::southEast) = -weightD * (-2 * Fx + 2 * Fy - 2 * FxUx + 3 * FyUx + 3 * FxUy - 2 * FyUy);

		// Calculate rest force component
		forcePopulations_.at(CellDirection::rest) = -weightR * (3 * FxUx + 3 * FyUy);

		//for (uint_t cellDirection = 0; cellDirection < nCellDirections; cellDirection++) {
		//	std::cout << "\n forcePopulations_[" << cellDirection << "] = " << forcePopulations_.at(cellDirection);
		//}
		//std::cout << "\n\n";
		//system("pause");
	}*/


	// Some other code related to the alternative body force method above
	/*field_t dotProduct(const std::array<field_t, nDimensions> leftVector, const std::array<field_t, nDimensions> rightVector) const{
		field_t result = 0.;
		for (uint_t spatialDirection = 0; spatialDirection < nDimensions; spatialDirection++) {
			result += leftVector.at(spatialDirection) * rightVector.at(spatialDirection);
		}
		return result;

		}

	void computeForcePopulations() {

		const field_t csSqr = 1. / 3;
		//const std::array<int, nCellDirections * nDimensions> c = { 1,1,0,-1,-1,-1,0,1,0, 0,1,1,1,0,-1,-1,-1,0 };
		const std::array<field_t, nDimensions * nCellDirections> c = { 1,0, 1,1, 0,1, -1,1, -1,0, -1,-1, 0,-1, 1,-1, 0,0 };
		const std::array<field_t, nCellDirections> w = { (1. / 9), (1. / 36), (1. / 9), (1. / 36), (1. / 9), (1. / 36), (1. / 9), (1. / 36), (4. / 9) };

		field_t cDotF = 0;
		field_t cDotu = 0;
		field_t uDotF = 0;

		for (uint_t cellDirection = 0; cellDirection < nCellDirections; cellDirection++) {
			std::array<field_t, nDimensions> cTemp = { c.at((2 * cellDirection) + 0) , c.at((2 * cellDirection) + 1) };
			cDotF = dotProduct(cTemp, bodyForce);
			cDotu = dotProduct(cTemp, velocity_);
			uDotF = dotProduct(velocity_, bodyForce);

			//std::cout << "\ncDotF = " << cDotF;
			//std::cout << "\ncDotu: ("
			//	<< c.at((2 * cellDirection) + 0) << " * " << velocity_.at(0) << ") + ("
			//	<< c.at((2 * cellDirection) + 1) << " * " << velocity_.at(1)
			//	<< ") = " << cDotu;
			////std::cout << "\tuDotF = " << uDotF;
			//std::cout << "\nuDotF: ("
			//	<< velocity_.at(0) << " * " << bodyForce.at(0) << ") + ("
			//	<< velocity_.at(1) << " * " << bodyForce.at(1)
			//	<< ") = " << uDotF;
			//std::cout << "\n popultionIndex = " << cellDirection;

			forcePopulations_.at(cellDirection) = w.at(cellDirection) * (1 - (dt/(2 * tau))) * ((cDotF / csSqr) + (((cDotF * cDotu) - (csSqr * uDotF)) / (csSqr * csSqr)));
			std::cout << "\n forcePopulations_[" << cellDirection << "] = " << forcePopulations_.at(cellDirection);
		}
		//std::cout << "\n forcePopulations_[" << 8 << "] = " << forcePopulations_.at(8);
		std::cout << "\n\n";
		system("pause");
	}*/