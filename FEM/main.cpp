#include <algorithm>
#include <iostream>
#include <functional>
#include <array>
#include <utility>
#include <iomanip>
#include <fstream>
#include <cassert>

#include "output.h"
#include "DataStructures.h"
#include "functions.h"
#include "SparseMatrix.h"
#include "GridData.h"

#include "FEMSolver.h"

double testFunction(double x, double y) { return x; }
// ��� � ��������� P
double rhsFunction(double x, double y) { return 0; }

int main()
{
	FEMSolver FemSolver("FEMMode.txt"); // BL - ����������, BQ - ����������������, L - ��������, Q - ��������������.
	//FemSolver.LoadPlastPressure("input/Plast.txt");


	GridData gridData;
	gridData.BoundaryFunction = [&](double x, double y) { return testFunction(x, y); };
	gridData.BoundaryFunction = [&](double x, double y) { return FemSolver.GetPlastPressure(); };
	gridData.RhsFunction = [](double x, double y) { return rhsFunction(x, y); };
	gridData.LoadNodes("input/grid.txt");
	gridData.LoadDirichletConditions("input/BC1.txt");
	gridData.LoadElements("input/outTriangle.txt");
	gridData.LoadNeumannConditions("input/BC2_Tri.txt");

	FemSolver.SetGridData(gridData);
	FemSolver.LoadDomainData("input/mat.txt");
	FemSolver.LoadViscosity("input/phaseprop.txt");
	FemSolver.CalculateStiffnesMatrix();
	FemSolver.CalculateLoadVector();
	FemSolver.ApplyNeumann();
	FemSolver.ApplyDirichlet();
	FemSolver.CalculateSolution();
	auto solution = FemSolver.GetSolution();


	std::ofstream resout("result.txt");

	// ���� �� ������� 4 �������� (��������� �������), �� ������ ��������� �����.
	prettyPrint(std::cout, solution, gridData, testFunction);
	prettyPrint(resout, solution, gridData);
}

