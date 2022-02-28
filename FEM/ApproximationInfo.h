#pragma once

#include <vector>

class ApproximationInfo
{
	using Matrix = std::vector<std::vector<double>>;
public:

	ApproximationInfo() = default;
	ApproximationInfo(const std::vector<Matrix>& stfMatrices, const Matrix& msMatrices, double stfCoeff, double msCoeff, const std::vector<double>& nmCoeffs);
	~ApproximationInfo() {}

	void AddLocalStiffnessMatrix(const Matrix& stiffnessMatrix) { m_LocalStiffnessMatrices.push_back(stiffnessMatrix); }
	void SetLocalMassMatrix(const Matrix& massMatrix) { m_LocalMassMatrix = massMatrix; }
	auto& GetLocalStiffnessMatrices() { return m_LocalStiffnessMatrices; }
	auto& GetLocalMassMatrix() { return m_LocalMassMatrix; }

	void SetStiffnessCoefficient(double coeff) { m_StiffnessCoefficient = coeff; }
	auto GetStiffnessCoefficient() { return m_StiffnessCoefficient; }
	auto GetMassCoefficient() { return m_MassCoefficient; }
	auto& GetNeumannCoefficients() { return m_NeumannCoefficients; }

private:
	std::vector<Matrix> m_LocalStiffnessMatrices;
	Matrix m_LocalMassMatrix;

	double m_StiffnessCoefficient = 1.0;
	double m_MassCoefficient = 1.0;
	std::vector<double> m_NeumannCoefficients;

};

namespace ApprInfo
{
	using Matrix = std::vector<std::vector<double>>;
	const ApproximationInfo LinearInfo(
		std::vector<Matrix>
		{
			{
				{  1.0, -1.0, 0.0 },
				{ -1.0,  1.0, 0.0 },
				{  0.0,  0.0, 0.0 }
			},
			{
				{  1.0, 0.0, -1.0 },
				{  0.0, 0.0,  0.0 },
				{ -1.0, 0.0,  1.0 }
			},
		},
		Matrix
		{
			{ 2.0, 1.0, 1.0 },
			{ 1.0, 2.0, 1.0 },
			{ 1.0, 1.0, 2.0 }
		},
		1.0 / 2.0, 1.0 / 24.0, { 3.0 / 6.0, 3.0 / 6.0 }
	);
}


