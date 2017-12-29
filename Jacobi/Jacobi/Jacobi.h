#pragma once
class Jacobi
{

private:
	double* coefficients = nullptr;
	double* freeMembers = nullptr;
	double* aproxValues = nullptr;
	double epsilon;
	int numberOfEquations;
	int factNumberOfEquations;
	//int newCountOfEquations;
	int procNum;
	int procRank;

public:
	void loadDataFromFiles(const std::string fnCoefficients,
									const std::string fnApproxValues,
									const double epsilon);
	void saveToFile(const std::string fnAnswer);
	void solve();
	Jacobi(int argc, char * argv[]);
	~Jacobi();
};

