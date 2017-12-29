#include "stdafx.h"
#include "Jacobi.h"

#define root 0
#define tag_epsilon 0
#define tag_size 1


Jacobi::Jacobi(int argc, char * argv[])
{
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);

	if (procRank == root) {
		loadDataFromFiles(argv[1], argv[2], atof(argv[3]));
		for (int i = 0; i < procNum; i++) {
			if (i != root)
			{
				MPI_Send(&epsilon, 1, MPI_INT, i, tag_epsilon, MPI_COMM_WORLD);
				MPI_Send(&numberOfEquations, 1, MPI_INT, i, tag_size, MPI_COMM_WORLD);
			}
		}
	}
	else {
		MPI_Recv(&epsilon, 1, MPI_INT, root, tag_epsilon, MPI_COMM_WORLD, &status);
		MPI_Recv(&numberOfEquations, 1, MPI_INT, root, tag_size, MPI_COMM_WORLD, &status);
		aproxValues = new double[numberOfEquations];
	}
}

void Jacobi::saveToFile(const std::string fnAnswer) {
	if (procRank == root) {
		fstream outputFile;
		outputFile.open(fnAnswer, fstream::out);

		outputFile << factNumberOfEquations << endl;
		for (int i = 0; i < factNumberOfEquations; i++) {
			outputFile << aproxValues[i] << endl;
		}

		outputFile.close();

		outputFile.open("checkAnswer", fstream::out);

		outputFile << numberOfEquations << endl;
		double rightAnswer;
		for (int i = 0; i < numberOfEquations; i++) {
			rightAnswer = 0;
			for (int j = 0; j < numberOfEquations; j++) {
				rightAnswer += coefficients[i*numberOfEquations + j] * aproxValues[j];
			}

			outputFile << rightAnswer << " = " << freeMembers[i] << endl;
		}

		outputFile.close();
	}
}

void Jacobi::solve() {
	int sizeBuffer = numberOfEquations / procNum;
	double dtStart, dtEnd, dtDelta;
	double norm = 0, maxNorm = 0;
	double *tempX = new double[sizeBuffer];
	double* threadCoef = new double[sizeBuffer * numberOfEquations];
	double* threadFreeMembers = new double[sizeBuffer];
	//int procRank;

//	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

	if (procRank == root) {
		cout << "Size buffer " << sizeBuffer << endl;
	}
	MPI_Bcast(aproxValues, numberOfEquations, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Scatter(coefficients, sizeBuffer*numberOfEquations, MPI_DOUBLE, threadCoef, sizeBuffer * numberOfEquations, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Scatter(freeMembers, sizeBuffer, MPI_DOUBLE, threadFreeMembers, sizeBuffer, MPI_DOUBLE, root, MPI_COMM_WORLD);
	dtStart = MPI_Wtime();
	do {
		for (int i = 0; i < sizeBuffer; i++) {
			tempX[i] = threadFreeMembers[i];
			for (int j = 0; j < numberOfEquations; j++) {
				if ((i + procRank*sizeBuffer) != j) {
					tempX[i] -= threadCoef[i*numberOfEquations + j] * aproxValues[j];
				}
			}
			tempX[i] /= threadCoef[i*numberOfEquations + i + procRank*sizeBuffer];
		}

		norm = fabs(aproxValues[procRank*sizeBuffer] - tempX[0]);
		for (int i = 0; i < sizeBuffer; i++) {
			if (fabs(aproxValues[procRank*sizeBuffer + i] - tempX[i]) > norm) {
				norm = fabs(aproxValues[procRank*sizeBuffer + i] - tempX[i]);
			}
		}

		MPI_Reduce(&norm, &maxNorm, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
		MPI_Bcast(&maxNorm, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
		MPI_Allgather(tempX, sizeBuffer, MPI_DOUBLE, aproxValues, sizeBuffer, MPI_DOUBLE, MPI_COMM_WORLD);

		if (procRank == root) {
			cout << maxNorm << endl;
		}

	} while (maxNorm > epsilon);

	dtEnd = MPI_Wtime();

	dtDelta = dtEnd - dtStart;



	if (procRank == root) {
		cout << dtDelta;
	}

	delete[] threadCoef;
	delete[] threadFreeMembers;
	delete[] tempX;
}

void Jacobi::loadDataFromFiles(const std::string fnCoefficients, const std::string fnApproxValues, const double epsilon) {
	try
	{
		fstream inputFile;

		int countRowsCoef = 0;
		int countColumnsCoef = 0;
		int countRowsApproxValues = 0;

		inputFile.open(fnCoefficients);

		if (!inputFile.is_open()) {
			throw exception("Can't open file with coefficients ");
		}

		inputFile >> countRowsCoef;
		inputFile >> countColumnsCoef;

		factNumberOfEquations = countRowsCoef;

		if (countRowsCoef % procNum != 0) {
			numberOfEquations = (countRowsCoef / procNum + 1) * procNum;
		}
		else {
			numberOfEquations = countRowsCoef;
		}


		this->coefficients = new double[numberOfEquations*numberOfEquations];
		this->freeMembers = new double[numberOfEquations];

		for (int i = 0; i < numberOfEquations; i++) {
			for (int j = 0; j < numberOfEquations; j++) {
				if (i < countRowsCoef && j < countRowsCoef) {
					inputFile >> coefficients[i*numberOfEquations + j];
				}
				else {
					if (i < countRowsCoef || j < countRowsCoef) {
						coefficients[i*numberOfEquations + j] = 0;
					}
					else {
						coefficients[i*numberOfEquations + j] = 1;
					}
				}
			}

			if (i < countRowsCoef) {
				inputFile >> freeMembers[i];
			}
			else {
				freeMembers[i] = 1;
			}
		}
		inputFile.close();

		/*for (int i = countRowsCoef; i < numberOfEquations; i++) {
			for (int j = countRowsCoef; j < numberOfEquations; j++) {
				if (i == j) {
					coefficients[i*numberOfEquations + j] = 1.0;
				}
				else {
					coefficients[i*numberOfEquations + j] = 0.0;
				}
			}
			freeMembers[i] = 0.0;
		}
		*/

		inputFile.open(fnApproxValues);

		if (!inputFile.is_open()) {
			delete[] coefficients;
			delete[] freeMembers;
			throw exception("Can't open file with aproximation values");
		}

		inputFile >> countRowsApproxValues;

		if (countRowsApproxValues != countRowsApproxValues) {
			delete[] coefficients;
			delete[] freeMembers;
			throw exception("Incorrect files with start values");
		}

		this->aproxValues = new double[numberOfEquations];

		for (int i = 0; i < countRowsApproxValues; i++) {
			inputFile >> aproxValues[i];
		}

		for (int i = countRowsApproxValues; i < numberOfEquations; i++) {
			aproxValues[i] = 0.0;
		}

		//this->numberOfEquations = countRowsApproxValues;

		fstream outputFile;
		outputFile.open("newAnswer.txt", fstream::out);

		outputFile << numberOfEquations << endl;
		for (int i = 0; i < numberOfEquations; i++) {
			outputFile << aproxValues[i] << endl;
		}

		outputFile.close();

		outputFile.open("newCoef.txt", fstream::out);

		outputFile << numberOfEquations << endl;

		for (int i = 0; i < numberOfEquations; i++) {
			for (int j = 0; j < numberOfEquations; j++) {
				outputFile << coefficients[i*numberOfEquations + j] << " ";
			}
			outputFile << freeMembers[i] << endl;
		}

		outputFile.close();
	}
	catch (const std::runtime_error re) {
		cout << "Runtime error: " << re.what() << endl;

	}
	catch (exception ex) {
		cout << "Exception: " << ex.what() << endl;
	}

	catch (...) {
		cout << ("Unknown exception");
	}
}

Jacobi::~Jacobi()
{
	if (procRank == root) {

		delete[] coefficients;
		delete[] this->freeMembers;
	}

	delete[] aproxValues;
}
