#include "stdafx.h"
#include "Jacobi.h"


int main(int argc, char * argv[])
{
	try {
		int procRank, procNum;
		MPI_Status status;
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &procNum);
		MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
		Jacobi jacobi = Jacobi(argc,argv);
		jacobi.solve();
		jacobi.saveToFile(argv[4]);
		MPI_Finalize();
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
	//system("pause");
    return 0;
}
