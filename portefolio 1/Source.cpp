/* 
COMMENT BOX:



*/

#include <Numerical code\nr3.h>
#include <Numerical code\svd.h>
#include <Numerical code\utilities.h>

#include <math.h>

int main() {

	// Import data from dataset
	VecDoub theta1(500), theta2(500), measuredX(500), measuredY(500);
	std::ifstream input("d1");

	for (int i = 0; i < 500; i++) {
		input >> theta1[i];
		input >> theta2[i];
		input >> measuredX[i];
		input >> measuredY[i];
	}


	// Create A matrix
	MatDoub A(1000, 4);
	for (int i = 0; i < 500; i++) {
		// Construct row for X value
		A[i * 2][0] = 1;
		A[i * 2][1] = 0;
		A[i * 2][2] = cos(theta1[i]);
		A[i * 2][3] = cos(theta1[i] + theta2[i]);

		// Construct row for Y value
		A[(i * 2) + 1][0] = 0;
		A[(i * 2) + 1][1] = 1;
		A[(i * 2) + 1][2] = sin(theta1[i]);
		A[(i * 2) + 1][3] = sin(theta1[i] + theta2[i]);

	}

	// Create z vector
	VecDoub z(1000);
	for (int i = 0; i < 500; i++) {
		z[i * 2] = measuredX[i];
		z[(i * 2) + 1] = measuredY[i];
	}

	// Find solution with SVD
	VecDoub q(4);
	SVD svd(A);

	svd.solve(z, q, svd.eps);
	
	util::print(q);

	return 1;
}