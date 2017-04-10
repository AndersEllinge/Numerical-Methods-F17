/* 
COMMENT BOX:



*/

#include <Numerical code\nr3.h>
#include <Numerical code\svd.h>
#include <Numerical code\utilities.h>

#include <math.h>

Doub norm(VecDoub input);
VecDoub var(MatDoub v, VecDoub w);

int main() {

	// Import data from dataset
	VecDoub theta1_d1(500), theta2_d1(500), measuredX_d1(500), measuredY_d1(500);
	std::ifstream d1("d1");

	for (int i = 0; i < 500; i++) {
		d1 >> theta1_d1[i];
		d1 >> theta2_d1[i];
		d1 >> measuredX_d1[i];
		d1 >> measuredY_d1[i];
	}

	VecDoub theta1_d2(500), theta2_d2(500), measuredX_d2(500), measuredY_d2(500);
	std::ifstream d2("d2");

	for (int i = 0; i < 500; i++) {
		d2 >> theta1_d2[i];
		d2 >> theta2_d2[i];
		d2 >> measuredX_d2[i];
		d2 >> measuredY_d2[i];
	}


	/// Create matrix A ///
	// d1
	MatDoub A_d1(1000, 4);
	for (int i = 0; i < 500; i++) {
		// Construct row for X value
		A_d1[i * 2][0] = 1;
		A_d1[i * 2][1] = 0;
		A_d1[i * 2][2] = cos(theta1_d1[i]);
		A_d1[i * 2][3] = cos(theta1_d1[i] + theta2_d1[i]);

		// Construct row for Y value
		A_d1[(i * 2) + 1][0] = 0;
		A_d1[(i * 2) + 1][1] = 1;
		A_d1[(i * 2) + 1][2] = sin(theta1_d1[i]);
		A_d1[(i * 2) + 1][3] = sin(theta1_d1[i] + theta2_d1[i]);

	}

	// d2
	MatDoub A_d2(1000, 4);
	for (int i = 0; i < 500; i++) {
		// Construct row for X value
		A_d2[i * 2][0] = 1;
		A_d2[i * 2][1] = 0;
		A_d2[i * 2][2] = cos(theta1_d2[i]);
		A_d2[i * 2][3] = cos(theta1_d2[i] + theta2_d2[i]);

		// Construct row for Y value
		A_d2[(i * 2) + 1][0] = 0;
		A_d2[(i * 2) + 1][1] = 1;
		A_d2[(i * 2) + 1][2] = sin(theta1_d2[i]);
		A_d2[(i * 2) + 1][3] = sin(theta1_d2[i] + theta2_d2[i]);

	}

	/// Create z vector ///
	// d1
	VecDoub z_d1(1000);
	for (int i = 0; i < 500; i++) {
		z_d1[i * 2] = measuredX_d1[i];
		z_d1[(i * 2) + 1] = measuredY_d1[i];
	}

	// d2
	VecDoub z_d2(1000);
	for (int i = 0; i < 500; i++) {
		z_d2[i * 2] = measuredX_d2[i];
		z_d2[(i * 2) + 1] = measuredY_d2[i];
	}

	/// Find solution with SVD ///
	// d1
	VecDoub q_d1(4);
	SVD svd_d1(A_d1);

	svd_d1.solve(z_d1, q_d1, svd_d1.eps);

	// d2
	VecDoub q_d2(4);
	SVD svd_d2(A_d2);

	svd_d2.solve(z_d2, q_d2, svd_d2.eps);

	/// Find the residuals ///
	// d1
	VecDoub residuals_d1(1000);
	residuals_d1 = A_d1*q_d1;
	for (int i = 0; i < residuals_d1.size(); i++)
		residuals_d1[i] = residuals_d1[i] - z_d1[i];

	// d2
	VecDoub residuals_d2(1000);
	residuals_d2 = A_d2*q_d2;
	for (int i = 0; i < residuals_d2.size(); i++)
		residuals_d2[i] = residuals_d2[i] - z_d2[i];

	/// Variance of parameters, eq(15.4.19) from NR3 ///
	// d1
	VecDoub var_d1(4);
	for (int j = 0; j < 4; j++) {
		var_d1[j] = pow(svd_d1.v[j][0] / svd_d1.w[0], 2);
		for (int i = 1; i < 4; i++)
			var_d1[j] = var_d1[j] + pow(svd_d1.v[j][i] / svd_d1.w[i], 2);
	}

	// d2
	VecDoub var_d2(4);
	for (int j = 0; j < 4; j++) {
		var_d2[j] = pow(svd_d2.v[j][0] / svd_d2.w[0], 2);
		for (int i = 1; i < 4; i++)
			var_d2[j] = var_d2[j] + pow(svd_d2.v[j][i] / svd_d2.w[i], 2);
	}

	/// Print results ///
	// Printing w
	std::cout << "w for d1:";
	util::print(svd_d1.w);
	std::cout << "w for d2:";
	util::print(svd_d2.w);
	std::cout << std::endl;

	// Printing v
	std::cout << "v for d1:\n";
	util::print(svd_d1.v);
	std::cout << "v for d2:\n";
	util::print(svd_d2.v);
	std::cout << std::endl;

	// Printing solution
	std::cout << "Solution for d1:";
	util::print(q_d1);
	std::cout << "Solution for d2:";
	util::print(q_d2);
	std::cout << std::endl;

	// Printing residual
	std::cout << "Residual for d1: " << norm(residuals_d1) << std::endl;
	std::cout << "Residual for d2: " << norm(residuals_d2) << std::endl;
	std::cout << std::endl;

	// Printing variance of x0
	std::cout << "Variance of x0 for d1: " << var_d1[0] << std::endl;
	std::cout << "Variance of x0 for d2: " << var_d2[0] << std::endl;
	std::cout << std::endl;

	// Printing variance of y0
	std::cout << "Variance of y0 for d1: " << var_d1[1] << std::endl;
	std::cout << "Variance of y0 for d2: " << var_d2[1] << std::endl;
	std::cout << std::endl;

	// Printing variance of a
	std::cout << "Variance of a for d1: " << var_d1[2] << std::endl;
	std::cout << "Variance of a for d2: " << var_d2[2] << std::endl;
	std::cout << std::endl;

	// Printing variance of b
	std::cout << "Variance of b for d1: " << var_d1[3] << std::endl;
	std::cout << "Variance of b for d2: " << var_d2[3] << std::endl;
	std::cout << std::endl;

	return 1;
}

Doub norm(VecDoub input) {
	Doub tempSum = 0;
	for (int i = 0; i < input.size(); i++)
		tempSum += pow(input[i], 2);
	tempSum = sqrt(tempSum);
	return tempSum;
}