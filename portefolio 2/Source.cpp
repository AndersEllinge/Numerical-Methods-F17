#include <Numerical code\nr3.h>
#include <Numerical code\quadrature.h>
#include <Numerical code\roots.h>
#include <Numerical code\mins.h>
#include <math.h>
#include <iostream>
#include <fstream>

Doub Y = 0.1;
Doub x;

template<class T>
Doub trapezoidalIntegral(T &func, const Doub a, const Doub b, int N) {
	Doub h = (b - a) / N;
	Doub res;

	res = 0.5 * (func(a) + func(b)); // Add ends to sum

	for (int i = 0; i < N; i++) { // Add inbetween parts
		res = res + func(a + h + (i*h));
	}

	res = res * h; // Multiply by h

	return res;
}

template<class T>
Doub RE(T &func, const Doub a, const Doub b, int N) {
	if (N < 1 || N % 2 == 1)
		return 1;

	Doub Ah1 = trapezoidalIntegral(func, a, b, (N / 2));
	Doub Ah2 = trapezoidalIntegral(func, a, b, N);

	return Ah2 + ((Ah2 - Ah1) / 3);
}

template<class T>
Doub estimateErrorTrap(T &func, const Doub a, const Doub b, int N) {
	if (N == 1) // Catch for bad N
		return 1;

	Doub h2 = (b - a) / N;
	Doub h1 = (b - a) / (N / 2);
	Doub Ah2 = trapezoidalIntegral(func, a, b, N);
	Doub Ah1 = trapezoidalIntegral(func, a, b, N / 2);
	Doub alpha = h1 / h2;


	Doub error = (Ah1 - Ah2) / (pow(alpha, 2) - 1);

	return error;
}

template<class T>
Doub estimateErrorRE(T &func, const Doub a, const Doub b, int N) {
	if (N == 1) // Catch for bad N
		return 1;

	Doub Ah2 = RE(func, a, b, N);
	Doub Ah1 = RE(func, a, b, N / 2);

	Doub error = (Ah1 - Ah2) / 15;

	return error;
}

template<class T>
Doub calcOrder(T &func, const Doub a, const Doub b, int N) {
	Doub Ah1 = trapezoidalIntegral(func, a, b, N);
	Doub Ah2 = trapezoidalIntegral(func, a, b, N * 2);
	Doub Ah3 = trapezoidalIntegral(func, a, b, N * 4);

	return log2((Ah1 - Ah2) / (Ah2 - Ah3));
}

Doub f(Doub gamma) {
	if (abs(gamma) < (Doub)1/sqrt((Doub)2))
		return 0.0;

	return (1 - (pow(1 - gamma, 2) / (((Doub)3 / (Doub)2) - sqrt(2)))) * cos(1.2 * pow(1 - gamma, 2));
}

Doub integralFunc(Doub t) {
	Doub bigGamma = pow((pow(1 * t, 2) / 0.0025) + (pow(Y, 2) / 0.04) + 1, -((Doub)1 / (Doub)2));

	return 0.744 * pow(bigGamma, 3) * f(bigGamma);
}

Doub cTrap(Doub inputY, int N) {
	if (abs(inputY) > (Doub)0.2)
		return (Doub)0.0;

	Y = inputY;

	Doub alpha = (Doub)0.05 * sqrt((Doub)1 - (pow(inputY, (Doub)2) / (Doub)0.04));

	Doub res = trapezoidalIntegral(integralFunc, -alpha, alpha, N);
	std::cout << "Estimated integral: " << res << std::endl;
	std::cout << "Estimated error: " << estimateErrorTrap(integralFunc, -alpha, alpha, N) << std::endl;
	std::cout << "Estimated order: " << calcOrder(integralFunc, -alpha, alpha, N) << std::endl;

	res = ((Doub)1 / (Doub)1200) * res;
	return res;
}

Doub cTrapNoCout(Doub inputY, int N) {
	if (abs(inputY) > (Doub)0.2)
		return (Doub)0.0;

	Y = inputY;

	Doub alpha = (Doub)0.05 * sqrt((Doub)1 - (pow(inputY, (Doub)2) / (Doub)0.04));

	Doub res = trapezoidalIntegral(integralFunc, -alpha, alpha, N);

	res = ((Doub)1 / (Doub)1200) * res;
	return res;
}


Doub cRE(Doub inputY, int N) {
	if (abs(inputY) > (Doub)0.2)
		return (Doub)0.0;

	Y = inputY;

	Doub alpha = (Doub)0.05 * sqrt((Doub)1 - (pow(inputY, (Doub)2) / (Doub)0.04));
	
	Doub res = RE(integralFunc, -alpha, alpha, N);
	std::cout << "Estimated integral: " << res << std::endl;
	std::cout << "Estimated error: " << estimateErrorRE(integralFunc, -alpha, alpha, N) << std::endl;

	res = ((Doub)1 / (Doub)1200) * res;
	return res;
}

Doub calcAlpha(Doub inputY) {
	Doub alpha = (Doub)0.05 * sqrt((Doub)1 - (pow(inputY, (Doub)2) / (Doub)0.04));
	return alpha;
}

Doub minfunc(Doub Yi) {
	return cTrap(Yi, 100) + cTrap(x - Yi, 100);
}

Doub minifunc(Doub d) {

	// Brute min
	Doub acc = 0.00001;
	Doub min = cTrapNoCout(0, 100) + cTrapNoCout(d - 0, 100);
	Doub newMin;
	Doub Yi = 0;
	for (int i = 1; i < d / acc; i++) {
		Yi = Yi + acc;
		newMin = cTrapNoCout(Yi, 100) + cTrapNoCout(d - Yi, 100);
		if (newMin < min)
			min = newMin;
	}

	// To see the iterative root find values, out comment below
	//std::cout << "d: " << setprecision(15) << d << std::endl;
	//std::cout << "min: " << setprecision(15) << min * 1000000 << std::endl;
	//std::cout << "res: " << setprecision(15) << (min - (30.0 / 1000000)) * 1000000 << std::endl;
	//std::cout << std::endl;

	return min - (30.0 / 1000000);
}


int main() {

	for (int i = 4; i < 33; i = i * 2) {
		std::cout << "Calculation of c with trap, N = " << i << ": " << std::endl;
		std::cout << "res with trap: " << cTrap(0.1, i) * 1000000 << " micrometer" << std::endl;
		std::cout << std::endl;

		std::cout << "Calculation of c with RE, N = " << i << ": " << std::endl;
		std::cout << "res with RE: " << cRE(0.1, i) * 1000000 << " micrometer" << std::endl;
		std::cout << std::endl;
	}
	
	

	Doub Yi;
	int N;
	ofstream outputFile;
	std::string data3 = "dataFor3.csv";
	outputFile.open(data3);

	for (Doub i = 0; i < 50; i++) {
		Yi = (i / 50.0)*0.2;
		N = round((2 * calcAlpha(Yi)) / 0.001);
		
		outputFile << i << "," << Yi << "," << cTrapNoCout(Yi, N) * 1000000 << std::endl;
	}

	outputFile.close();

	std::string data4 = "dataFor4.csv";
	outputFile.open(data4);

	for (Doub i = 0; i < 50; i++) {
		Yi = (i / 50.0) * 0.2;
		N = round((2 * calcAlpha(Yi)) / 0.001);
		
		outputFile << i << "," << Yi << "," << ((cTrapNoCout(Yi, 100) + cTrapNoCout(0.234389134310186 - Yi, 100)) * 1000000) << std::endl;
	}

	outputFile.close();

	std::cout << "Calculating, this might take a minute..." << std::endl;
	std::cout << "The minimum distance d is: " << setprecision(15) << rtbis(minifunc, 0.2, 0.3, pow(10, -10)) << std::endl;
	


	return 1;
}