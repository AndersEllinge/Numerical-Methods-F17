#include <Numerical code\nr3.h>
#include <Numerical code\quadrature.h>
#include <math.h>

Doub Y = 0.1;

template<class T>
struct MyTrapzd : Quadrature {
	Doub a, b, s;
	T &func;
	MyTrapzd() {};
	MyTrapzd(T &funcc, const Doub aa, const Doub bb) :
		func(funcc), a(aa), b(bb) {
		n = 0;
	}
	Doub next() {
		Doub x, tnm, sum, del;
		Int it, j;
		n++;
		if (n == 1) {
			//std::cout << "Number of traps: " << 1 << std::endl;
			return (s = 0.5*(b - a)*(func(a) + func(b)));
		}
		else {
			for (it = 1, j = 1;j<n - 1;j++) it <<= 1;
			tnm = it;
			//std::cout << "Number of traps: " << tnm*2 << std::endl;
			del = (b - a) / tnm;
			x = a + 0.5*del;
			for (sum = 0.0, j = 0;j<it;j++, x += del) sum += func(x);
			s = 0.5*(s + (b - a)*sum / tnm);
			return s;
		}
	}
};

template<class T>
Doub Myqtrap(T &func, const Doub a, const Doub b, const Doub eps = 1.0e-10) {
	const Int JMAX = 20;
	Doub s, olds = 0.0;
	MyTrapzd<T> t(func, a, b);
	for (Int j = 0;j<JMAX;j++) {
		s = t.next();
		//std::cout << "I ran again with j = " << j << " and s = " << s << std::endl;
		if (j > 5)
			if (abs(s - olds) < eps*abs(olds) ||
				(s == 0.0 && olds == 0.0)) return s;
		olds = s;
	}
	throw("Too many steps in routine qtrap");
}

Doub f(Doub gamma) {
	if (abs(gamma) < (Doub)1/sqrt((Doub)2))
		return 0.0;

	return (1 - (pow(1 - gamma, 2) / (((Doub)3 / (Doub)2) - sqrt(2)))) * cos(1.2 * pow(1 - gamma, 2));
}

Doub integralFunc(Doub t) {
	Doub bigGamma = pow((pow(1 * t, 2) / 0.0025) + (pow(Y, 2) / 0.04) + 1, -((Doub)1 / (Doub)2));
	//std::cout << "bigGamma: " << bigGamma << std::endl;

	return 0.744 * pow(bigGamma, 3) * f(bigGamma);
}

Doub c(Doub inputY) {
	if (abs(inputY) > (Doub)0.2)
		return (Doub)0.0;

	//std::cout << "I made it past catch" << std::endl;

	Y = inputY;

	Doub alpha = (Doub)0.05 * sqrt((Doub)1 - (pow(inputY, (Doub)2) / (Doub)0.04));
	std::cout << "alpha: " << alpha << std::endl;

	Doub res = Myqtrap(integralFunc, -alpha, alpha, 0.001);
	std::cout << "qtrap: " << res << std::endl;

	res = ((Doub)1 / (Doub)1200) * res;
	//std::cout << "res: " << res << std::endl;

	return res;
}

template<class T>
Doub estimateIntegral(int N, Doub inputY, T func) {
	Doub alpha = (Doub)0.05 * sqrt((Doub)1 - (pow(inputY, (Doub)2) / (Doub)0.04));
	MyTrapzd<T> tmp(func, -alpha, alpha);
	Doub res = 0;

	Y = inputY;

	for (int i = 0; i < N; i++) {
		res = tmp.next();
	}
	return res;
}

Doub calcOrder(int startN, Doub inputY) {
	Doub h1 = estimateIntegral(startN, inputY, integralFunc);
	Doub h2 = estimateIntegral(startN+1, inputY, integralFunc);
	Doub h3 = estimateIntegral(startN+2, inputY, integralFunc);

	return sqrt((h1 - h2) / (h2 - h3));
}

Doub calcRE(Doub h1, Doub h2) {
	return h2 + ((h2 - h1) / 3);
}


int main() {
	std::cout << "Calculations with trap:" << std::endl;
	std::cout << "res with trap: " << c(0.1) * 1000000 << " micrometer" << std::endl;
	std::cout << std::endl;

	std::cout << "Calculations with RE:" << std::endl;
	std::cout << "RE: " << calcRE(estimateIntegral(3, 0.1, integralFunc), estimateIntegral(4, 0.1, integralFunc)) << std::endl;
	std::cout << "res with RE: " << ((Doub)1 / (Doub)1200) * calcRE(estimateIntegral(3, 0.1, integralFunc), estimateIntegral(4, 0.1, integralFunc)) * 1000000 << " micrometer" << std::endl;
	std::cout << std::endl;

	std::cout << "Convergence of integral from trap and RE:" << std::endl;
	std::cout << "N=1   Trap: " << estimateIntegral(1, (Doub)0.1, integralFunc) << std::endl;
	std::cout << "N=2   Trap: " << estimateIntegral(2, (Doub)0.1, integralFunc) << "	RE: " << calcRE(estimateIntegral(1, 0.1, integralFunc), estimateIntegral(2, 0.1, integralFunc)) << std::endl;
	std::cout << "N=3   Trap: " << estimateIntegral(3, (Doub)0.1, integralFunc) << "	RE: " << calcRE(estimateIntegral(2, 0.1, integralFunc), estimateIntegral(3, 0.1, integralFunc)) << std::endl;
	std::cout << "N=4   Trap: " << estimateIntegral(4, (Doub)0.1, integralFunc) << "	RE: " << calcRE(estimateIntegral(3, 0.1, integralFunc), estimateIntegral(4, 0.1, integralFunc)) << std::endl;
	std::cout << "N=5   Trap: " << estimateIntegral(5, (Doub)0.1, integralFunc) << "	RE: " << calcRE(estimateIntegral(4, 0.1, integralFunc), estimateIntegral(5, 0.1, integralFunc)) << std::endl;
	std::cout << "N=6   Trap: " << estimateIntegral(6, (Doub)0.1, integralFunc) << "	RE: " << calcRE(estimateIntegral(5, 0.1, integralFunc), estimateIntegral(6, 0.1, integralFunc)) << std::endl;
	std::cout << std::endl;

	std::cout << "Convergence of order:" << std::endl;
	std::cout << "order estimate: " << calcOrder(1, 0.1) << std::endl;
	std::cout << "order estimate: " << calcOrder(2, 0.1) << std::endl;
	std::cout << "order estimate: " << calcOrder(3, 0.1) << std::endl;
	std::cout << std::endl;




	return 1;
}