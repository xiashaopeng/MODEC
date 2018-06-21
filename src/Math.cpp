#include "Math.h"

Complex CMath::exp_C(const Complex &rhs)
{
	double coeff = exp(rhs.real_);
	return Complex(coeff * cos(rhs.imag_), coeff * sin(rhs.imag_));
};

Complex CMath::sin_C(const Complex &rhs)
{
	return Complex(sin(rhs.real_) * (exp(rhs.imag_) + exp(-rhs.imag_)) / 2.,
		cos(rhs.real_) * (exp(rhs.imag_) - exp(-rhs.imag_)) / 2.);
}

Complex CMath::cos_C(const Complex &rhs)
{
	return Complex(cos(rhs.real_) * (exp(rhs.imag_) + exp(-rhs.imag_)) / 2.,
		sin(rhs.real_) * (exp(-rhs.imag_) - exp(rhs.imag_)) / 2.);
}

Complex CMath::pow_C(const Complex &lhs, int pow)
{
	if (pow == 0) return 1.0;
	else if (pow == 1) return lhs;
	else
		return lhs * pow_C(lhs, pow - 1);
};

int CMath::Factorial(int N)
{
	if (N > 0)
		return N * Factorial(N - 1);
	else
		return 1;
}