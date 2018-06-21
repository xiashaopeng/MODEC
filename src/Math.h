#ifndef CMATH_H
#define CMATH_H

#include "Complex.h"
#include <math.h>

namespace CMath
{
	Complex exp_C(const Complex &rhs);

	Complex sin_C(const Complex &rhs);

	Complex cos_C(const Complex &rhs);

	Complex pow_C(const Complex &lhs, int pow);

	int Factorial(int N);
}

#endif