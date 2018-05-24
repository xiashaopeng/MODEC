#ifndef COMPLEX_H
#define COMPLEX_H
#include <iostream>
#include <cmath>

using namespace std;

class Complex {
  public:
    double real_;
    double imag_;

  public:
    Complex() {
        real_ = 0.0;
        imag_ = 0.0;
    }
    Complex(const double &r) {
        real_ = r;
        imag_ = 0.0;
    }
    Complex(const int &r) {
        real_ = r;
        imag_ = 0.0;
    }
    Complex(const double &r, const double &i) {
        real_ = r;
        imag_ = i;
    }
    Complex(const int &r, const double &i) {
        real_ = r;
        imag_ = i;
    }
    Complex(const double &r, const int &i) {
        real_ = r;
        imag_ = i;
    }
    Complex(const int &r, const int &i) {
        real_ = r;
        imag_ = i;
    }

    Complex(const Complex &number) {
        real_ = number.real_;
        imag_ = number.imag_;
    }

//-----------------------------------------运算符重载------------------------------------------//
    // 加减乘除运算符重载 //
    Complex operator+(const Complex &number) {
        return Complex(real_ + number.real_,imag_ + number.imag_);
    }
    Complex operator+(const double &number) {
        return Complex(real_ + number,imag_);
    }
    Complex operator+(const int &number) {
        return Complex(real_ + number,imag_);
    }

    Complex operator-(const Complex &number) {
        return Complex(real_ - number.real_,imag_ - number.imag_);
    }
    Complex operator-(const double &number) {
        return Complex(real_ - number,imag_);
    }
    Complex operator-(const int &number) {
        return Complex(real_ - number,imag_);
    }

    Complex operator*(const Complex &number) {
        return Complex(real_ * number.real_ - imag_ * number.imag_ , real_ * number.imag_ + imag_ * number.real_);
    }
    Complex operator*(const double &number) {
        return Complex(real_ * number , imag_ * number);
    }
    Complex operator*(const int &number) {
        return Complex(real_ * number , imag_ * number);
    }

    // In order to allow your class to appear as the right hand operand //
    friend Complex operator*(double lhs, const Complex &rhs) {
        return Complex(lhs * rhs.real_, lhs * rhs.imag_);
    }
    friend Complex operator*(int lhs, const Complex &rhs) {
        return Complex(lhs * rhs.real_, lhs * rhs.imag_);
    }

    friend Complex operator*(const Complex &lhs, const Complex &rhs) {
        return Complex(lhs.real_ * rhs.real_ - lhs.imag_ * rhs.imag_, lhs.real_ * rhs.imag_ + lhs.imag_ * rhs.real_);
    }
    //

    Complex operator/(const Complex &number) {
        double demominator = number.real_ * number.real_  + number.imag_  *  number.imag_;
        return Complex((real_ * number.real_ + imag_ * number.imag_)/demominator , (- real_ * number.imag_ + imag_ * number.real_)/demominator);
    }
    Complex operator/(const double &number) {
        return Complex(real_ / number , imag_ / number);
    }
    Complex operator/(const int &number) {
        return Complex(real_ / number , imag_ / number);
    }

	// In order to allow your class to appear as the right hand operand //
    friend Complex operator/(const double lhs, const Complex &rhs) {
        return Complex(lhs * rhs.real_  / (rhs.real_ * rhs.real_ + rhs.imag_ * rhs.imag_), - lhs * rhs.imag_ / (rhs.real_ * rhs.real_ + rhs.imag_ * rhs.imag_));
    }
    friend Complex operator/(const int lhs, const Complex &rhs) {
        return Complex(lhs * rhs.real_  / (rhs.real_ * rhs.real_ + rhs.imag_ * rhs.imag_), - lhs * rhs.imag_ / (rhs.real_ * rhs.real_ + rhs.imag_ * rhs.imag_));
    }

    friend Complex operator/(const Complex &lhs, const Complex &rhs) {
        return Complex((lhs.real_ * rhs.real_ + lhs.imag_ * rhs.imag_)/ (rhs.real_ * rhs.real_ + rhs.imag_ * rhs.imag_), ( - lhs.real_ * rhs.imag_ + lhs.imag_ * rhs.real_) / (rhs.real_ * rhs.real_ + rhs.imag_ * rhs.imag_));
    }
    //
	
	// In order to allow your class to appear as the right hand operand //
    friend Complex operator+(const double lhs, const Complex &rhs) {
        return Complex(lhs + rhs.real_ , rhs.imag_ );
    }
    friend Complex operator+(const int lhs, const Complex &rhs) {
        return Complex(lhs + rhs.real_ , rhs.imag_);
    }

    friend Complex operator+(const Complex &lhs, const Complex &rhs) {
        return Complex(lhs.real_ + rhs.real_ , lhs.imag_ + rhs.imag_);
    }
    //
	
	
    // ----------------- //

    // 自增减运算符重载 //
    Complex& operator++() { // 前置++
        ++ real_;
        return *this;
    }

    Complex operator++(int) { // 后置++
        Complex temp = *this;
        ++ *this;
        return temp;
    }

    Complex& operator--() { // 前置--
        -- real_;
        return *this;
    }

    Complex operator--(int) { // 后置--
        Complex temp = *this;
        -- *this;
        return temp;
    }
    // ----------------- //

    // += 和 -= 运算符重载 //
    Complex& operator+=(const Complex &number) {
        real_ += number.real_;
        imag_ += number.imag_;
        return *this;
    }
    Complex& operator+=(const double &number) {
        real_ += number;
        return *this;
    }
    Complex& operator+=(const int &number) {
        real_ += number;
        return *this;
    }

    Complex& operator-=(const Complex &number) {
        real_ -= number.real_;
        imag_ -= number.imag_;
        return *this;
    }
    Complex& operator-=(const double &number) {
        real_ -= number;
        return *this;
    }
    Complex& operator-=(const int &number) {
        real_ -= number;
        return *this;
    }

    Complex& operator*=(const Complex &number) {
        double temp_r = real_;
        double temp_i = imag_;
        real_ = temp_r * number.real_ - temp_i * number.imag_;
        imag_ = temp_r * number.imag_ + temp_i * number.real_;
        return *this;
    }
    Complex& operator*=(const double &number) {
        real_ *= number;
        imag_ *= number;
        return *this;
    }
    Complex& operator*=(const int &number) {
        real_ *= number;
        imag_ *= number;
        return *this;
    }

    Complex& operator=(const Complex &number) {
        if(this != &number) {
            real_ = number.real_;
            imag_ = number.imag_;
        }

        return *this;
    }
    Complex& operator=(const double &number) {
        real_ = number;
        imag_ = 0.0;
        return *this;
    }
    Complex& operator=(const int &number) {
        real_ = number;
        imag_ = 0.0;
        return *this;
    }
    //---------------------------------------//


    // == ， != ，逻辑运算符重载 //
    bool operator==(const Complex &number) {
        if((real_ == number.real_) && (imag_ == number.imag_)) {
            return true;
        } else {
            return false;
        }
    }
    bool operator==(const double &number) {
        if((real_ == number) && (imag_ == 0.0)) {
            return true;
        } else {
            return false;
        }
    }
    bool operator==(const int &number) {
        if((real_ == number) && (imag_ == 0.0)) {
            return true;
        } else {
            return false;
        }
    }

    bool operator!=(const Complex &number) {
        if((real_ == number.real_) && (imag_ == number.imag_)) {
            return false;
        } else {
            return true;
        }
    }
    bool operator!=(const double &number) {
        if((real_ == number) && (imag_ == 0.0)) {
            return false;
        } else {
            return true;
        }
    }
    bool operator!=(const int &number) {
        if((real_ == number) && (imag_ == 0.0)) {
            return false;
        } else {
            return true;
        }
    }
    // -------------------------------------------------------------------//
    //     cout, cin chongzai  //
    friend ostream& operator<<(ostream& os, const Complex& number) {
        os << "(" << number.real_ <<","<< number.imag_<<")";
        return os;
    }
};

#endif
