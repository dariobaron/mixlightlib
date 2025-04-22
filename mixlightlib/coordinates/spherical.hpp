#pragma once

#include <vector>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <numeric>
#include "../concepts.hpp"


template<unsigned DIM>
class Cartesian;


std::vector<double> convertSpherStandardOrder(std::vector<double> other){
	std::vector<double> result(other.size());
	result[0] = other[0];
	std::copy(other.rbegin(), other.rend()-1, result.begin()+1);
	return result;
}


template<unsigned DIM>
class Spherical{

private:
	std::vector<double> x_;

public:
	friend Cartesian<DIM>;
	friend std::vector<double> convertSpherStandardOrder(std::vector<double> other);

	// Constructors
	Spherical() : x_(DIM, 0) {};

	Spherical(const std::vector<double> & x) : x_(convertSpherStandardOrder(x)) {
		checkDimensionality();
	};
	
	Spherical(std::vector<double> && x) : x_(convertSpherStandardOrder(x)) {
		checkDimensionality();
	};
	
	Spherical(const Cartesian<DIM> & cart) : x_(DIM) {
		auto & cartesians = cart.x_;
		double r2 = 0;
		for (auto i : cartesians){
			r2 += i * i;
		}
		double r = std::sqrt(r2);
		x_[0] = r;
		if (r > 1e-10){
			double component = r;
			for (unsigned i = 1; i < DIM; ++i){
				x_[i] = std::acos(cartesians[i-1] / component);
				component *= std::sin(x_[i]);
			}
		}
	};

	// assignment operators
	Spherical<DIM>& operator=(Spherical<DIM> other){
		this->swap(other);
		return *this;
	};

	// output stream operator
	template<unsigned D>
	friend std::ostream& operator<<(std::ostream & output, const Spherical<D> & obj);

	// additive operators
	template<unsigned D>
	friend Spherical<D> operator+(Spherical<D> a, const Spherical<D> & b);

	template<unsigned D>
	friend Spherical<D> operator-(Spherical<D> a, const Spherical<D> & b);

	Spherical<DIM>& operator+=(const Spherical<DIM> & other){
		Cartesian<DIM> tmp(*this);
		tmp += other;
		*this = tmp;
		return *this;
	};

	Spherical<DIM>& operator-=(const Spherical<DIM> & other){
		Cartesian<DIM> tmp(*this);
		tmp -= other;
		*this = tmp;
		return *this;
	};

	// scalar product
	template<unsigned D>
	friend double operator*(const Spherical<D> & a, const Spherical<D> & b);
	
	// rescaling
	template<Numeric T>
	Spherical<DIM>& operator*=(T coeff){
		x_[0] *= coeff;
		return *this;
	};
	
	template<unsigned D,Numeric T>
	friend Spherical<D> operator*(Spherical<D> a, T b);

	template<unsigned D,Numeric T>
	friend Spherical<D> operator*(T a, Spherical<D> b);

	// norm
	double norm() const{
		return x_[0];
	};

private:
	void checkDimensionality() const{
		if (x_.size() != DIM){
			std::runtime_error("The coordinate provided does not match the dimensionality of the system");
		}
	}
};


template<unsigned DIM>
std::ostream& operator<<(std::ostream & output, const Spherical<DIM> & obj){
	output << "[";
	auto components = convertSpherStandardOrder(obj.x_);
	output << components[0];
	for (unsigned i = 1; i < DIM; ++i){
		output << "," << components[i];
	}
	output << "]\n";
	return output;
}

template<unsigned DIM>
Spherical<DIM> operator+(Spherical<DIM> a, const Spherical<DIM> & b){
	return Spherical(Cartesian<DIM>(a) + Cartesian<DIM>(b));
}

template<unsigned DIM>
Spherical<DIM> operator-(Spherical<DIM> a, const Spherical<DIM> & b){
	return Spherical(Cartesian<DIM>(a) - Cartesian<DIM>(b));
}

template<unsigned DIM>
double operator*(const Spherical<DIM> & a, const Spherical<DIM> & b){
	return Cartesian<DIM>(a) * Cartesian<DIM>(b);
}

template<unsigned DIM,Numeric T>
Spherical<DIM> operator*(Spherical<DIM> a, T b){
	return a *= b;
}

template<unsigned DIM,Numeric T>
Spherical<DIM> operator*(T a, Spherical<DIM> b){
	return b *= a;
}