#pragma once

#include <vector>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <numeric>
#include "../concepts.hpp"


template<unsigned DIM>
class Spherical;


std::vector<double> convertCartStandardOrder(std::vector<double> other){
	return std::vector<double>(other.rbegin(), other.rend());
}


template<unsigned DIM>
class Cartesian{

private:
	std::vector<double> x_;

public:
	friend Spherical<DIM>;
	friend std::vector<double> convertCartStandardOrder(std::vector<double> other);

	// Constructors
	Cartesian() : x_(DIM, 0) {};

	Cartesian(double val) : x_(DIM, val) {};

	Cartesian(const std::vector<double> & x) : x_(convertCartStandardOrder(x)) {
		checkDimensionality();
	};
	
	Cartesian(std::vector<double> && x) : x_(convertCartStandardOrder(x)) {
		checkDimensionality();
	};
	
	Cartesian(const Spherical<DIM> & sph) : x_(DIM) {
		auto & sphericals = sph.x_;
		for (unsigned i = 0; i < DIM; ++i){
			x_[i] = sphericals[0];
			for (unsigned j = 1; j <= i; ++j){
				x_[i] *= std::sin(sphericals[j]);
			}
			if (i != DIM - 1){
				x_[i] *= std::cos(sphericals[i+1]);
			}
		}
	};

	// assignment operators
	Cartesian<DIM>& operator=(Cartesian<DIM> other){
		this->swap(other);
		return *this;
	};

	// output stream operator
	template<unsigned D>
	friend std::ostream& operator<<(std::ostream & output, const Cartesian<D> & obj);

	// additive operators
	template<unsigned D>
	friend Cartesian<D> operator+(Cartesian<D> a, const Cartesian<D> & b);

	template<unsigned D>
	friend Cartesian<D> operator-(Cartesian<D> a, const Cartesian<D> & b);

	Cartesian<DIM>& operator+=(const Cartesian<DIM> & other){
		for (unsigned i = 0; i < DIM; ++i){
			x_[i] += other.x_[i];
		}
		return *this;
	};

	Cartesian<DIM>& operator-=(const Cartesian<DIM> & other){
		for (unsigned i = 0; i < DIM; ++i){
			x_[i] -= other.x_[i];
		}
		return *this;
	};

	// scalar product
	template<unsigned D>
	friend double operator*(const Cartesian<D> & a, const Cartesian<D> & b);
	
	// rescaling
	template<Numeric T>
	Cartesian<DIM>& operator*=(T coeff){
		for (auto & i : x_){
			i *= coeff;
		}
		return *this;
	};
	
	template<unsigned D,Numeric T>
	friend Cartesian<D> operator*(Cartesian<D> a, T b);

	template<unsigned D,Numeric T>
	friend Cartesian<D> operator*(T a, Cartesian<D> b);

	// norm
	double norm() const{
		return std::inner_product(x_.begin(), x_.end(), x_.begin(), 0.);
	};

	// low-level access
	const double * data() const{
		return x_.data();
	}

private:
	void checkDimensionality() const{
		if (x_.size() != DIM){
			std::runtime_error("The coordinate provided does not match the dimensionality of the system");
		}
	}
};


template<unsigned DIM>
std::ostream& operator<<(std::ostream & output, const Cartesian<DIM> & obj){
	output << "[";
	auto components = convertCartStandardOrder(obj.x_);
	output << components[0];
	for (unsigned i = 1; i < DIM; ++i){
		output << "," << components[i];
	}
	output << "]\n";
	return output;
}

template<unsigned DIM>
Cartesian<DIM> operator+(Cartesian<DIM> a, const Cartesian<DIM> & b){
	return a += b;
}

template<unsigned DIM>
Cartesian<DIM> operator-(Cartesian<DIM> a, const Cartesian<DIM> & b){
	return a -= b;
}

template<unsigned DIM>
double operator*(const Cartesian<DIM> & a, const Cartesian<DIM> & b){
	return std::inner_product(a.x_.begin(), a.x_.end(), b.x_.begin(), 0.);
}

template<unsigned DIM,Numeric T>
Cartesian<DIM> operator*(Cartesian<DIM> a, T b){
	return a *= b;
}

template<unsigned DIM,Numeric T>
Cartesian<DIM> operator*(T a, Cartesian<DIM> b){
	return b *= a;
}