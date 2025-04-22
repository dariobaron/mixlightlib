#pragma once

#include <vector>
#include <iostream>
#include <stdexcept>
#include "coordinates.hpp"


template<unsigned NROW,unsigned NCOL=NROW>
class LinearMap{

private:
	std::vector<double> M_;

public:

	// Constructors
	LinearMap() : M_(NROW*NCOL) {};

	LinearMap(double val) : M_(NROW*NCOL, val) {};

	LinearMap(const std::vector<double> & M) : M_(convertCartStandardOrder(M)) {
		checkDimensionality();
	};

	LinearMap(const std::vector<std::vector<double>> & M){
		for (auto & v : M){
			for (auto i : v){
				M_.push_back(i);
			}
		}
		checkDimensionality();
	};
	
	LinearMap(std::vector<double> && M) : M_(convertCartStandardOrder(M)) {
		checkDimensionality();
	};
	
	LinearMap(std::vector<std::vector<double>> && M){
		for (auto & v : M){
			for (auto i : v){
				M_.push_back(i);
			}
		}
		checkDimensionality();
	};

	// named constructors
	static LinearMap<NROW,NROW> eye(){
		LinearMap<NROW,NROW> m(0);
		for (unsigned i = 0; i < NROW; ++i){
			m(i,i) = 1;
		}
		return m;
	};

	// assignment operators
	LinearMap<NROW,NCOL>& operator=(LinearMap<NROW,NCOL> other){
		this->swap(other);
		return *this;
	};

	// output stream operator
	template<unsigned D,unsigned SD>
	friend std::ostream& operator<<(std::ostream & output, const LinearMap<D,SD> & obj);

	// additive operators
	template<unsigned D,unsigned SD>
	friend LinearMap<D,SD> operator+(LinearMap<D,SD> a, const LinearMap<D,SD> & b);

	template<unsigned D,unsigned SD>
	friend LinearMap<D,SD> operator-(LinearMap<D,SD> a, const LinearMap<D,SD> & b);

	LinearMap<NROW,NCOL>& operator+=(const LinearMap<NROW,NCOL> & other){
		for (unsigned i = 0; i < NROW; ++i){
			M_[i] += other.M_[i];
		}
		return *this;
	};

	LinearMap<NROW,NCOL>& operator-=(const LinearMap<NROW,NCOL> & other){
		for (unsigned i = 0; i < NROW; ++i){
			M_[i] -= other.M_[i];
		}
		return *this;
	};
	
	// rescaling
	template<Numeric T>
	LinearMap<NROW,NCOL>& operator*=(T coeff){
		for (auto & i : M_){
			i *= coeff;
		}
		return *this;
	};
	
	template<unsigned D,unsigned SD,Numeric T>
	friend LinearMap<D,SD> operator*(LinearMap<D,SD> a, T b);

	template<unsigned D,unsigned SD,Numeric T>
	friend LinearMap<D,SD> operator*(T a, LinearMap<D,SD> b);

	// matrix multiplication
	template<unsigned INNERDIM>
	friend LinearMap<NROW,NCOL> operator*(LinearMap<NROW,INNERDIM> a, const LinearMap<INNERDIM,NCOL> & b);

	template<unsigned NR,unsigned NC>
	friend Cartesian<NR> operator*(const LinearMap<NR,NC> & M, const Cartesian<NC> & v);

	template<unsigned NR,unsigned NC>
	friend Cartesian<NC> operator*(const Cartesian<NR> & v, const LinearMap<NR,NC> & M);

	// getters
	std::pair<unsigned,unsigned> shape() const{	return std::make_pair(NROW,NCOL);	};

	const double * data() const{	return M_.data();	};

	double operator()(unsigned r, unsigned c) const{	return M_[c+r*NCOL];	};

	std::vector<double> row(unsigned r) const{
		if (r >= NROW){	throw std::runtime_error("Requested row "+std::to_string(r)+" from matrix with "+std::to_string(NROW)+" rows!");	}
		return std::vector<double>(M_.begin()+r*NCOL, M_.begin()+(r+1)*NCOL);
	};

	std::vector<double> col(unsigned c) const{
		if (c >= NCOL){	throw std::runtime_error("Requested column "+std::to_string(c)+" from matrix with "+std::to_string(NCOL)+" columns!");	}
		std::vector<double> thecolumn(NROW);
		for (unsigned i = 0; i < NROW; ++i){
			thecolumn[i] = M_[i*NCOL+c];
		}
		return thecolumn;
	};

	// setters
	double & operator()(unsigned r, unsigned c){	return M_[c+r*NCOL];	};

	void row(unsigned r, const std::vector<double> & values){
		if (r >= NROW){	throw std::runtime_error("Requested row "+std::to_string(r)+" from matrix with "+std::to_string(NROW)+" rows!");	}
		else if (values.size() != NCOL){	throw std::runtime_error("Requested row assignment with vector of wrong length!");	}
		for (unsigned i = 0; i < NCOL; ++i){
			M_[r*NCOL+i] = values[i];
		}
	};

	void col(unsigned c, const std::vector<double> & values){
		if (c >= NCOL){	throw std::runtime_error("Requested column "+std::to_string(c)+" from matrix with "+std::to_string(NCOL)+" columns!");	}
		else if (values.size() != NROW){	throw std::runtime_error("Requested column assignment with vector of wrong length!");	}
		for (unsigned i = 0; i < NROW; ++i){
			M_[i*NCOL+c] = values[i];
		}
	};

	// transposition
	LinearMap<NCOL,NROW> transpose() const{
		LinearMap<NCOL,NROW> transposed;
		for (unsigned i = 0; i < NROW; ++i){
			for (unsigned j = 0; j < NCOL; ++j){
				transposed(j,i) = (*this)(i,j);
			}
		}
		return transposed;
	};

private:
	void checkDimensionality(){
		if (M_.size() != NROW*NCOL){
			std::runtime_error("The coordinate provided does not match the dimensionality of the system");
		}
	};
};


template<unsigned D,unsigned SD>
std::ostream& operator<<(std::ostream & output, const LinearMap<D,SD> & obj){
	auto entries = convertCartStandardOrder(obj.M_);
	output << "[ ";
	for (unsigned i = 0; i < D; ++i){
		output << "[" << entries[i*SD];
		for (unsigned j = 1; j < SD; ++j){	output << "," << entries[i*SD+j];	}
		output << "]";
		if (i != D-1) {	output << ",\n";	}
	}
	output << " ]\n";
	return output;
}

template<unsigned NROW,unsigned NCOL>
LinearMap<NROW,NCOL> operator+(LinearMap<NROW,NCOL> a, const LinearMap<NROW,NCOL> & b){
	return a += b;
}

template<unsigned NROW,unsigned NCOL>
LinearMap<NROW,NCOL> operator-(LinearMap<NROW,NCOL> a, const LinearMap<NROW,NCOL> & b){
	return a -= b;
}

template<unsigned NROW,unsigned NCOL,Numeric T>
LinearMap<NROW,NCOL> operator*(LinearMap<NROW,NCOL> a, T b){
	return a *= b;
}

template<unsigned NROW,unsigned NCOL,Numeric T>
LinearMap<NROW,NCOL> operator*(T a, LinearMap<NROW,NCOL> b){
	return b *= a;
}

template<unsigned NROW,unsigned INNERDIM,unsigned NCOL>
LinearMap<NROW,NCOL> operator*(LinearMap<NROW,INNERDIM> a, const LinearMap<INNERDIM,NCOL> & b){
	LinearMap<NROW,NCOL> result(0);
	for (unsigned i = 0; i < NROW; ++i){
		for (unsigned k = 0; k < INNERDIM; ++k){
			for (unsigned j = 0; j < NCOL; ++j){
				result(i,j) += a(i,k) * b(k,j);
			}
		}
	}
	return result;
}

template<unsigned NROW,unsigned NCOL>
Cartesian<NROW> operator*(const LinearMap<NROW,NCOL> & M, const Cartesian<NCOL> & v){
	std::vector<double> result(NROW,0);
	const double * v_ptr = v.data();
	const double * M_ptr = M.data();
	for (unsigned i = 0; i < NROW; ++i){
		for (unsigned j = 0; j < NCOL; ++j){
			result[NROW-1-i] += M_ptr[i*NCOL+j] * v_ptr[j];
		}
	}
	return Cartesian<NROW>(result);
}

template<unsigned NROW,unsigned NCOL>
Cartesian<NCOL> operator*(const Cartesian<NROW> & v, const LinearMap<NROW,NCOL> & M){
	return M.transpose() * v;
}