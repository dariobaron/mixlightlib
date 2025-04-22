#pragma once

#include <algorithm>
#include "gausskronrodtables.hpp"
#include "mapper.hpp"


template <typename Func, typename T=double>
class GaussKronrod{
private:
	Mapper<Func,T> f_;
	double tolerance_;
	double result_;
	double error_;
	bool success_;
public:
	GaussKronrod(const Mapper<Func,T> & f, double tolerance);
	void integrate();
	Mapper<Func,T> mapper() const;
	double a() const;
	double b() const;
	double tolerance() const;
	double result() const;
	double error() const;
	bool success() const;
};


template <typename Func, typename T>
GaussKronrod<Func,T>::GaussKronrod(const Mapper<Func,T> & f, double tolerance)
	: f_(f), tolerance_(tolerance), result_(0), error_(0), success_(false) {}


template <typename Func, typename T>
void GaussKronrod<Func,T>::integrate(){
	result_ = 0;
	double result7 = 0;
	for (unsigned i = 0; i < 15; ++i){
		double evaluation = f_(gctable::xgk15[i]);
		result_ += gctable::wgk15[i] * evaluation;
		if (i % 2){
			result7 += gctable::wgk7[i/2] * evaluation;
		}
	}
	error_ = std::abs(result7 - result_);
	success_ = error_ < tolerance_;
}


template <typename Func, typename T>
Mapper<Func,T> GaussKronrod<Func,T>::mapper() const{
	return f_;
}


template <typename Func, typename T>
double GaussKronrod<Func,T>::a() const{
	return f_.a();
}


template <typename Func, typename T>
double GaussKronrod<Func,T>::b() const{
	return f_.b();
}


template <typename Func, typename T>
double GaussKronrod<Func,T>::tolerance() const{
	return tolerance_;
}


template <typename Func, typename T>
double GaussKronrod<Func,T>::result() const{
	return result_;
}


template <typename Func, typename T>
double GaussKronrod<Func,T>::error() const{
	return error_;
}


template <typename Func, typename T>
bool GaussKronrod<Func,T>::success() const{
	return success_;
}