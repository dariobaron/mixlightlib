#pragma once

#include <list>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include "gausskronrod.hpp"
#include "mapper.hpp"


template <typename Func, typename T=double>
class Integrator{
public:
	using Solver = GaussKronrod<Func,T>;
	using FuncWrapper = Mapper<Func,T>;
private:
	const Func & f_;
	double a_, b_;
	double tolerance_;
	double scale_, loc_;
	double result_;
	double error_;
	unsigned maxiter_;
	unsigned nSubintervals_;
	bool success_;
public:
	Integrator(const Func & f, double a, double b, double tolerance, double scale=1, double loc=0, unsigned maxiter=100);
	double a() const;
	double b() const;
	double tolerance() const;
	double result() const;
	double error() const;
	unsigned nSubintervals() const;
	bool success() const;
private:
	bool reverse_;
	void integrate();
};


template <typename Func, typename T>
Integrator<Func,T>::Integrator(const Func & f, double a, double b, double tolerance, double scale, double loc, unsigned maxiter)
												: f_(f), a_(std::min(a,b)), b_(std::max(a,b)), tolerance_(tolerance), scale_(scale), loc_(loc),
												result_(0), error_(0), maxiter_(maxiter), nSubintervals_(0), success_(false), reverse_(a>b)
{
	integrate();
	if (reverse_){
		result_ = -result_;
	}
}


template <typename Func, typename T>
double Integrator<Func,T>::a() const{
	return reverse_ ? b_ : a_;
}


template <typename Func, typename T>
double Integrator<Func,T>::b() const{
	return reverse_ ? a_ : b;
}


template <typename Func, typename T>
double Integrator<Func,T>::tolerance() const{
	return tolerance_;
}


template <typename Func, typename T>
double Integrator<Func,T>::result() const{
	return result_;
}


template <typename Func, typename T>
double Integrator<Func,T>::error() const{
	return error_;
}


template <typename Func, typename T>
unsigned Integrator<Func,T>::nSubintervals() const{
	return nSubintervals_;
}


template <typename Func, typename T>
bool Integrator<Func,T>::success() const{
	return success_;
}


template <typename Func, typename T>
void Integrator<Func,T>::integrate(){
	std::list<Solver> subintegrals;
	subintegrals.emplace_back(FuncWrapper(f_, a_, b_, scale_, loc_), tolerance_);
	auto current = subintegrals.begin();
	while (current != subintegrals.end()){
		current->integrate();
		if (current->success() && subintegrals.size() > 1){
			++current;
		}
		else{
			auto to_remove = current;
			double a = current->a();
			double b = current->b();
			double m = FuncWrapper(f_, a, b, scale_, loc_).transformVariable(0);
			double new_tolerance = 0.5 * current->tolerance();
			current = subintegrals.insert(current, {Solver(FuncWrapper(f_, a, m, scale_, loc_), new_tolerance), Solver(FuncWrapper(f_, m, b, scale_, loc_), new_tolerance)});
			subintegrals.erase(to_remove);
		}
		if (subintegrals.size() > maxiter_){
			while (current != subintegrals.end()){
				current->integrate();
				++current;
			}
			break;
		}
	}
	subintegrals.sort([](const Solver & left, const Solver & right){return left.result() < right.result();});
	result_ = std::accumulate(subintegrals.begin(), subintegrals.end(), 0.0, [](double sum, const Solver & integral){return sum + integral.result();});
	subintegrals.sort([](const Solver & left, const Solver & right){return left.error() < right.error();});
	error_ = std::accumulate(subintegrals.begin(), subintegrals.end(), 0.0, [](double sum, const Solver & integral){return sum + integral.error();});
	nSubintervals_ = subintegrals.size();
	success_ = std::all_of(subintegrals.begin(), subintegrals.end(), [](const Solver & integral){return integral.success();});
}