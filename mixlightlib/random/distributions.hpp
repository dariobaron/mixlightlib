#pragma once

#include <random>


class GammaDistribution{
	std::gamma_distribution<> distr_;
	double x0;
public:
	GammaDistribution() : distr_(), x0(0) {};
	GammaDistribution(double k, double theta, double x0=0) : distr_(k, theta), x0(x0) {};
	template<typename RndEng>
	double operator()(RndEng & eng){
		return distr_(eng) - x0;
	}
};


class BetaDistribution{
	std::gamma_distribution<> distr_a_;
	std::gamma_distribution<> distr_b_;
	double x0_;
	double scale_;
public:
	BetaDistribution() : distr_a_(), distr_b_(), x0_(0), scale_(1) {};
	BetaDistribution(double alpha, double beta, double x0=0, double scale=1) : distr_a_(alpha, 1), distr_b_(beta, 1), x0_(x0), scale_(scale) {};
	template<typename RndEng>
	double operator()(RndEng & eng){
		double a = distr_a_(eng);
		double b = distr_b_(eng);
		return a / (a + b) * scale_ - x0_;
	}
};


template<typename T>
class TwoPointDistribution{
	std::bernoulli_distribution distr_;
	T val1_;
	T val2_;
public:
	TwoPointDistribution(double p=0.5, T val1=0, T val2=1) : distr_(p), val1_(val1), val2_(val2) {};
	template<typename RndEng>
	T operator()(RndEng & eng){
		if (distr_(eng))	{	return val1_;	}
		else				{	return val2_;	}
	}
};


template<typename T>
class DiscreteDistribution{
	std::discrete_distribution<unsigned> distr_;
	std::vector<T> values_;
public:
	DiscreteDistribution() = default;
	DiscreteDistribution(const std::vector<double> & probs) : distr_(probs.begin(), probs.end()), values_() {};
	DiscreteDistribution(const std::vector<double> & probs, const std::vector<T> & values) : distr_(probs.begin(), probs.end()), values_(values) {};
	DiscreteDistribution(std::initializer_list<double> probs) : distr_(probs), values_() {};
	DiscreteDistribution(std::initializer_list<double> probs, std::initializer_list<double> values) : distr_(probs), values_(values) {};
	template<typename RndEng>
	T operator()(RndEng & eng){
		unsigned idx = distr_(eng);
		if (values_.size())	{	return values_[idx];	}
		else				{	return idx;				}
	}
};
