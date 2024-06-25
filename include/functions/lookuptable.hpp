#pragma once

#include <vector>
#include <exception>
#include <algorithm>
#include <limits>


class LookupTable{
	std::vector<double> xs_;
	std::vector<double> fs_;
public:
	LookupTable(const std::vector<double> & xs, const std::vector<double> & fs);
	double evaluate(double x) const;
	double interpolate(double x, unsigned a, unsigned b) const;
};


LookupTable::LookupTable(const std::vector<double> & xs, const std::vector<double> & fs) : xs_(xs), fs_(fs) {
	if (xs_.size() != fs_.size()){
		throw std::runtime_error("Constructor arrays have different size");
	}
	for (unsigned i = 0; i < xs_.size(); ++i){
		if ((i > 0) && (xs_[i] <= xs_[i-1])){
			throw std::runtime_error("Domain vector is unsorted or has repeated values");
		}
	}
}


double LookupTable::evaluate(double x) const{
	auto where = std::lower_bound(xs_.begin(), xs_.end(), x);
	unsigned idx = std::distance(xs_.begin(), where);
	if		(where == xs_.end())						{	return std::numeric_limits<double>::signaling_NaN();	}
	else if	((where == xs_.begin()) && (x != xs_[0]))	{	return std::numeric_limits<double>::signaling_NaN();	}
	else if	(*where == x)								{	return fs_[idx];	}
	else												{	return interpolate(x, idx-1, idx);	}
}


double LookupTable::interpolate(double x, unsigned a, unsigned b) const{
	double p = (x - xs_[a]) / (xs_[b] - xs_[a]);
	return p * fs_[b] + (1-p) * fs_[a];
}
