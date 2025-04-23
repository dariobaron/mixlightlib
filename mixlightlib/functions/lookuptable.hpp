#pragma once

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <limits>


class LookupTable{
	struct Element{
		double x;
		double f;
		Element() = default;
		Element(double x, double f) : x(x), f(f) {};
		bool operator<(const Element & other) const{ return x < other.x; };
	};
	std::vector<Element> elements_;
public:
	LookupTable() = default;
	LookupTable(const std::vector<double> & xs, const std::vector<double> & fs);
	double a() const;
	double b() const;
	double evaluate(double x) const;
	double operator()(double x) const;
	void extend(const LookupTable & other);
private:
	double interpolate(double x, unsigned a, unsigned b) const;

};


LookupTable::LookupTable(const std::vector<double> & xs, const std::vector<double> & fs) : elements_(xs.size()) {
	if (xs.size() != fs.size()){
		throw std::runtime_error("Constructor arrays have different size");
	}
	elements_[0].x = xs[0];
	elements_[0].f = fs[0];
	for (unsigned i = 1; i < elements_.size(); ++i){
		if ((xs[i] <= xs[i-1])){
			throw std::runtime_error("Domain vector is unsorted or has repeated values");
		}
		elements_[i].x = xs[i];
		elements_[i].f = fs[i];
	}
}


double LookupTable::a() const{
	return elements_.front().x;
}
double LookupTable::b() const{
	return elements_.back().x;
}


double LookupTable::evaluate(double x) const{
	auto where = std::lower_bound(elements_.begin(), elements_.end(), Element(x, 0));
	unsigned idx = std::distance(elements_.begin(), where);
	if		(where == elements_.end())									{	return std::numeric_limits<double>::signaling_NaN();	}
	else if	((where == elements_.begin()) && (x != elements_[0].x))		{	return std::numeric_limits<double>::signaling_NaN();	}
	else																{	return interpolate(x, idx-1, idx);	}
}

double LookupTable::operator()(double x) const{
	return evaluate(x);
}

void LookupTable::extend(const LookupTable & other){
	std::vector<Element> new_elements(elements_.size() + other.elements_.size());
	std::merge(elements_.begin(), elements_.end(), other.elements_.begin(), other.elements_.end(), new_elements.begin());
	elements_ = std::move(new_elements);
	auto where = std::unique(elements_.begin(), elements_.end(), [](const Element & a, const Element & b){return a.x == b.x;});
	elements_.erase(where, elements_.end());
}

double LookupTable::interpolate(double x, unsigned a, unsigned b) const{
	double p = (x - elements_[a].x) / (elements_[b].x - elements_[a].x);
	return p * elements_[b].f + (1-p) * elements_[a].f;
}
