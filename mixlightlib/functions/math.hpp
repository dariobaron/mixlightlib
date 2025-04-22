#pragma once

#include <algorithm>
#include <stdexcept>
#include <string>


unsigned long binomCoeff(unsigned n, unsigned k){
	k = std::min(k, n-k);
	if (k > n){
		std::string error = "Binomial coefficient with k > n (";
		error += std::to_string(k) + " > " + std::to_string(n) + ")";
		throw std::runtime_error(error);
	}
	unsigned long result = n - k + 1;
	for (unsigned i = 1; i < k; ++i){
		result *= (n - k + 1 + i) / (1 + i);
	}
	return result;
}