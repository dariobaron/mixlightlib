#pragma once

#include <type_traits>


template<typename T>
concept Numeric = std::is_arithmetic_v<T>;