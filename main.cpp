/**
 * @file
 * @brief Evaluation of implementation of Rodrigues formula's trigonometric a_i coefficients using
 * hyper-dual numeric derivation.
 *
 * This file implements the a_i expressions (first equations of section 2.2 of the paper "On the
 * differentiation of the Rodrigues formula and its significance for the vector-like
 * parameterization of Reissner–Simo beam theory") using Fike's hyper-dual numerical derivation
 * method.
 *
 * For comparison, implements the expressions as well directly using the symbolic expressions for
 * the first order and second order derivatives.
 *
 * This code is written in C++11 and only makes uses of std library functions and the included
 * hyper-dual class as implemented by Fike. Build system is CMake.
 */

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>
#include "hyperdual.hpp"

namespace rodrigues_formula
{
     
     enum class CalculationMode { Direct, NumericHyperDual };
    
     namespace detail
     {
	  template < typename T, CalculationMode mode >
	  class DependentFalse : std::false_type
	  { };

	  template < typename T, CalculationMode mode > class TrigonometricCoeffsImpl
	  {
	  public:
	       static T a0(T theta) {
		    static_assert(DependentFalse<T, mode>::value, "no default implementation");
	       }
	       
	       static T a1(T theta) {
		    static_assert(DependentFalse<T, mode>::value, "no default implementation");
	       }
	       
	       static T a2(T theta) {
		    static_assert(DependentFalse<T, mode>::value, "no default implementation");
	       }
	  };
     }

     template < typename T, CalculationMode mode > class TrigonometricCoeffs
     {
     protected:
	  typedef detail::TrigonometricCoeffsImpl<T, mode> Impl;
	  Impl m_impl;

     public:
	  TrigonometricCoeffs() :
	       a0(m_impl),
	       a1(m_impl),
	       a2(m_impl) {
	  }

	  Impl &impl() {
	       return m_impl;
	  }

	  class A0
	  {
	  public:
	       A0(const TrigonometricCoeffs::Impl &impl) : m_impl(impl) { }
	       T operator()(T theta) const {
		    return m_impl.a0(theta);
	       }

	  protected:
	       const TrigonometricCoeffs::Impl &m_impl;
	  };

	  class A1
	  {
	  public:
	       A1(const TrigonometricCoeffs::Impl &impl) : m_impl(impl) { }
	       T operator()(T theta) const {
		    return m_impl.a1(theta);
	       }

	  protected:
	       const TrigonometricCoeffs::Impl &m_impl;
	  };

	  class A2
	  {
	  public:
	       A2(const TrigonometricCoeffs::Impl &impl) : m_impl(impl) { }
	       T operator()(T theta) const {
		    return m_impl.a2(theta);
	       }

	  protected:
	       const TrigonometricCoeffs::Impl &m_impl;
	  };

	  const A0 a0;
	  const A1 a1;
	  const A2 a2;

	  T d(const A0 &, T theta) const {
	       return m_impl.da0(theta);
	  }

	  T d(const A1 &, T theta) const {
	       return m_impl.da1(theta);
	  }

	  T d(const A2 &, T theta) const {
	       return m_impl.da2(theta);
	  }

	  T d2(const A0 &, T theta) const {
	       return m_impl.d2a0(theta);
	  }

	  T d2(const A1 &, T theta) const {
	       return m_impl.d2a1(theta);
	  }

	  T d2(const A2 &, T theta) const {
	       return m_impl.d2a2(theta);
	  }
     };

     namespace detail
     {
	  template <typename T>
	  class TrigonometricCoeffsImpl<T, CalculationMode::Direct>
	  {
	  public:
	       static T a0(T theta) {
		    return cos(theta);
	       }
	       
	       static T a1(T theta) {
		    return sin(theta) / theta;
	       }

	       static T a2(T theta) {
		    return (T(1) - cos(theta)) / pow(theta, 2);
	       }

	       static T da0(T theta) {
		    return -sin(theta);
	       }

	       static T da1(T theta) {
		    return (theta * cos(theta) - sin(theta)) / pow(theta, 2);
	       }

	       static T da2(T theta) {
		    return (theta * sin(theta) + T(2) * cos(theta) - T(2)) / pow(theta, 3);
	       }

	       static T d2a0(T theta) {
		    return -cos(theta);
	       }

	       static T d2a1(T theta) {
		    return -((pow(theta, 2) - 2) * sin(theta) + 2 * theta * cos(theta)) / pow(theta, 3);
	       }

	       static T d2a2(T theta) {
		    return ((pow(theta, 2) - 6) * cos(theta) - 4 * theta * sin(theta) + 6) / pow(theta, 4);
	       }
	  };

	  template <>
	  class TrigonometricCoeffsImpl< double, CalculationMode::NumericHyperDual >
	  {
	  public:
	       typedef double RealType;

	       TrigonometricCoeffsImpl() :
		    m_h1(1e-10),
		    m_h2(1e-10) {
	       }
	       
	       void set_steps(RealType h1, RealType h2) {
		    m_h1 = h1;
		    m_h2 = h2;
	       }

	       static RealType a0(RealType theta) {
		    return cos(theta);
	       }
	       
	       static RealType a1(RealType theta) {
		    return sin(theta) / theta;
	       }

	       static RealType a2(RealType theta) {
		    return (RealType(1) - cos(theta)) / pow(theta, 2);
	       }

	       RealType da0(RealType theta) const {
		    return _a0(theta).eps1() / m_h1;
	       }

	       RealType da1(RealType theta) const {
		    return _a1(theta).eps1() / m_h1;
	       }
	       
	       RealType da2(RealType theta) const {
		    return _a2(theta).eps1() / m_h1;
	       }

	       RealType d2a0(RealType theta) const {
		    return _a0(theta).eps1eps2() / (m_h1 * m_h2);
	       }

	       RealType d2a1(RealType theta) const {
		    return _a1(theta).eps1eps2() / (m_h1 * m_h2);
	       }

	       RealType d2a2(RealType theta) const {
		    return _a2(theta).eps1eps2() / (m_h1 * m_h2);
	       }

	  protected:
	       RealType m_h1, m_h2;

	       hyperdual _a0(RealType theta) const {
		    hyperdual theta_hat(theta, m_h1, m_h2, 0);
		    auto res = cos(theta_hat);
		    return res;
	       }

	       hyperdual _a1(RealType theta) const {
		    hyperdual theta_hat(theta, m_h1, m_h2, 0);
		    auto v = sin(theta_hat);
		    return v / theta_hat;
	       }
	       
	       hyperdual _a2(RealType theta) const {
		    hyperdual theta_hat(theta, m_h1, m_h2, 0);
		    auto v = hyperdual(1, 0, 0, 0) - cos(theta_hat);
		    return v / pow(theta_hat, 2);
	       }

	  };

     }

}

namespace rf = rodrigues_formula;

int main(int argc, char *argv[])
{
     using namespace std::placeholders;

     typedef rf::TrigonometricCoeffs< double, rf::CalculationMode::Direct > TCsDir;
     typedef rf::TrigonometricCoeffs< double, rf::CalculationMode::NumericHyperDual > TCsHD;

     const double STEP = 1e-7;
     const int N_EVAL_PTS = 21;
     std::vector<double> eval_pts;
     int m;
     unsigned int i;
     for (m = - N_EVAL_PTS / 2, i = 0;
	  i < N_EVAL_PTS;
	  m++, i++)
     {
	  eval_pts.push_back(m * STEP);
     }

     TCsDir tcs_dir;
     TCsHD tcs_hd;
     auto tcs_hd_impl = tcs_hd.impl();
     tcs_hd_impl.set_steps(1e-14, 1e-14);

     // a_0(0.0) -> tcs_dir.a0(0.0);

     std::map<std::string, std::function<double(double)>> derivs;
     derivs["d(a0)/dtheta"] = [&](double v) { return tcs_dir.d(tcs_dir.a0, v); };
     derivs["d(a1)/dtheta"] = [&](double v) { return tcs_dir.d(tcs_dir.a1, v); };
     derivs["d(a2)/dtheta"] = [&](double v) { return tcs_dir.d(tcs_dir.a2, v); };
     derivs["d^2(a0)/dtheta^2"] = [&](double v) { return tcs_dir.d2(tcs_dir.a0, v); };
     derivs["d^2(a1)/dtheta^2"] = [&](double v) { return tcs_dir.d2(tcs_dir.a1, v); };
     derivs["d^2(a2)/dtheta^2"] = [&](double v) { return tcs_dir.d2(tcs_dir.a2, v); };

     std::map<std::string, decltype(derivs)> all_derivs;
     all_derivs["direct"] = std::move(derivs);

     derivs.clear();
     derivs["d(a0)/dtheta"] = [&](double v) { return tcs_hd.d(tcs_hd.a0, v); };
     derivs["d(a1)/dtheta"] = [&](double v) { return tcs_hd.d(tcs_hd.a1, v); };
     derivs["d(a2)/dtheta"] = [&](double v) { return tcs_hd.d(tcs_hd.a2, v); };
     derivs["d^2(a0)/dtheta^2"] = [&](double v) { return tcs_hd.d2(tcs_hd.a0, v); };
     derivs["d^2(a1)/dtheta^2"] = [&](double v) { return tcs_hd.d2(tcs_hd.a1, v); };
     derivs["d^2(a2)/dtheta^2"] = [&](double v) { return tcs_hd.d2(tcs_hd.a2, v); };

     size_t max_name_len = 0;
     for (auto &deriv : derivs)
     {
	  max_name_len = std::max(max_name_len, deriv.first.length());
     }

     all_derivs["hyperdual"] = std::move(derivs);
     derivs.clear();

     std::map<std::string, std::map<std::string, std::vector<double>>> results;

     for (auto &derivs : all_derivs)
     {
	  for (auto &f : derivs.second)
	  {
	       std::vector<double> res;
	       std::transform(std::begin(eval_pts), std::end(eval_pts), std::back_inserter(res), f.second);
	       results[derivs.first][f.first] = std::move(res);
	  }
     }

     const int WIDTH = 10;
     std::cout << std::scientific;
     std::cout << std::setprecision(3);
     const std::string SEPARATOR(" | ");
     for (unsigned int i = 0; i < max_name_len; ++i)
     {
	  std::cout << " ";
     }
     for (auto v : eval_pts)
     {
	  std::cout << SEPARATOR << std::setw(WIDTH) << v;
     }
     std::cout << "\n";
     for (unsigned int i = 0; i < max_name_len + (WIDTH + SEPARATOR.length()) * eval_pts.size(); ++i)
     {
	  std::cout << "-";
     }
     std::cout << "\n";
     for (auto &group : results)
     {
	  for (auto &results : group.second)
	  {
	       std::cout << std::setw(max_name_len) << results.first;
	       for (auto &v : results.second)
	       {
		    std::cout << SEPARATOR << std::setw(WIDTH) << v;
	       }
	       std::cout << "\n";
	  }
     }

     return 0;
}
