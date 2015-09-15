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
#include <cassert>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>
#include "hyperdual.hpp"

constexpr unsigned long int factorial(unsigned long int n)
{
    return n <= 1 ? 1 : (n * factorial(n - 1));
}

/**
 *
 * b_i = \frac{1}{\theta} \diff{a_i(\theta)}{\theta}
 * c_i = \frac{1}{\theta} \diff{b_i(\theta)}{\theta}
 */
namespace rodrigues_formula
{
     
     enum class CalculationMode { Direct, NumericHyperDual, SeriesExpansion };
    
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
	       a2(m_impl),
	       b0(m_impl),
	       b1(m_impl),
	       b2(m_impl) {
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

	  class B0
	  {
	  public:
	       B0(const TrigonometricCoeffs::Impl &impl) : m_impl(impl) { }
	       T operator()(T theta) const {
		    return m_impl.b0(theta);
	       }

	  protected:
	       const TrigonometricCoeffs::Impl &m_impl;
	  };

	  class B1
	  {
	  public:
	       B1(const TrigonometricCoeffs::Impl &impl) : m_impl(impl) { }
	       T operator()(T theta) const {
		    return m_impl.b1(theta);
	       }

	  protected:
	       const TrigonometricCoeffs::Impl &m_impl;
	  };

	  class B2
	  {
	  public:
	       B2(const TrigonometricCoeffs::Impl &impl) : m_impl(impl) { }
	       T operator()(T theta) const {
		    return m_impl.b2(theta);
	       }

	  protected:
	       const TrigonometricCoeffs::Impl &m_impl;
	  };


	  const A0 a0;
	  const A1 a1;
	  const A2 a2;
	  const B0 b0;
	  const B1 b1;
	  const B2 b2;

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
		    return (T(1) - cos(theta)) / (theta * theta);
	       }

	       static T da0(T theta) {
		    return -sin(theta);
	       }

	       static T da1(T theta) {
		    return (theta * cos(theta) - sin(theta)) / (theta * theta);
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

	       static T b0(T theta) {
		    return -sin(theta) / theta;
	       }

	       /**
		* b_1 = \frac{1}{\theta} \diff{a_1(\theta)}{\theta}
		*/
	       static T b1(T theta) {
		    return (theta * cos(theta) - sin(theta)) / pow(theta, 3);
	       }

	       /**
		* b_2 = \frac{1}{\theta} \diff{a_2(\theta)}{\theta}
		*/
	       static T b2(T theta) {
		    return (theta * sin(theta) + T(2) * cos(theta) - T(2)) / pow(theta, 4);
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

	       RealType b0(RealType theta) const {
		    return da0(theta) / theta;
	       }

	       RealType b1(RealType theta) const {
		    return da1(theta) / theta;
	       }

	       RealType b2(RealType theta) const {
		    return da2(theta) / theta;
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

	  template <typename T>
	  class TrigonometricCoeffsImpl<T, CalculationMode::SeriesExpansion>
	  {
	  public:
	       static T a0(T theta) {
		    return s_direct.a0(theta);
	       }
	       
	       static T a1(T theta) {
		    if (theta > S_THRESHOLD) return s_direct.a1(theta);
		    return ai(1, theta);
		    
	       }

	       static T a2(T theta) {
		    if (theta > S_THRESHOLD) return s_direct.a2(theta);
		    return ai(2, theta);
	       }

	       static T b0(T theta) {
		    if (theta > S_THRESHOLD) return s_direct.b0(theta);
		    return bi(0, theta);
	       }

	       static T b1(T theta) {
		    if (theta > S_THRESHOLD) return s_direct.b1(theta);
		    return bi(1, theta);
	       }

	       static T b2(T theta) {
		    if (theta > S_THRESHOLD) return s_direct.b2(theta);
		    return bi(2, theta);
	       }

	  protected:
	       static constexpr T S_ONE = 1.0;
	       static constexpr T S_THRESHOLD = 0.25;
	       static const int N_FACTORIALS = 15;
	       static const std::array<T,N_FACTORIALS> S_INV_FACTORIALS;
	       static class TrigonometricCoeffsImpl<T, CalculationMode::Direct> s_direct;

#define theta_powers(theta)					\
	       T theta2, theta4, theta6, theta8, theta10;	\
	       theta2 = (theta) * (theta);			\
	       theta4 = theta2 * theta2;			\
	       theta6 = theta2 * theta4;			\
	       theta8 = theta4 * theta4;			\
	       theta10 = theta8 * theta2;			\

	       static T ai(unsigned int i, T theta) {
		    assert(i < 4);
		    constexpr int N_STEPS = 6;
		    theta_powers(theta);
		    T s[N_STEPS] = { 1, -theta2, theta4, -theta6, theta8, -theta10 };
		    for (unsigned int j = 0; j < N_STEPS; j++) {
			 s[j] *= S_INV_FACTORIALS[2*j + i];
		    }
		    T res = 0.;
		    for (unsigned int j = 0; j < N_STEPS; j++)
		    {
			 res += s[j];
		    }
		    return res;
	       }
	       
	       static T bi(unsigned int i, T theta) {
		    assert(i < 3);
		    constexpr int N_STEPS = 6;
		    theta_powers(theta);
		    T s[N_STEPS] = { -2, 4*theta2, -6*theta4, 8*theta6, -10*theta8, 12*theta10 };
		    for (unsigned int j = 0; j < N_STEPS; j++)
		    {
			 s[j] *= S_INV_FACTORIALS[2 + 2*j + i];
			 // Max factorial idx: 2 + 2*5 + 3 = 15 -> fits
		    }
		    T res = 0.;
		    for (unsigned int j = 0; j < N_STEPS; j++)
		    {
			 res += s[j];
		    }
		    return res;
	       }
	  };

	  template <typename T> std::array<T,TrigonometricCoeffsImpl<T, CalculationMode::SeriesExpansion>::N_FACTORIALS> const
	  TrigonometricCoeffsImpl<T, CalculationMode::SeriesExpansion>::S_INV_FACTORIALS =
	  { S_ONE / factorial(0), S_ONE / factorial(1), S_ONE / factorial(2), S_ONE / factorial(3), S_ONE / factorial(4),
	    S_ONE / factorial(5), S_ONE / factorial(6), S_ONE / factorial(7), S_ONE / factorial(8), S_ONE / factorial(9),
	    S_ONE / factorial(10), S_ONE / factorial(11), S_ONE / factorial(12), S_ONE / factorial(13), S_ONE / factorial(14) };

     }

}

namespace rf = rodrigues_formula;

int main(int argc, char *argv[])
{
     using namespace std::placeholders;

     typedef float RealType;

     typedef rf::TrigonometricCoeffs<RealType, rf::CalculationMode::Direct> TCsDir;
     typedef rf::TrigonometricCoeffs<double, rf::CalculationMode::NumericHyperDual> TCsHD;
     typedef rf::TrigonometricCoeffs<RealType, rf::CalculationMode::SeriesExpansion> TCsSE;

     const RealType STEP = 1e-2;
     const int N_EVAL_PTS = 101;
     std::vector<RealType> eval_pts;
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
     TCsSE tcs_se;
     auto tcs_hd_impl = tcs_hd.impl();
     tcs_hd_impl.set_steps(1e-14, 1e-14);

     // a_0(0.0) -> tcs_dir.a0(0.0);

     std::map<std::string, std::function<RealType(RealType)>> derivs;
     derivs["a0"] = tcs_dir.a0;
     derivs["a1"] = tcs_dir.a1;
     derivs["a2"] = tcs_dir.a2;
     derivs["b0"] = tcs_dir.b0;
     derivs["b1"] = tcs_dir.b1;
     derivs["b2"] = tcs_dir.b2;

     std::map<std::string, decltype(derivs)> all_derivs;
     all_derivs["direct"] = std::move(derivs);

#if 0
     derivs.clear();
     // derivs["d(a0)/dtheta"] = [&](double v) { return tcs_hd.d(tcs_hd.a0, v); };
     // derivs["d(a1)/dtheta"] = [&](double v) { return tcs_hd.d(tcs_hd.a1, v); };
     // derivs["d(a2)/dtheta"] = [&](double v) { return tcs_hd.d(tcs_hd.a2, v); };
     // derivs["d^2(a0)/dtheta^2"] = [&](double v) { return tcs_hd.d2(tcs_hd.a0, v); };
     // derivs["d^2(a1)/dtheta^2"] = [&](double v) { return tcs_hd.d2(tcs_hd.a1, v); };
     // derivs["d^2(a2)/dtheta^2"] = [&](double v) { return tcs_hd.d2(tcs_hd.a2, v); };
     derivs["a0"] = tcs_hd.a0;
     derivs["a1"] = tcs_hd.a1;
     derivs["a2"] = tcs_hd.a2;
     derivs["b0"] = tcs_hd.b0;
     derivs["b1"] = tcs_hd.b1;
     derivs["b2"] = tcs_hd.b2;
     all_derivs["hyperdual"] = std::move(derivs);
#endif

     derivs.clear();
     derivs["a0"] = tcs_se.a0;
     derivs["a1"] = tcs_se.a1;
     derivs["a2"] = tcs_se.a2;
     derivs["b0"] = tcs_se.b0;
     derivs["b1"] = tcs_se.b1;
     derivs["b2"] = tcs_se.b2;

     size_t max_name_len = 0;
     for (auto &deriv : derivs)
     {
	  max_name_len = std::max(max_name_len, deriv.first.length());
     }

     all_derivs["series"] = std::move(derivs);
     derivs.clear();

     std::map<std::string, std::map<std::string, std::vector<RealType>>> results;

     for (auto &derivs : all_derivs)
     {
	  for (auto &f : derivs.second)
	  {
	       std::vector<RealType> res;
	       std::transform(std::begin(eval_pts), std::end(eval_pts), std::back_inserter(res), f.second);
	       results[derivs.first][f.first] = std::move(res);
	  }
     }

     const int WIDTH = 14;
     std::cout << std::scientific;
     std::cout << std::setprecision(7);
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

#define print_line()							\
     for (unsigned int i = 0; i < max_name_len + (WIDTH + SEPARATOR.length()) * eval_pts.size(); ++i) \
     {									\
	  std::cout << "-";						\
     }									\
     std::cout << "\n";

     print_line();
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
	  print_line();
     }

     return 0;
}
