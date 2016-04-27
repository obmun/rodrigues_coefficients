/*
 * Written by: Jeffrey A. Fike
 * Stanford University, Department of Aeronautics and Astronautics
 * 
 * Copyright (c) 2006 Jeffrey A. Fike
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */

#ifndef _hyperdual_h
#define _hyperdual_h

#include <iostream>
#include <math.h>

/**
 * @brief Implementation of hyper-dual numbers
 */
template<typename Real>
class Hyperdual
{
     Real f0, f1, f2, f12;

public:
     Hyperdual();
     Hyperdual(Real x1, Real x2, Real x3, Real x4);
     Hyperdual(Real x1);
     void setvalues(Real x1, Real x2, Real x3, Real x4);

     //examine values
     void view();
     Real real();
     Real eps1();
     Real eps2();
     Real eps1eps2();
     template<class R> friend std::ostream& operator<<(std::ostream& output, const Hyperdual<R>& rhs);

     // Basic manipulation
     Hyperdual<Real> operator+() const;
     Hyperdual<Real> operator+(const Hyperdual<Real> &rhs) const;
     template<class R> friend Hyperdual operator+(const R lhs, const Hyperdual<R> &rhs);
     Hyperdual<Real> operator-() const;
     Hyperdual<Real> operator-(const Hyperdual<Real> &rhs) const;
     template<class R> friend Hyperdual<R> operator-(const R lhs, const Hyperdual<R> &rhs);
     Hyperdual<Real> operator*(const Hyperdual<Real> &rhs)const;
     template<class R> friend Hyperdual<R> operator*(const R lhs, const Hyperdual<R> &rhs);
     template<class R> friend Hyperdual<R> operator/(const Hyperdual<R> &lhs, const Hyperdual<R> &rhs);
     template<class R> friend Hyperdual<R> operator/(const R lhs, const Hyperdual<R> &rhs);
     template<class R> friend Hyperdual<R> operator/(const Hyperdual<R> &lhs, const R rhs);
     Hyperdual<Real>& operator+=(const Hyperdual<Real> &rhs);
     Hyperdual<Real>& operator-=(const Hyperdual<Real> &rhs);
     Hyperdual<Real>& operator*=(const Hyperdual<Real> &rhs);
     Hyperdual<Real>& operator*=(Real rhs);
     Hyperdual<Real>& operator/=(Real rhs);

     // math.h functions
     template<class R> friend Hyperdual<R> pow(const Hyperdual<R> &x, R a);
     template<class R> friend Hyperdual<R> pow(const Hyperdual<R> &x, const Hyperdual<R> &a);
     template<class R> friend Hyperdual<R> exp(const Hyperdual<R> &x);
     template<class R> friend Hyperdual<R> log(const Hyperdual<R> &x);
     template<class R> friend Hyperdual<R> sin(const Hyperdual<R> &x);
     template<class R> friend Hyperdual<R> cos(const Hyperdual<R> &x);
     template<class R> friend Hyperdual<R> tan(const Hyperdual<R> &x);
     template<class R> friend Hyperdual<R> asin(const Hyperdual<R> &x);
     template<class R> friend Hyperdual<R> acos(const Hyperdual<R> &x);
     template<class R> friend Hyperdual<R> atan(const Hyperdual<R> &x);
     template<class R> friend Hyperdual<R> sqrt(const Hyperdual<R> &x);
     template<class R> friend Hyperdual<R> fabs(const Hyperdual<R> &x);
     template<class R> friend Hyperdual<R> max(const Hyperdual<R> &x1, const Hyperdual<R> &x2);
     template<class R> friend Hyperdual<R> max(const Hyperdual<R> &x1, R x2);
     template<class R> friend Hyperdual<R> max(R x1, const Hyperdual<R> &x2);
     template<class R> friend Hyperdual<R> min(const Hyperdual<R> &x1, const Hyperdual<R> &x2);
     template<class R> friend Hyperdual<R> min(const Hyperdual<R> &x1, R x2);
     template<class R> friend Hyperdual<R> min(R x1, const Hyperdual<R> &x2);

     // comparisons
     template<class R> friend bool operator>(const Hyperdual<R> lhs, const Hyperdual<R> &rhs);
     template<class R> friend bool operator>(R lhs, const Hyperdual<R> &rhs);
     template<class R> friend bool operator>(const Hyperdual<R> &lhs, R rhs);
     template<class R> friend bool operator>=(const Hyperdual<R> &lhs, const Hyperdual<R> &rhs);
     template<class R> friend bool operator>=(R lhs, const Hyperdual<R> &rhs);
     template<class R> friend bool operator>=(Hyperdual<R> lhs, R rhs);
     template<class R> friend bool operator<(Hyperdual<R> lhs, Hyperdual<R> rhs);
     template<class R> friend bool operator<(R lhs, Hyperdual<R> rhs);
     template<class R> friend bool operator<(Hyperdual<R> lhs, R rhs);
     template<class R> friend bool operator<=(Hyperdual<R> lhs, Hyperdual<R> rhs);
     template<class R> friend bool operator<=(R lhs, Hyperdual<R> rhs);
     template<class R> friend bool operator<=(Hyperdual<R> lhs, R rhs);
     template<class R> friend bool operator==(Hyperdual<R> lhs, Hyperdual<R> rhs);
     template<class R> friend bool operator==(R lhs, Hyperdual<R> rhs);
     template<class R> friend bool operator==(Hyperdual<R> lhs, R rhs);
     template<class R> friend bool operator!=(Hyperdual<R> lhs, Hyperdual<R> rhs);
     template<class R> friend bool operator!=(R lhs, Hyperdual<R> rhs);
     template<class R> friend bool operator!=(Hyperdual<R> lhs, R rhs);
};

#include "Hyperdual.ipp"

#endif

