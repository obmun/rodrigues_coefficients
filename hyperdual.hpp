/*
 * Class: hyperdual
 * 
 * Implementation of hyper-dual numbers
 *
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

class hyperdual
{
     double f0, f1, f2, f12;

public:
     //creation operators and function to manually set values
     hyperdual();
     hyperdual(double x1, double x2, double x3, double x4);
     hyperdual(double x1);
     void setvalues(double x1, double x2, double x3, double x4);

     //examine values
     void view(void);
     double real(void);
     double eps1(void);
     double eps2(void);
     double eps1eps2(void);
     friend std::ostream& operator<<(std::ostream& output, const hyperdual& rhs);

     //basic manipulation
     hyperdual operator+ () const;
     hyperdual operator+ (const hyperdual rhs) const;
     friend hyperdual operator+ (const double lhs, const hyperdual rhs);
     hyperdual operator- () const;
     hyperdual operator- (const hyperdual rhs) const;
     friend hyperdual operator- (const double lhs, const hyperdual rhs);
     hyperdual operator* (const hyperdual rhs)const;
     friend hyperdual operator* (const double lhs, const hyperdual rhs);
     friend hyperdual operator/ (const hyperdual lhs, const hyperdual rhs);
     friend hyperdual operator/ (const double lhs, const hyperdual rhs);
     friend hyperdual operator/ (const hyperdual lhs, const double rhs);
     hyperdual& operator+= (hyperdual rhs);
     hyperdual& operator-= (hyperdual rhs);
     hyperdual& operator*= (hyperdual rhs);
     hyperdual& operator*= (double rhs);
     hyperdual& operator/= (double rhs);

     // math.h functions
     friend hyperdual pow (hyperdual x, double a);
     friend hyperdual pow (hyperdual x, hyperdual a);
     friend hyperdual exp(hyperdual x);
     friend hyperdual log(hyperdual x);
     friend hyperdual sin(hyperdual x);
     friend hyperdual cos(hyperdual x);
     friend hyperdual tan(hyperdual x);
     friend hyperdual asin(hyperdual x);
     friend hyperdual acos(hyperdual x);
     friend hyperdual atan(hyperdual x);
     friend hyperdual sqrt(hyperdual x);
     friend hyperdual fabs(hyperdual x);
     friend hyperdual max(hyperdual x1, hyperdual x2);
     friend hyperdual max(hyperdual x1, double x2);
     friend hyperdual max(double x1, hyperdual x2);
     friend hyperdual min(hyperdual x1, hyperdual x2);
     friend hyperdual min(hyperdual x1, double x2);
     friend hyperdual min(double x1, hyperdual x2);

     // comparisons
     friend bool operator> (hyperdual lhs, hyperdual rhs);
     friend bool operator> (double lhs, hyperdual rhs);
     friend bool operator> (hyperdual lhs, double rhs);
     friend bool operator>= (hyperdual lhs, hyperdual rhs);
     friend bool operator>= (double lhs, hyperdual rhs);
     friend bool operator>= (hyperdual lhs, double rhs);
     friend bool operator< (hyperdual lhs, hyperdual rhs);
     friend bool operator< (double lhs, hyperdual rhs);
     friend bool operator< (hyperdual lhs, double rhs);
     friend bool operator<= (hyperdual lhs, hyperdual rhs);
     friend bool operator<= (double lhs, hyperdual rhs);
     friend bool operator<= (hyperdual lhs, double rhs);
     friend bool operator== (hyperdual lhs, hyperdual rhs);
     friend bool operator== (double lhs, hyperdual rhs);
     friend bool operator== (hyperdual lhs, double rhs);
     friend bool operator!= (hyperdual lhs, hyperdual rhs);
     friend bool operator!= (double lhs, hyperdual rhs);
     friend bool operator!= (hyperdual lhs, double rhs);
};

#endif

