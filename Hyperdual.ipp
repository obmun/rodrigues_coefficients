#include <cstdio>

template<class Real>
Hyperdual<Real>::Hyperdual()
{
     f0 = 0.0;
     f1 = 0.0;
     f2 = 0.0;
     f12 = 0.0;
}

template<class Real>
Hyperdual<Real>::Hyperdual(Real x1, Real x2, Real x3, Real x4)
{
     f0 = x1;
     f1 = x2;
     f2 = x3;
     f12 = x4;
}

template<class Real>
Hyperdual<Real>::Hyperdual(Real x1)
{
     f0 = x1;
     f1 = 0.0;
     f2 = 0.0;
     f12 = 0.0;
}

template<class Real>
void Hyperdual<Real>::setvalues(Real x1, Real x2, Real x3, Real x4)
{
     f0 = x1;
     f1 = x2;
     f2 = x3;
     f12 = x4;
}

template<class Real>
void Hyperdual<Real>::view()
{
     printf("%g + %g epsilon1 + %g epsilon2 + %g epsilon1 epsilon2\n", f0, f1, f2, f12);
}

template<class Real>
Real Hyperdual<Real>::real()
{
     return f0;
}

template<class Real>
Real Hyperdual<Real>::eps1()
{
     return f1;
}

template<class Real>
Real Hyperdual<Real>::eps2()
{
     return f2;
}

template<class Real>
Real Hyperdual<Real>::eps1eps2()
{
     return f12;
}

template<class Real>
std::ostream& operator<<(std::ostream& output, const Hyperdual<Real>& rhs)
{
     output << "(" << rhs.f0 << ","<< rhs.f1 << ","<< rhs.f2 << ","<< rhs.f12 << ")";
     return output;
}

template<class Real>
Hyperdual<Real> Hyperdual<Real>::operator+() const
{
     return *this;
}

template<class Real>
Hyperdual<Real> Hyperdual<Real>::operator+(const Hyperdual<Real> &rhs) const
{
     Hyperdual<Real> temp;
     temp.f0 = f0 + rhs.f0;
     temp.f1 = f1 + rhs.f1;
     temp.f2 = f2 + rhs.f2;
     temp.f12 = f12 + rhs.f12;
     return temp;
}

template<class Real>
Hyperdual<Real> operator+(const Real lhs, const Hyperdual<Real> &rhs)
{
     Hyperdual<Real> temp;
     temp.f0 = lhs + rhs.f0;
     temp.f1 = rhs.f1;
     temp.f2 = rhs.f2;
     temp.f12 = rhs.f12;
     return temp;
}

template<class Real>
Hyperdual<Real> Hyperdual<Real>::operator-() const
{
     Hyperdual<Real> temp;
     temp.f0 = -f0;
     temp.f1 = -f1;
     temp.f2 = -f2;
     temp.f12 = -f12;
     return temp;
}

template<class Real>
Hyperdual<Real> Hyperdual<Real>::operator-(const Hyperdual<Real> &rhs) const
{
     Hyperdual<Real> temp;
     temp.f0 = f0 - rhs.f0;
     temp.f1 = f1 - rhs.f1;
     temp.f2 = f2 - rhs.f2;
     temp.f12 = f12 - rhs.f12;
     return temp;
}

template<class Real>
Hyperdual<Real> operator-(const Real lhs, const Hyperdual<Real> &rhs)
{
	Hyperdual<Real> temp;
	temp.f0 = lhs - rhs.f0;
	temp.f1 = -rhs.f1;
	temp.f2 = -rhs.f2;
	temp.f12 = -rhs.f12;
	return temp;
}

template<class Real>
Hyperdual<Real> Hyperdual<Real>::operator* (const Hyperdual<Real> &rhs) const
{
	Hyperdual temp;
	temp.f0 = f0*rhs.f0;
	temp.f1 = f0*rhs.f1 + f1*rhs.f0;
	temp.f2 = f0*rhs.f2 + f2*rhs.f0;
	temp.f12 = f0*rhs.f12 + f1*rhs.f2 + f2*rhs.f1 + f12*rhs.f0;
	return temp;
}

template<class Real>
Hyperdual<Real> operator*(const Real lhs, const Hyperdual<Real> &rhs)
{
	Hyperdual<Real> temp;
	temp.f0 = lhs*rhs.f0;
	temp.f1 = lhs*rhs.f1;
	temp.f2 = lhs*rhs.f2;
	temp.f12 = lhs*rhs.f12;
	return temp;
}

template<class Real>
Hyperdual<Real> operator/(const Hyperdual<Real> &lhs, const Hyperdual<Real> &rhs)
{
	Hyperdual<Real> temp, inv;
	inv = pow(rhs, Real(-1.0));
	temp = lhs*inv;
	return temp;
}

template<class Real>
Hyperdual<Real> operator/(const Real lhs, const Hyperdual<Real> &rhs)
{
	Hyperdual<Real> temp, inv;
	inv = pow(rhs,-1);
	temp = lhs*inv;
	return temp;
}

template<class Real>
Hyperdual<Real> operator/(const Hyperdual<Real> &lhs, const Real rhs)
{
	Hyperdual<Real> temp;
	Real inv;
	inv = 1.0/rhs;
	temp.f0 = inv*lhs.f0;
	temp.f1 = inv*lhs.f1;
	temp.f2 = inv*lhs.f2;
	temp.f12 = inv*lhs.f12;
	return temp;
}

template<class Real>
Hyperdual<Real>& Hyperdual<Real>::operator+=(const Hyperdual<Real> &rhs)
{
	f0 += rhs.f0;
	f1 += rhs.f1;
	f2 += rhs.f2;
	f12 += rhs.f12;
	return *this;
}

template<class Real>
Hyperdual<Real>& Hyperdual<Real>::operator-=(const Hyperdual<Real> &rhs)
{
	f0 -= rhs.f0;
	f1 -= rhs.f1;
	f2 -= rhs.f2;
	f12 -= rhs.f12;
	return *this;
}

template<class Real>
Hyperdual<Real>& Hyperdual<Real>::operator*=(const Hyperdual<Real> &rhs)
{
	Real tf0,tf1,tf2,tf12;
	tf0 = f0;
	tf1 = f1;
	tf2 = f2;
	tf12 = f12;
	f0 = tf0*rhs.f0;
	f1 = tf0*rhs.f1 + tf1*rhs.f0;
	f2 = tf0*rhs.f2 + tf2*rhs.f0;
	f12 = tf0*rhs.f12 + tf1*rhs.f2 + tf2*rhs.f1 + tf12*rhs.f0;
	return *this;
}

template<class Real>
Hyperdual<Real>& Hyperdual<Real>::operator*=(Real rhs)
{
	f0 *= rhs;
	f1 *= rhs;
	f2 *= rhs;
	f12 *= rhs;
	return *this;
}

template<class Real>
Hyperdual<Real>& Hyperdual<Real>::operator/=(Real rhs)
{
	f0 /= rhs;
	f1 /= rhs;
	f2 /= rhs;
	f12 /= rhs;
	return *this;
}

template<class Real>
Hyperdual<Real> pow(const Hyperdual<Real> &x, Real a)
{
	Hyperdual<Real> temp;
	Real deriv, xval, tol;
	xval = x.f0;
	tol = 1e-15;
	if (fabs(xval) < tol)
	{
		if (xval >= 0)
			xval = tol;
		if (xval < 0)
			xval = -tol;
	}
	deriv = a*pow(xval,(a-1));
	// temp.f0 = pow(xval,a);
	temp.f0 = pow(x.f0,a); // Use actual x value, only use tol for derivs
	temp.f1 = x.f1*deriv;
	temp.f2 = x.f2*deriv;
	temp.f12 = x.f12*deriv + a*(a-1)*x.f1*x.f2*pow(xval,(a-2));
	
	return temp;
}

template<class Real>
Hyperdual<Real> pow(const Hyperdual<Real> &x, const Hyperdual<Real> &a)
{
	return exp(a*log(x));
}

template<class Real>
Hyperdual<Real> exp(const Hyperdual<Real> &x)
{
	Hyperdual<Real> temp;
	Real deriv;
	deriv = exp(x.f0);
	temp.f0 = deriv;
	temp.f1 = deriv*x.f1;
	temp.f2 = deriv*x.f2;
	temp.f12 = deriv*(x.f12 + x.f1*x.f2);
	return temp;
}

template<class Real>
Hyperdual<Real> log(const Hyperdual<Real> &x)
{
	Hyperdual<Real> temp;
	Real deriv1, deriv2;
	deriv1 = x.f1/x.f0;
	deriv2 = x.f2/x.f0;
	temp.f0 = log(x.f0);
	temp.f1 = deriv1;
	temp.f2 = deriv2;
	temp.f12 = x.f12/x.f0 - (deriv1*deriv2);
	return temp;
}

template<class Real>
Hyperdual<Real> sin(const Hyperdual<Real> &x)
{
	Hyperdual<Real> temp;
	Real funval, deriv;
	funval = sin(x.f0);
	deriv = cos(x.f0);
	temp.f0 = funval;
	temp.f1 = deriv*x.f1;
	temp.f2 = deriv*x.f2;
	temp.f12 = deriv*x.f12 - funval*x.f1*x.f2;
	return temp;
}

template<class Real>
Hyperdual<Real> cos(const Hyperdual<Real> &x)
{
	Hyperdual<Real> temp;
	Real funval,deriv;
	funval = cos(x.f0);
	deriv = -sin(x.f0);
	temp.f0 = funval;
	temp.f1 = deriv*x.f1;
	temp.f2 = deriv*x.f2;
	temp.f12 = deriv*x.f12 - funval*x.f1*x.f2;
	return temp;
}

template<class Real>
Hyperdual<Real> tan(const Hyperdual<Real> &x)
{
	Hyperdual<Real> temp;
	Real funval,deriv;
	funval = tan(x.f0);
	deriv  = funval*funval + 1.0;
	temp.f0 = funval;
	temp.f1 = deriv*x.f1;
	temp.f2 = deriv*x.f2;
	temp.f12 = deriv*x.f12 + x.f1*x.f2*(2*funval*deriv);
	return temp;
}

template<class Real>
Hyperdual<Real> asin(const Hyperdual<Real> &x)
{
	Hyperdual<Real> temp;
	Real funval,deriv1,deriv;
	funval = asin(x.f0);
	deriv1 = 1.0-x.f0*x.f0;
	deriv = 1.0/sqrt(deriv1);
	temp.f0 = funval;
	temp.f1 = deriv*x.f1;
	temp.f2 = deriv*x.f2;
	temp.f12 = deriv*x.f12 + x.f1*x.f2*(x.f0*pow(deriv1,-1.5));
	return temp;
}

template<class Real>
Hyperdual<Real> acos(const Hyperdual<Real> &x)
{
	Hyperdual<Real> temp;
	Real funval,deriv1,deriv;
	funval = acos(x.f0);
	deriv1 = 1.0-x.f0*x.f0;
	deriv = -1.0/sqrt(deriv1);
	temp.f0 = funval;
	temp.f1 = deriv*x.f1;
	temp.f2 = deriv*x.f2;
	temp.f12 = deriv*x.f12 + x.f1*x.f2*(-x.f0*pow(deriv1,-1.5));
	return temp;
}

template<class Real>
Hyperdual<Real> atan(const Hyperdual<Real> &x)
{
	Hyperdual<Real> temp;
	Real funval,deriv1,deriv;
	funval = atan(x.f0);
	deriv1 = 1.0+x.f0*x.f0;
	deriv = 1.0/deriv1;
	temp.f0 = funval;
	temp.f1 = deriv*x.f1;
	temp.f2 = deriv*x.f2;
	temp.f12 = deriv*x.f12 + x.f1*x.f2*(-2*x.f0/(deriv1*deriv1));
	return temp;
}

template<class Real>
Hyperdual<Real> sqrt(const Hyperdual<Real> &x)
{
	return pow(x,0.5);
}

template<class Real>
Hyperdual<Real> fabs(const Hyperdual<Real> &x)
{
	Hyperdual<Real> temp;
	if (x < 0.0)
		temp = -x;
	else
		temp = x;
	return temp;
}

template<class Real>
Hyperdual<Real> max(const Hyperdual<Real> &x1, const Hyperdual<Real> &x2)
{
	if (x1 > x2)
		return x1;
	return x2;
}

template<class Real>
Hyperdual<Real> max(const Hyperdual<Real> &x1, Real x2)
{
	if (x1 > x2)
		return x1;
	return x2;
}

template<class Real>
Hyperdual<Real> max(Real x1, const Hyperdual<Real> &x2)
{
	if (x1 > x2)
		return x1;
	return x2;
}

template<class Real>
Hyperdual<Real> min(const Hyperdual<Real> &x1, const Hyperdual<Real> &x2)
{
	if (x1 < x2)
          return x1;
     return x2;
}

template<class Real>
Hyperdual<Real> min(const Hyperdual<Real> &x1, Real x2)
{
	if (x1 < x2)
		return x1;
     return x2;
}

template<class Real>
Hyperdual<Real> min(Real x1, const Hyperdual<Real> &x2)
{
	if (x1 < x2)
		return x1;
	return x2;
}

template<class Real>
bool operator>(const Hyperdual<Real> &lhs, const Hyperdual<Real> &rhs)
{
	return (lhs.f0 > rhs.f0);
}

template<class Real>
bool operator>(Real lhs, const Hyperdual<Real> &rhs)
{
	return (lhs > rhs.f0);
}

template<class Real>
bool operator>(const Hyperdual<Real> &lhs, Real rhs)
{
	return (lhs.f0 > rhs);
}

template<class Real>
bool operator>=(const Hyperdual<Real> &lhs, const Hyperdual<Real> &rhs)
{
	return (lhs.f0 >= rhs.f0);
}

template<class Real>
bool operator>=(Real lhs, const Hyperdual<Real> &rhs)
{
	return (lhs >= rhs.f0);
}

template<class Real>
bool operator>=(const Hyperdual<Real> &lhs, Real rhs)
{
	return (lhs.f0 >= rhs);
}

template<class Real>
bool operator<(const Hyperdual<Real> &lhs, const Hyperdual<Real> &rhs)
{
	return (lhs.f0 < rhs.f0);
}

template<class Real>
bool operator<(Real lhs, const Hyperdual<Real> &rhs)
{
	return (lhs < rhs.f0);
}

template<class Real>
bool operator<(const Hyperdual<Real> &lhs, Real rhs)
{
	return (lhs.f0 < rhs);
}

template<class Real>
bool operator<=(const Hyperdual<Real> &lhs, const Hyperdual<Real> &rhs)
{
	return (lhs.f0 <= rhs.f0);
}

template<class Real>
bool operator<=(Real lhs, const Hyperdual<Real> &rhs)
{
	return (lhs <= rhs.f0);
}

template<class Real>
bool operator<=(const Hyperdual<Real> &lhs, Real rhs)
{
	return (lhs.f0 <= rhs);
}

template<class Real>
bool operator==(const Hyperdual<Real> &lhs, const Hyperdual<Real> &rhs)
{
	return (lhs.f0 == rhs.f0);
}

template<class Real>
bool operator==(Real lhs, Hyperdual<Real> rhs)
{
	return (lhs == rhs.f0);
}

template<class Real>
bool operator== (Hyperdual<Real> lhs, Real rhs)
{
	return (lhs.f0 == rhs);
}

template<class Real>
bool operator!= (Hyperdual<Real> lhs, Hyperdual<Real> rhs)
{
	return (lhs.f0 != rhs.f0);
}

template<class Real>
bool operator!= (Real lhs, Hyperdual<Real> rhs)
{
	return (lhs != rhs.f0);
}

template<class Real>
bool operator!= (Hyperdual<Real> lhs, Real rhs)
{
	return (lhs.f0 != rhs);
}
