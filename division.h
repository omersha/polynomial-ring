// division.h

///////////////////////////////////////////////////////////////////////////////////////////////
// (1) "divides" - Boolean methods that check divisability a one entity by another,
// (2) safelyDivide - Perform a simple division of one term by another (assuming no reminder).
// (3) divide - A fully division algorithm.
//              (Polynomial is being divided by a sequence of polynomials).
///////////////////////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef division_H__
#define division_H__

#include <tuple>

#include "monomials.h"
#include "polynomials.h"


// Divisability (binary relations)
////////////////////////////////////////////////////////////////////////////

// Relation: monomial | monomial
template<typename PolyRing>
bool divides(Monomial<PolyRing> const &divisor, Monomial<PolyRing> const &dividend)
{
   for (size_t i = 0; i < PolyRing::VARIABLES; ++i)
      if (divisor[i] > dividend[i]) return false;
   return true;
}

// Relation: term | term
template<typename PolyRing>
bool divides(Term<PolyRing> const &divisor, Term<PolyRing> const &dividend)
{
   if (PolyRing::isZero(dividend.getCoeff())) return true;
   if (PolyRing::isZero(divisor.getCoeff())) return false;
   return divides(divisor.getMonomial(), dividend.getMonomial());
}
 
// Relation: monomial | polynomial
template<typename PolyRing, typename MonomialOrdering>
bool divides(Monomial<PolyRing> const &divisor, Polynomial<PolyRing, MonomialOrdering> const &dividend)
{
   for (size_t i = 0; i < dividend.terms(); ++i)
      if (!divides(divisor, dividend.getMonomial(i))) return false;
   return true;
}

// Relation: term | polynomial
template<typename PolyRing, typename MonomialOrdering>
bool divides(Term<PolyRing> const &divisor, Polynomial<PolyRing, MonomialOrdering> const &dividend)
{
   for (size_t i = 0; i < dividend.terms(); ++i)
      if (!divides(divisor, dividend[i])) return false;
   return true;
}


/// Division (binary operations)
////////////////////////////////////////////////////////////////////////////

// Divides 2 terms, assuming the dvisior divides the dividend (without a reminder).
template<typename PolyRing>
Term<PolyRing> safelyDivide(Term<PolyRing> const &divisor, Term<PolyRing> const &dividend)
{
   Term<PolyRing> res;
   for (size_t i = 0; i < PolyRing::VARIABLES; ++i)
      res.getMonomial().set(i, dividend.getMonomial()[i]-divisor.getMonomial()[i]);
   res.getCoeff() = dividend.getCoeff()/divisor.getCoeff();
   return res;
}

// Dvidies a polynomial with respect the a sequence of other polynomials.
// Returns a tuple with the reminder (a polynomial) and sequence of coefficients (also polynomials)
template<typename PolynomialType, typename DivisorsContainer>
std::tuple<PolynomialType, std::vector<PolynomialType>> divide(PolynomialType dividend, DivisorsContainer&& divisors)
{
   PolynomialType r;
   std::vector<PolynomialType> coeffs(divisors.size());

   while (dividend.terms() != 0)
   {
      bool division_occurred = false;
      for (size_t i = 0; (i < divisors.size()) && (!division_occurred); ++i)
      {
         auto &curr_divisor = *(divisors.begin()+i);
         if (divides(LT(curr_divisor), LT(dividend)))
         {
            auto d = safelyDivide(LT(curr_divisor), LT(dividend)); 
            coeffs[i] += d;
            dividend -= (d*curr_divisor);
            division_occurred = true;
         }
      }
      if (!division_occurred)
      {
         r += LT(dividend);
         dividend -= LT(dividend);
      }
   }
   return std::make_tuple(r, coeffs);
}


template<typename PolynomialType>
std::tuple<PolynomialType, std::vector<PolynomialType>> divide(PolynomialType dividend, std::initializer_list<PolynomialType> divisors)
{
   return divide(dividend, std::vector<PolynomialType>(divisors));
}

#endif

