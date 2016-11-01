// python.h

#pragma once
#ifndef python_H__
#define python_H__

#include "monomials.h"
#include "polynomials.h"
#include "division.h"
#include "buchbergers.h"


using PythonPolyRing = PolynomialRing<double, 3>;
using PythonOrdering = LexOrder;
//using PythonOrdering = GrevlexOrder;


template<typename Ctor, size_t... S>
std::array<std::result_of_t<Ctor(size_t)>, sizeof...(S)> makeArray(Ctor&& ctor, std::index_sequence<S...>)
{
   return std::array<std::result_of_t<Ctor(size_t)>, sizeof...(S)> {ctor(S)...};
}

template<size_t N, typename Ctor>
std::array<std::result_of_t<Ctor(size_t)>, N> makeArray(Ctor&& ctor)
{
   return makeArray(std::forward<Ctor>(ctor), std::make_index_sequence<N>());
}



// Import / Export
//////////////////////////////////////////////////////////////////////////


template<typename PolyRing, class MonomialOrdering>
Polynomial<PolyRing, MonomialOrdering> importPolynomial(unsigned int terms, double const * const coeffs, unsigned int const * const powers)
{
   Polynomial<PolyRing, PythonOrdering> p;
   for (unsigned int i = 0; i < terms; ++i) {
      p += Term<PolyRing>(coeffs[i], makeArray<PolyRing::VARIABLES>([i, powers](size_t j){return powers[i*PolyRing::VARIABLES+j];}));
   }
   return p;
}

template<typename PolyRing, class MonomialOrdering>
unsigned int exportPolynomial(Polynomial<PolyRing, MonomialOrdering> polynomial, double * out_coeffs, unsigned int * out_powers)
{
   for (unsigned int i = 0; i < polynomial.terms(); ++i)
   {
      out_coeffs[i] = polynomial[i].getCoeff();
      auto& monomial = polynomial[i].getMonomial();
      for (unsigned int j = 0; j < PolyRing::VARIABLES; ++j) {
         out_powers[PolyRing::VARIABLES*i+j] = monomial[j];
      }
   }
   return polynomial.terms();
}


// SparseMultiplication
//////////////////////////////////////////////////////////////////////////
template<typename PolyRing, class MonomialOrdering>
class SparseMultiplication
{
public:
   void addMultiplicand(Polynomial<PolyRing, MonomialOrdering> multiplicand)
   {
      if (m_result.terms() == 0) {
         m_result = multiplicand;
      } else {
         m_result = m_result * multiplicand;
      }
   }

   Polynomial<PolyRing, MonomialOrdering> result()
   {
      return m_result;
   }

private:
   Polynomial<PolyRing, MonomialOrdering> m_result;
}; // SparseMultiplication


// Division
//////////////////////////////////////////////////////////////////////////
template<typename PolyRing, class MonomialOrdering>
class Division
{
public:
   Division(Polynomial<PolyRing, MonomialOrdering> dividend)
      : m_dividend(dividend) {}

   void addDivisor(Polynomial<PolyRing, MonomialOrdering> divisor)
   {
      m_divisors.push_back(divisor);
   }

   void calculate()
   {
      std::tie(m_remainder, m_quotients) = divide(m_dividend, m_divisors);
   }

   size_t quotients()
   {
      return m_quotients.size();
   }

   Polynomial<PolyRing, MonomialOrdering> quotient(size_t i)
   {
      return m_quotients[i];
   }

   Polynomial<PolyRing, MonomialOrdering> remainder()
   {
      return m_remainder;
   }

private:
   Polynomial<PolyRing, MonomialOrdering> m_dividend;
   std::deque<Polynomial<PolyRing, MonomialOrdering>> m_divisors;
   
   Polynomial<PolyRing, MonomialOrdering> m_remainder;
   std::vector<Polynomial<PolyRing, MonomialOrdering>> m_quotients;
}; // Division


// Addition
//////////////////////////////////////////////////////////////////////////
template<typename PolyRing, class MonomialOrdering>
class Addition
{
public:
   void addSummand(Polynomial<PolyRing, MonomialOrdering> summand)
   {
      m_result += summand;
   }

   Polynomial<PolyRing, MonomialOrdering> result()
   {
      return m_result;
   }

private:
   Polynomial<PolyRing, MonomialOrdering> m_result;
}; // Addition


// Subtraction
//////////////////////////////////////////////////////////////////////////
template<typename PolyRing, class MonomialOrdering>
class Subtraction
{
public:
   Subtraction(Polynomial<PolyRing, MonomialOrdering> minuend, Polynomial<PolyRing, MonomialOrdering> subtrahend)
      : m_result(minuend) {m_result -= subtrahend;}

   Polynomial<PolyRing, MonomialOrdering> result()
   {
      return m_result;
   }

private:
   Polynomial<PolyRing, MonomialOrdering> m_result;
}; // Subtraction




// Buchbergers
//////////////////////////////////////////////////////////////////////////
template<typename PolyRing, class MonomialOrdering>
class Buchbergers
{
public:
   Buchbergers() : m_minimal(false) {}

   size_t generators()
   {
      return m_ideal_generators.size();
   }

   void addIdealGenerator(Polynomial<PolyRing, MonomialOrdering> polynomial)
   {
      m_ideal_generators.push_back(polynomial);
   }

   void calculate()
   {
      m_minimal = false;
      m_groebner = runBuchbergers(m_ideal_generators);
   }

   void reduce()
   {
      if (!m_minimal) {
         minimize();
      }
      makeReducedGroebner(m_groebner);
   }

   void minimize()
   {
      makeMinimalGroebner(m_groebner);
      m_minimal = true;
   }

   size_t basisSize()
   {
      return m_groebner.size();
   }

   Polynomial<PolyRing, MonomialOrdering> basisElement(size_t i)
   {
      return m_groebner[i];
   }

private:
   bool m_minimal;
   std::deque<Polynomial<PolyRing, MonomialOrdering>> m_ideal_generators;
   std::deque<Polynomial<PolyRing, MonomialOrdering>> m_groebner;
};



#endif
