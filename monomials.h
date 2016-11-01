// monomials.h

///////////////////////////////////////////////////////////////////////////////////////////
// (1) PolyRing is a policy describes the underlying polynomial ring K[x1, x2, ..., xn].
//     PolyRing::Coefficient - A type for the coefficients field (K).
//     PolyRing::VARAIBLES   - A constant. The number of varaibles of the polynomials (n).
// (2) class Monomial<PolyRing> is a multiplication of (powers of) varaibles.
///////////////////////////////////////////////////////////////////////////////////////////


#pragma once
#ifndef monomials_H__
#define monomials_H__

#include <cmath>
#include <array>
#include <string>
#include <sstream>
#include <cassert>
#include <numeric>
#include <algorithm>

// ** struct PolynomialRing
////////////////////////////////////////////////////////////////////////////

// Default policy for K[x_1, x_2, ..., x_n]. Assumes field-elements are native types.
template<typename FieldType, size_t VARIABLES_COUNT>
struct PolynomialRing
{
   using Coefficient = FieldType;
   static const size_t VARIABLES = VARIABLES_COUNT;
   static bool isZero(Coefficient a) {return std::fabs(a)<=1e-14;}; // Arbitrary choice.
};


// ** class Monomial
////////////////////////////////////////////////////////////////////////////
template <typename PolyRing>
class Monomial
{
public:
   Monomial();
   Monomial(std::array<unsigned int, PolyRing::VARIABLES> powers);

   std::string toString() const;

   bool operator==(Monomial<PolyRing> const& other) const;
   bool operator!=(Monomial<PolyRing> const& other) const;

   unsigned int operator[](size_t i) const;
   void set(size_t i, unsigned int v);
   unsigned int powersSum() const; // The sum of all the powers (used by some orderings)

private:
   std::array<unsigned int, PolyRing::VARIABLES> m_powers;
   unsigned int m_powers_sum;
};


// Least Common Multplication
template<typename PolyRing>
Monomial<PolyRing> LCM(Monomial<PolyRing> const &m1, Monomial<PolyRing> const &m2);


// Monomial ordering policies
/////////////////////////////////////////////////////////////////////////////

// Lexicographical Ordering
struct LexOrder
{
   template<typename PolyRing>
   static bool lessThen(Monomial<PolyRing> const &m1, Monomial<PolyRing> const &m2)
   {
      for (size_t i = 0; i < PolyRing::VARIABLES; ++i)
         if (m1[i] != m2[i]) return m1[i] < m2[i];
      return false;
   }
};

// Graded Lexicographical Ordering
struct GrlexOrder
{ 
   template<typename PolyRing>
   static bool lessThen(Monomial<PolyRing> const &m1, Monomial<PolyRing> const &m2)
   {
      if (m2.powersSum() != m1.powersSum())
         return (m2.powersSum() > m1.powersSum());
      return LexOrder::lessThen(m1, m2);
   }
};


// Graded Reverse Lexicographical Ordering
struct GrevlexOrder
{ 
   template<typename PolyRing>
   static bool lessThen(Monomial<PolyRing> const &m1, Monomial<PolyRing> const &m2)
   {
      if (m2.powersSum() != m1.powersSum())
         return (m2.powersSum() > m1.powersSum());
      for (int i = PolyRing::VARIABLES-1; i >= 0; --i)
         if (m1[i] != m2[i]) return m1[i] > m2[i];
      return false;
   }
};




// class Monomials - implementationa
/////////////////////////////////////////////////////////////////////////////
template <typename PolyRing>
Monomial<PolyRing>::Monomial()
   : m_powers_sum(0)
{
   std::fill(m_powers.begin(), m_powers.end(), 0);
}


template <typename PolyRing>
Monomial<PolyRing>::Monomial(std::array<unsigned int, PolyRing::VARIABLES> powers)
   : m_powers(powers), m_powers_sum(0)
{
   assert(powers.size() == PolyRing::VARIABLES);
   for (auto c: m_powers) m_powers_sum += c;
}


template <typename PolyRing>
std::string Monomial<PolyRing>::toString() const
{
   std::stringstream str;
   str << "<" << m_powers[0];
   for (size_t i = 1; i < PolyRing::VARIABLES; ++i)         
      str << "," << m_powers[i];
   str << ">";
   return str.str();
}

template <typename PolyRing>
bool Monomial<PolyRing>::operator==(Monomial<PolyRing> const& other) const
{
   for (size_t i = 0; i < PolyRing::VARIABLES; ++i)
      if (m_powers[i] != other[i]) return false;
   return true;
}

template <typename PolyRing>
bool Monomial<PolyRing>::operator!=(Monomial<PolyRing> const& other) const
{
   return !((*this) == other);
}

template <typename PolyRing>
unsigned int Monomial<PolyRing>::operator[](size_t i) const
{
   return m_powers[i];
}

template <typename PolyRing>
void Monomial<PolyRing>::set(size_t i, unsigned int v)
{
   m_powers_sum -= m_powers[i];
   m_powers_sum += v;
   m_powers[i] = v;
}

template <typename PolyRing>
unsigned int Monomial<PolyRing>::powersSum() const
{
   return m_powers_sum;
}

template<typename PolyRing>
Monomial<PolyRing> LCM(Monomial<PolyRing> const &m1, Monomial<PolyRing> const &m2)
{
   Monomial<PolyRing> res;
   for (size_t i = 0; i < PolyRing::VARIABLES; ++i)
      res.set(i, std::max(m1[i], m2[i]));
   return res;
}


#endif
