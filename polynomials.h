// polynomials.h

///////////////////////////////////////////////////////////////////////////////////////////
// (1) class Term<PolyRing> is a  polynomial term (coefficient from the field + a monomial).
// (2) class Polynomial<PolyRing, MonomialOrdering> is an ordered sequence of terms.
// Methods:
//   * LC:polymonial-->coefficient             : Leading coefficient.
//   * LM:polymonial-->monomial                : Leading monomial.
//   * LT:polymonial-->term                    : Leading term.
//   * operator* (coefficient, polynomial)     : Multiplication.
//   * operator* (term, polynomial)            : Multiplication.
//   * operator* (polynomial, term)            : Multiplication.
//   * operator* (polynomial, polynomial)      : Multiplication.

///////////////////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef polynomials_H__
#define polynomials_H__

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <initializer_list>

#include "monomials.h"


// Terms
////////////////////////////////////////////////////////////////////////////

// ** class Term
template<typename PolyRing>
class Term
{
public:
   Term();
   Term(typename PolyRing::Coefficient c, Monomial<PolyRing> m);
   Term(std::tuple<typename PolyRing::Coefficient, Monomial<PolyRing>> tup);

   std::string toString() const;
   
   bool operator==(Term<PolyRing> const &other) const;

   Monomial<PolyRing>& getMonomial();
   Monomial<PolyRing> const& getMonomial() const;

   typename PolyRing::Coefficient& getCoeff();
   typename PolyRing::Coefficient const& getCoeff() const;

   void operator*=(typename PolyRing::Coefficient factor);
   void operator*=(Term<PolyRing> const &other);

private:
   typename PolyRing::Coefficient m_coeff;
   Monomial<PolyRing> m_monomial;
};

template<typename PolyRing>
Term<PolyRing> operator*(typename PolyRing::Coefficient, Term<PolyRing> term);


// Polynomials
////////////////////////////////////////////////////////////////////////////

// ** class Polynomial
template<typename PolyRing, typename MonomialOrdering>
class Polynomial
{
public:
   typedef Term<PolyRing> TermType;

   Polynomial(size_t terms_preallocation=10);
   Polynomial(std::initializer_list<Term<PolyRing>> terms);

   std::string toString() const;

   bool operator==(Polynomial<PolyRing, MonomialOrdering> const &other) const;
   bool operator!=(Polynomial<PolyRing, MonomialOrdering> const &other) const;

   void operator+=(TermType term);
   void operator+=(Polynomial<PolyRing, MonomialOrdering> polynomial);
   void operator-=(TermType term);
   void operator-=(Polynomial<PolyRing, MonomialOrdering> polynomial);
   void operator*=(typename PolyRing::Coefficient factor);
   void operator*=(TermType const &m);
 
   size_t terms() const;
   TermType const& operator[](size_t i) const;
   typename PolyRing::Coefficient const& getCoeff(size_t i) const;
   Monomial<PolyRing> const& getMonomial(size_t i) const;
 
   void normalize(); // Factors so the leading coefficient is 1.

private:
   void collectTerms(); 
   void removeZeros(); // Removes terms whose coefficient is 0.
   void sortSelf();
   
private:
   std::vector<TermType> m_terms;
};

// Leading coefficient
template<typename PolyRing, typename MonomialOrdering>
typename PolyRing::Coefficient LC(Polynomial<PolyRing, MonomialOrdering> const &p);

// Leading monomial
template<typename PolyRing, typename MonomialOrdering>
Monomial<PolyRing> const& LM(Polynomial<PolyRing, MonomialOrdering> const &p);

// Leading term
template<typename PolynomialType>
typename PolynomialType::TermType const& LT(PolynomialType const &p);


template<typename PolyRing, typename MonomialOrdering>
Polynomial<PolyRing, MonomialOrdering> operator*(Term<PolyRing> const &m, Polynomial<PolyRing, MonomialOrdering> p);

template<typename PolyRing, typename MonomialOrdering>
Polynomial<PolyRing, MonomialOrdering> operator*(Polynomial<PolyRing, MonomialOrdering> p, Term<PolyRing> const &m);

template<typename PolyRing, typename MonomialOrdering>
Polynomial<PolyRing, MonomialOrdering> operator*(typename PolyRing::Coefficient factor, Polynomial<PolyRing, MonomialOrdering> p);

template<typename PolyRing, typename MonomialOrdering>
Polynomial<PolyRing, MonomialOrdering> operator*(Polynomial<PolyRing, MonomialOrdering> const &p1, Polynomial<PolyRing, MonomialOrdering> const &p2);


// Terms - Implementation
////////////////////////////////////////////////////////////////////////////

template<typename PolyRing>
Term<PolyRing>::Term() {}

template<typename PolyRing>
Term<PolyRing>::Term(typename PolyRing::Coefficient c, Monomial<PolyRing> m)
   : m_coeff(c), m_monomial(std::move(m)) {}

template<typename PolyRing>
Term<PolyRing>::Term(std::tuple<typename PolyRing::Coefficient, Monomial<PolyRing>> tup)
   : m_coeff(std::get<0>(tup)), m_monomial(std::move(std::get<1>(tup))) {}

template<typename PolyRing>
std::string Term<PolyRing>::toString() const
{
   std::stringstream str;
   str << getCoeff() << "*" << getMonomial().toString();
   return str.str();
}

template<typename PolyRing>
bool Term<PolyRing>::operator==(Term<PolyRing> const &other) const
{
   return ((getCoeff() == other.getCoeff()) && (getMonomial() == other.getMonomial()));
}

template<typename PolyRing>
Monomial<PolyRing>& Term<PolyRing>::getMonomial()
{
   return m_monomial;
}

template<typename PolyRing>
Monomial<PolyRing> const& Term<PolyRing>::getMonomial() const
{
   return m_monomial;
}

template<typename PolyRing>
typename PolyRing::Coefficient& Term<PolyRing>::getCoeff()
{
   return m_coeff;
}

template<typename PolyRing>
typename PolyRing::Coefficient const& Term<PolyRing>::getCoeff() const
{
   return m_coeff;
}

template<typename PolyRing>
void Term<PolyRing>::operator*=(typename PolyRing::Coefficient factor)
{
   m_coeff *= factor;
}

template<typename PolyRing>
void Term<PolyRing>::operator*=(Term<PolyRing> const &other)
{
   for (size_t i = 0; i < PolyRing::VARIABLES; ++i)
      m_monomial.set(i, m_monomial[i]+other.getMonomial()[i]);
   m_coeff *= other.getCoeff();
}

template<typename PolyRing>
Term<PolyRing> operator*(typename PolyRing::Coefficient factor, Term<PolyRing> term)
{
   term *= factor;
   return term;
}



// Polynomials - Implementation
////////////////////////////////////////////////////////////////////////////

template<typename PolyRing, typename MonomialOrdering>
Polynomial<PolyRing, MonomialOrdering>::Polynomial(size_t terms_preallocation) 
{
   m_terms.reserve(terms_preallocation);
}

template<typename PolyRing, typename MonomialOrdering>
Polynomial<PolyRing, MonomialOrdering>::Polynomial(std::initializer_list<Term<PolyRing>> terms)
   : m_terms(terms)
{
   sortSelf();
}

template<typename PolyRing, typename MonomialOrdering>
std::string Polynomial<PolyRing, MonomialOrdering>::toString() const
{
   if (terms() == 0) return "";
   
   std::stringstream str;
   str << getCoeff(0) << "*" << getMonomial(0).toString();
   for (size_t i = 1; i < terms(); ++i)
      str << " + " << getCoeff(i) << "*" << getMonomial(i).toString();
   return str.str();
}

template<typename PolyRing, typename MonomialOrdering>
bool Polynomial<PolyRing, MonomialOrdering>::operator==(Polynomial<PolyRing, MonomialOrdering> const &other) const
{
   if (terms() != other.terms()) return false;
   for (size_t i = 0; i < terms(); ++i)
      if (getMonomial(i) != other.getMonomial(i)) return false;
   return true;
}

template<typename PolyRing, typename MonomialOrdering>
bool Polynomial<PolyRing, MonomialOrdering>::operator!=(Polynomial<PolyRing, MonomialOrdering> const &other) const
{
   return !((*this) == other);
}

template<typename PolyRing, typename MonomialOrdering>
void Polynomial<PolyRing, MonomialOrdering>::operator+=(TermType term)
{
   m_terms.push_back(term);
   sortSelf();
}

template<typename PolyRing, typename MonomialOrdering>
void Polynomial<PolyRing, MonomialOrdering>::operator+=(Polynomial<PolyRing, MonomialOrdering> polynomial)
{
   for (size_t i = 0; i < polynomial.terms(); ++i)
      m_terms.push_back(polynomial[i]);
   sortSelf();
}

template<typename PolyRing, typename MonomialOrdering>
void Polynomial<PolyRing, MonomialOrdering>::operator-=(TermType term)
{
   (*this) += (-1*term);
}

template<typename PolyRing, typename MonomialOrdering>
void Polynomial<PolyRing, MonomialOrdering>::operator-=(Polynomial<PolyRing, MonomialOrdering> polynomial)
{
   for (size_t i = 0; i < polynomial.terms(); ++i)
      m_terms.push_back(-1*polynomial[i]);
   sortSelf();
}

template<typename PolyRing, typename MonomialOrdering>
void Polynomial<PolyRing, MonomialOrdering>::operator*=(typename PolyRing::Coefficient factor)
{
   for (auto &t: m_terms) t *= factor;
   sortSelf();
}

template<typename PolyRing, typename MonomialOrdering>
void Polynomial<PolyRing, MonomialOrdering>::operator*=(TermType const &m)
{
   for (auto &t: m_terms) t *= m;
   sortSelf();
}

template<typename PolyRing, typename MonomialOrdering>
size_t Polynomial<PolyRing, MonomialOrdering>::terms() const
{
   return m_terms.size();
}

template<typename PolyRing, typename MonomialOrdering>
typename Polynomial<PolyRing, MonomialOrdering>::TermType const& Polynomial<PolyRing, MonomialOrdering>::operator[](size_t i) const
{
   return m_terms[i];
}

template<typename PolyRing, typename MonomialOrdering>
typename PolyRing::Coefficient const& Polynomial<PolyRing, MonomialOrdering>::getCoeff(size_t i) const
{
   return m_terms[i].getCoeff();
}

template<typename PolyRing, typename MonomialOrdering>
Monomial<PolyRing> const& Polynomial<PolyRing, MonomialOrdering>::getMonomial(size_t i) const
{
   return m_terms[i].getMonomial();
}

template<typename PolyRing, typename MonomialOrdering>
void Polynomial<PolyRing, MonomialOrdering>::normalize()
{
   auto factor = m_terms[0].getCoeff();
   for (auto &t: m_terms)
      t.getCoeff() /= factor;
}

template<typename PolyRing, typename MonomialOrdering>
void Polynomial<PolyRing, MonomialOrdering>::collectTerms()
{
   // Assumes *this is sorted.
   auto prev = m_terms.begin();
   auto it = m_terms.begin()+1;
   while (it != m_terms.end())
   {
      if (it->getMonomial() == prev->getMonomial())
      {
         prev->getCoeff() += it->getCoeff();
         it = m_terms.erase(it);
         prev = it-1;
      }
      else
      {
        prev = it;
         ++it;
      }
  }
}

template<typename PolyRing, typename MonomialOrdering>
void Polynomial<PolyRing, MonomialOrdering>::removeZeros()
{
   auto it = m_terms.begin();
   while (it != m_terms.end())
   {
      if (PolyRing::isZero(it->getCoeff()))
         it = m_terms.erase(it);
      else ++it;
   }
}

template<typename PolyRing, typename MonomialOrdering>
void Polynomial<PolyRing, MonomialOrdering>::sortSelf()
{
   if (terms() == 0) return;
   std::sort(m_terms.rbegin(), m_terms.rend(),
             [](TermType const &m1, TermType const &m2) {return MonomialOrdering::lessThen(m1.getMonomial(), m2.getMonomial());});
   collectTerms();
   removeZeros();
}


template<typename PolyRing, typename MonomialOrdering>
typename PolyRing::Coefficient LC(Polynomial<PolyRing, MonomialOrdering> const &p)
{
   return p.getCoeff(0);
}


template<typename PolyRing, typename MonomialOrdering>
Monomial<PolyRing> const& LM(Polynomial<PolyRing, MonomialOrdering> const &p)
{
   return p.getMonomial(0);
}

template<typename PolynomialType>
typename PolynomialType::TermType const& LT(PolynomialType const &p)
{
   return p[0];
}

template<typename PolyRing, typename MonomialOrdering>
Polynomial<PolyRing, MonomialOrdering> operator*(Term<PolyRing> const &m, Polynomial<PolyRing, MonomialOrdering> p)
{
   Polynomial<PolyRing, MonomialOrdering> res(std::move(p));
   res *= m;
   return res;
}

template<typename PolyRing, typename MonomialOrdering>
Polynomial<PolyRing, MonomialOrdering> operator*(Polynomial<PolyRing, MonomialOrdering> p, Term<PolyRing> const &m)
{
   return m*p;
}

template<typename PolyRing, typename MonomialOrdering>
Polynomial<PolyRing, MonomialOrdering> operator*(typename PolyRing::Coefficient factor, Polynomial<PolyRing, MonomialOrdering> p)
{
   Polynomial<PolyRing, MonomialOrdering> res(std::move(p));
   res *= factor;
   return res;
}



// Multiplication of 2 polynomials (direct).
template<typename PolyRing, typename MonomialOrdering>
Polynomial<PolyRing, MonomialOrdering> operator*(Polynomial<PolyRing, MonomialOrdering> const &p1, Polynomial<PolyRing, MonomialOrdering> const &p2)
{
   Polynomial<PolyRing, MonomialOrdering> product;
   for (size_t i = 0; i < p1.terms(); ++i) {
         product += p2*p1[i];
   }
   return product;
}



#endif
