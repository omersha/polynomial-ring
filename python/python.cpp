#include "python.h"


extern "C"
{
   // Sparse Multiplication
   //////////////////////////////////////////////////////////////////////////

   void* sparseMultiplicationCtor()
   {
      return new SparseMultiplication<PythonPolyRing, PythonOrdering>();
   }

   void sparseMultiplicationDtor(void *handler)
   {
      delete static_cast<SparseMultiplication<PythonPolyRing, PythonOrdering>*>(handler);
   }

   void sparseMultiplicationAddMultiplicand(void *handler, unsigned int terms, double const * const coeffs, unsigned int const * const powers)
   {
      static_cast<SparseMultiplication<PythonPolyRing, PythonOrdering>*>(handler)->addMultiplicand(importPolynomial<PythonPolyRing, PythonOrdering>(terms, coeffs, powers));
   }

   unsigned int sparseMultiplicationTerms(void *handler)
   {
      return static_cast<SparseMultiplication<PythonPolyRing, PythonOrdering>*>(handler)->result().terms();
   }

   unsigned int sparseMultiplicationResult(void *handler, double * out_coeffs, unsigned int * out_powers)
   {
      return exportPolynomial(static_cast<SparseMultiplication<PythonPolyRing, PythonOrdering>*>(handler)->result(), out_coeffs, out_powers);
   }



   // Polynomial Division
   //////////////////////////////////////////////////////////////////////////
   void* divisionCtor(unsigned int terms, double const * const coeffs, unsigned int const * const powers)
   {
      return new Division<PythonPolyRing, PythonOrdering>(importPolynomial<PythonPolyRing, PythonOrdering>(terms, coeffs, powers));
   }

   void divisionDtor(void *handler)
   {
      delete static_cast<Division<PythonPolyRing, PythonOrdering>*>(handler);
   }

   void divisionAddDivisor(void *handler, unsigned int terms, double const * const coeffs, unsigned int const * const powers)
   {
      static_cast<Division<PythonPolyRing, PythonOrdering>*>(handler)->addDivisor(importPolynomial<PythonPolyRing, PythonOrdering>(terms, coeffs, powers));
   }

   unsigned int divisionCalculate(void *handler)
   {
      static_cast<Division<PythonPolyRing, PythonOrdering>*>(handler)->calculate();
      return static_cast<Division<PythonPolyRing, PythonOrdering>*>(handler)->quotients();
   }

   unsigned int divisionQuotients(void *handler)
   {
      return static_cast<Division<PythonPolyRing, PythonOrdering>*>(handler)->quotients();
   }

   unsigned int divisionQuotientTerms(void *handler, unsigned int i)
   {
      return static_cast<Division<PythonPolyRing, PythonOrdering>*>(handler)->quotient(i).terms();
   }

   unsigned int divisionQuotient(void *handler, unsigned int i, double * out_coeffs, unsigned int * out_powers)
   {
      return exportPolynomial(static_cast<Division<PythonPolyRing, PythonOrdering>*>(handler)->quotient(i), out_coeffs, out_powers);
   }

   unsigned int divisionRemainderTerms(void *handler)
   {
      return static_cast<Division<PythonPolyRing, PythonOrdering>*>(handler)->remainder().terms();
   }

   unsigned int divisionRemainder(void *handler, double * out_coeffs, unsigned int * out_powers)
   {
      return exportPolynomial(static_cast<Division<PythonPolyRing, PythonOrdering>*>(handler)->remainder(), out_coeffs, out_powers);
   }


   // Addition
   //////////////////////////////////////////////////////////////////////////

   void* additionCtor()
   {
      return new Addition<PythonPolyRing, PythonOrdering>();
   }

   void additionDtor(void *handler)
   {
      delete static_cast<Addition<PythonPolyRing, PythonOrdering>*>(handler);
   }

   void additionAddSummand(void *handler, unsigned int terms, double const * const coeffs, unsigned int const * const powers)
   {
      static_cast<Addition<PythonPolyRing, PythonOrdering>*>(handler)->addSummand(importPolynomial<PythonPolyRing, PythonOrdering>(terms, coeffs, powers));
   }

   unsigned int additionTerms(void *handler)
   {
      return static_cast<Addition<PythonPolyRing, PythonOrdering>*>(handler)->result().terms();
   }

   unsigned int additionResult(void *handler, double * out_coeffs, unsigned int * out_powers)
   {
      return exportPolynomial(static_cast<Addition<PythonPolyRing, PythonOrdering>*>(handler)->result(), out_coeffs, out_powers);
   }


   // Subtraction
   //////////////////////////////////////////////////////////////////////////

   void* subtractionCtor(unsigned int subtrahend_terms, double const * const subtrahend_coeffs, unsigned int const * const subtrahend_powers,
                         unsigned int minuend_terms, double const * const minuend_coeffs, unsigned int const * const minuend_powers)
   {
      return new Subtraction<PythonPolyRing, PythonOrdering>(importPolynomial<PythonPolyRing, PythonOrdering>(subtrahend_terms, subtrahend_coeffs, subtrahend_powers),
                                                             importPolynomial<PythonPolyRing, PythonOrdering>(minuend_terms, minuend_coeffs, minuend_powers));
   }

   void subtractionDtor(void *handler)
   {
      delete static_cast<Subtraction<PythonPolyRing, PythonOrdering>*>(handler);
   }

   unsigned int subtractionTerms(void *handler)
   {
      return static_cast<Subtraction<PythonPolyRing, PythonOrdering>*>(handler)->result().terms();
   }

   unsigned int subtractionResult(void *handler, double * out_coeffs, unsigned int * out_powers)
   {
      return exportPolynomial(static_cast<Subtraction<PythonPolyRing, PythonOrdering>*>(handler)->result(), out_coeffs, out_powers);
   }



   // Buchbergers
   //////////////////////////////////////////////////////////////////////////
   void* buchbergersCtor()
   {
      return new Buchbergers<PythonPolyRing, PythonOrdering>();
   }

   void buchbergersDtor(void *handler)
   {
      delete static_cast<Buchbergers<PythonPolyRing, PythonOrdering>*>(handler);
   }

   void buchbergersAddGenerator(void *handler, unsigned int terms, double const * const coeffs, unsigned int const * const powers)
   {
      static_cast<Buchbergers<PythonPolyRing, PythonOrdering>*>(handler)->addIdealGenerator(importPolynomial<PythonPolyRing, PythonOrdering>(terms, coeffs, powers));
   }

   unsigned int buchbergersCalculate(void *handler)
   {
      static_cast<Buchbergers<PythonPolyRing, PythonOrdering>*>(handler)->calculate();
      return static_cast<Buchbergers<PythonPolyRing, PythonOrdering>*>(handler)->basisSize();
   }

   unsigned int buchbergersBasisSize(void *handler)
   {
      return static_cast<Buchbergers<PythonPolyRing, PythonOrdering>*>(handler)->basisSize();
   }

   unsigned int buchbergersBasisElementTerms(void *handler, unsigned int i)
   {
      return static_cast<Buchbergers<PythonPolyRing, PythonOrdering>*>(handler)->basisElement(i).terms();
   }

   unsigned int buchbergersBasisElement(void *handler, unsigned int i, double * out_coeffs, unsigned int * out_powers)
   {
      return exportPolynomial(static_cast<Buchbergers<PythonPolyRing, PythonOrdering>*>(handler)->basisElement(i), out_coeffs, out_powers);
   }



} // extern "C"

