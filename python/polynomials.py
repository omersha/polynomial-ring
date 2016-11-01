import itertools
import numpy as np
import ctypes


class Polynomial(object):
    def __init__(self, terms, powers, coefficients):
        self._terms = terms
        self._powers = powers
        self._coefficients = coefficients
        
    def __repr__(self):
        return str(self._terms)
    
    def powers(self):
        return self._powers

    def coefficients(self):
        return self._coefficients
    
    def terms(self):
        return self._terms

class PolynomialRing(object):
    def __init__(self, sofile):
        self._lib = ctypes.cdll.LoadLibrary(sofile)
        # Addition
        self._lib.additionCtor.restype = ctypes.c_void_p
        self._lib.additionTerms.restype = ctypes.c_uint32
        self._lib.additionResult.restype = ctypes.c_uint32
        # Subtraction
        self._lib.subtractionCtor.restype = ctypes.c_void_p
        self._lib.subtractionTerms.restype = ctypes.c_uint32
        self._lib.subtractionResult.restype = ctypes.c_uint32        
        # Sparse Multiplication
        self._lib.sparseMultiplicationCtor.restype = ctypes.c_void_p
        self._lib.sparseMultiplicationTerms.restype = ctypes.c_uint32
        self._lib.sparseMultiplicationResult.restype = ctypes.c_uint32
        # Polynomial Division
        self._lib.divisionCtor.restype = ctypes.c_void_p
        self._lib.divisionCalculate.restype = ctypes.c_uint32
        self._lib.divisionQuotientTerms.restype = ctypes.c_uint32
        self._lib.divisionRemainderTerms.restype = ctypes.c_uint32
        self._lib.divisionQuotients.restype = ctypes.c_uint32
        self._lib.divisionQuotient.restype = ctypes.c_uint32
        self._lib.divisionRemainder.restype = ctypes.c_uint32        
        # Buchberger's Algorithm
        self._lib.buchbergersCtor.restype = ctypes.c_void_p
        self._lib.buchbergersBasisSize.restype = ctypes.c_uint32
        self._lib.buchbergersCalculate.restype = ctypes.c_uint32
        self._lib.buchbergersReduce.restype = ctypes.c_uint32
        self._lib.buchbergersMinimize.restype = ctypes.c_uint32
        self._lib.buchbergersBasisElementTerms.restype = ctypes.c_uint32
        self._lib.buchbergersBasisElement.restype = ctypes.c_uint32        
        
    def polynomial_from_terms(self, terms):
        powers = np.zeros((len(terms), 3)).astype(np.uint32)
        coefficients = np.zeros(len(terms)).astype(np.float64)
        for i in xrange(len(terms)):
            coefficients[i] = terms[i][0]
            powers[i, :] = [terms[i][1][0], terms[i][1][1], terms[i][1][2]]
        return Polynomial(terms, powers, coefficients)
            
    def polynomial_from_numpy(self, coefficients, powers):
        terms = [(coefficients[i], [powers[i, j] for j in xrange(powers.shape[1])]) for i in xrange(len(coefficients))]
        return Polynomial(terms, powers, coefficients)

    def add(self, *polynomials):
        handler = self._lib.additionCtor()
        for polynomial in polynomials:
            self._lib.additionAddSummand(ctypes.c_voidp(handler),
                                         ctypes.c_uint32(len(polynomial.coefficients())),
                                         polynomial.coefficients().ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                         polynomial.powers().ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)))
        terms = self._lib.additionTerms(ctypes.c_voidp(handler))
        out_coeffs = np.zeros(terms, dtype=np.float64)
        out_powers = np.zeros((terms, 3), dtype=np.uint32)
        self._lib.additionResult(ctypes.c_voidp(handler),
                                 out_coeffs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                 out_powers.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        self._lib.additionDtor(ctypes.c_voidp(handler))
        return self.polynomial_from_numpy(out_coeffs, out_powers)
 
    def sub(self, minuend, subtrahend):
        handler = self._lib.subtractionCtor(ctypes.c_uint32(len(minuend.coefficients())),
                                            minuend.coefficients().ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                            minuend.powers().ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                                            ctypes.c_uint32(len(subtrahend.coefficients())),
                                            subtrahend.coefficients().ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                            subtrahend.powers().ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)))
        terms = self._lib.subtractionTerms(ctypes.c_voidp(handler))
        out_coeffs = np.zeros(terms, dtype=np.float64)
        out_powers = np.zeros((terms, 3), dtype=np.uint32)
        self._lib.subtractionResult(ctypes.c_voidp(handler),
                                    out_coeffs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                    out_powers.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        self._lib.subtractionDtor(ctypes.c_voidp(handler))
        return self.polynomial_from_numpy(out_coeffs, out_powers)

    def mul(self, *polynomials):
        handler = self._lib.sparseMultiplicationCtor()
        for polynomial in polynomials:
            self._lib.sparseMultiplicationAddMultiplicand(ctypes.c_voidp(handler),
                                                         ctypes.c_uint32(len(polynomial.coefficients())),
                                                         polynomial.coefficients().ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                                         polynomial.powers().ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)))
        terms = self._lib.sparseMultiplicationTerms(ctypes.c_voidp(handler))
        out_coeffs = np.zeros(terms, dtype=np.float64)
        out_powers = np.zeros((terms, 3), dtype=np.uint32)
        self._lib.sparseMultiplicationResult(ctypes.c_voidp(handler),
                                            out_coeffs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                            out_powers.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        self._lib.sparseMultiplicationDtor(ctypes.c_voidp(handler))
        return self.polynomial_from_numpy(out_coeffs, out_powers)
    
    def div(self, dividend, divisors):
        handler = self._lib.divisionCtor(ctypes.c_uint32(len(dividend.coefficients())),
                                         dividend.coefficients().ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                         dividend.powers().ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)))
        for divisor in divisors:
            self._lib.divisionAddDivisor(ctypes.c_voidp(handler),
                                         ctypes.c_uint32(len(divisor.coefficients())),
                                         divisor.coefficients().ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                         divisor.powers().ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)))
        quotients_count = self._lib.divisionCalculate(ctypes.c_voidp(handler))
        
        quotients = []
        for i in xrange(quotients_count):
            terms = self._lib.divisionQuotientTerms(ctypes.c_voidp(handler), ctypes.c_uint32(i))
            if terms > 0:
                out_coeffs = np.zeros(terms, dtype=np.float64)
                out_powers = np.zeros((terms, 3), dtype=np.uint32)
                self._lib.divisionQuotient(ctypes.c_voidp(handler),
                                           ctypes.c_uint32(i),
                                           out_coeffs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                           out_powers.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
                quotients.append(self.polynomial_from_numpy(out_coeffs, out_powers))
        terms = self._lib.divisionRemainderTerms(ctypes.c_voidp(handler))
        out_coeffs = np.zeros(terms, dtype=np.float64)
        out_powers = np.zeros((terms, 3), dtype=np.uint32)
        self._lib.divisionRemainder(ctypes.c_voidp(handler),
                                    out_coeffs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                    out_powers.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        remainder = self.polynomial_from_numpy(out_coeffs, out_powers)
        self._lib.divisionDtor(ctypes.c_voidp(handler))
        return quotients, remainder
    
    
    def buchbergers(self, reduce_flag, *generators):
        handler = self._lib.buchbergersCtor()
        for generator in generators:
            self._lib.buchbergersAddGenerator(ctypes.c_voidp(handler),
                                              ctypes.c_uint32(len(generator.coefficients())),
                                              generator.coefficients().ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                              generator.powers().ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)))
        basis_size = self._lib.buchbergersCalculate(ctypes.c_voidp(handler))
        if reduce_flag:
            basis_size = self._lib.buchbergersReduce(ctypes.c_voidp(handler))
        
        groebner = []
        for i in xrange(basis_size):
            terms = self._lib.buchbergersBasisElementTerms(ctypes.c_voidp(handler),ctypes.c_uint32(i))
            out_coeffs = np.zeros(terms, dtype=np.float64)
            out_powers = np.zeros((terms, 3), dtype=np.uint32)
            self._lib.buchbergersBasisElement(ctypes.c_voidp(handler),
                                              ctypes.c_uint32(i),
                                              out_coeffs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                              out_powers.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
            groebner.append(self.polynomial_from_numpy(out_coeffs, out_powers))
        self._lib.buchbergersDtor(ctypes.c_voidp(handler))
        return groebner

    
    
def find_standard_monomial_basis(groebner_basis):
    LT = set([tuple(p.terms()[0][1]) for p in groebner_basis])
    return set([z for y in [list(itertools.product(xrange(i+1),xrange(j+1),xrange(k+1))) for i,j,k in LT] for z in y]).difference(LT)
