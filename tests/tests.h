// tests.h

#pragma once
#ifndef tests_H__
#define tests_H__

#include <cassert>

#include "monomials.h"
#include "polynomials.h"
#include "division.h"
#include "buchbergers.h"

namespace Tests
{
   using PolyRing2 = PolynomialRing<double, 2>;
   using PolyRing3 = PolynomialRing<double, 3>;
   using PolyRing4 = PolynomialRing<double, 4>;

   template<class PolyRing>
   Monomial<PolyRing> X(size_t i) {Monomial<PolyRing> x; x.set(i, 1); return x;}
   
   void testMonomial()
   {  
      Monomial<PolyRing3> x1({1,0,0});
      Monomial<PolyRing3> x2({0,1,0});
      Monomial<PolyRing3> x3({0,0,1});
      Monomial<PolyRing3> m1({1,2,0});
      Monomial<PolyRing3> m2({0,3,4});
      Monomial<PolyRing3> m3({3,2,4});
      Monomial<PolyRing3> m4({3,2,1});

      assert(x1.powersSum() == 1);
      assert(x2.powersSum() == 1);
      assert(x3.powersSum() == 1);
      assert(m1.powersSum() == 3);
      assert(m2.powersSum() == 7);
      assert(m3.powersSum() == 9);
      assert(m4.powersSum() == 6);
   }

   void testLex()
   {
      Monomial<PolyRing3> x1({1,0,0});
      Monomial<PolyRing3> x2({0,1,0});
      Monomial<PolyRing3> x3({0,0,1});
      Monomial<PolyRing3> m1({1,2,0});
      Monomial<PolyRing3> m2({0,3,4});
      Monomial<PolyRing3> m3({3,2,4});
      Monomial<PolyRing3> m4({3,2,1});

      assert(LexOrder::lessThen(x3, x2));
      assert(LexOrder::lessThen(x2, x1));
      assert(!LexOrder::lessThen(x1, x1));
      assert(!LexOrder::lessThen(x3, x3));

      assert(LexOrder::lessThen(m2, m1));
      assert(!LexOrder::lessThen(m1, m2));

      assert(LexOrder::lessThen(m4, m3));
      assert(!LexOrder::lessThen(m3, m4));
   }

   void testGrlex()
   {
      Monomial<PolyRing3> x1({1,0,0});
      Monomial<PolyRing3> x2({0,1,0});
      Monomial<PolyRing3> x3({0,0,1});
      Monomial<PolyRing3> m1({1,2,3});
      Monomial<PolyRing3> m2({3,2,0});
      Monomial<PolyRing3> m3({1,2,4});
      Monomial<PolyRing3> m4({1,1,5});

      assert(GrlexOrder::lessThen(x3, x2));
      assert(GrlexOrder::lessThen(x2, x1));
      assert(!GrlexOrder::lessThen(x1, x1));
      assert(!GrlexOrder::lessThen(x3, x3));

      assert(GrlexOrder::lessThen(m2, m1));
      assert(!GrlexOrder::lessThen(m1, m2));

      assert(GrlexOrder::lessThen(m4, m3));
      assert(!GrlexOrder::lessThen(m3, m4));
   }

   void testPolynomial()
   {
      Polynomial<PolyRing3, LexOrder> p1;
      p1 += Term<PolyRing3>(2, Monomial<PolyRing3>({0, 1, 0}));
      p1 += Term<PolyRing3>(3, {{0, 0, 1}});
      p1 += Term<PolyRing3>(1, {{1, 0, 0}});

      for (size_t i = 0; i < p1.terms(); ++i)
      {
         assert(p1.getCoeff(i) == i+1);
         assert(p1.getMonomial(i)[i] == 1);
      }
      assert(LC(p1) == 1);
      assert(LM(p1)[0] == 1);
      assert(LT(p1).getCoeff() == LC(p1));

      Polynomial<PolyRing3, LexOrder> p2;
      p2 += p1;
      assert(p1 == p2);
      p2 += Term<PolyRing3>(-1, {{1, 0, 0}});
      p2 -= Term<PolyRing3>(2, {{0, 1, 0}});
      assert(p2.terms() == 1);
      assert(p2.getMonomial(0) == p1.getMonomial(2));
      assert(p2.getCoeff(0) == p1.getCoeff(2));
      p2 += Term<PolyRing3>(-3, {{0, 0, 1}});
      assert(p2.terms() == 0);
   }

   void testDivides()
   { 
      Monomial<PolyRing3> x1({1,0,0});
      Monomial<PolyRing3> x2({0,1,0});
      Monomial<PolyRing3> x3({0,0,1});

      Monomial<PolyRing3> m1({1,1,1});
      Monomial<PolyRing3> m2({4,2,1});
      Monomial<PolyRing3> m3({2,2,0});

      Polynomial<PolyRing3, LexOrder> p1;
      p1 += Term<PolyRing3>(1, m1);
      p1 += Term<PolyRing3>(2, m2);

      assert(!divides(x1, x2));
      assert(!divides(x2, x3));
      assert(!divides(x3, x1));

      assert(divides(x1, m1));
      assert(divides(x2, m1));
      assert(divides(x3, m1));

      assert(!divides(m1, x1));
      assert(!divides(m1, x2));
      assert(!divides(m1, x3));
 
      assert(divides(m1, m1));

      assert(divides(m1, m2));
      assert(!divides(m2, m1));
      assert(divides(m3, m2));
      assert(!divides(m2, m3));

      assert(divides(x1, p1));
      assert(divides(x2, p1));
      assert(divides(x3, p1));
      assert(divides(m1, p1));
      assert(!divides(m2, p1));
      assert(!divides(m3, p1));
   }

   void testDivision()
   {
      {
         Polynomial<PolyRing2, LexOrder> f;
         f += Term<PolyRing2>(1, {{2, 1}});
         f += Term<PolyRing2>(1, {{1, 2}});
         f += Term<PolyRing2>(1, {{0, 2}});

         Polynomial<PolyRing2, LexOrder> f1;
         f1 += Term<PolyRing2>(1, {{0, 2}});
         f1 += Term<PolyRing2>(-1, {{0, 0}});

         Polynomial<PolyRing2, LexOrder> f2;
         f2 += Term<PolyRing2>(1, {{1, 1}});
         f2 += Term<PolyRing2>(-1, {{0, 0}});

         auto res1 = divide(f, {f1, f2});

         // r = 2x+1
         assert(std::get<0>(res1).terms() == 2);
         assert(std::get<0>(res1)[0] == Term<PolyRing2>(2, {{1,0}}));
         assert(std::get<0>(res1)[1] == Term<PolyRing2>(1, {{0,0}}));
          // a1 = x+1
         assert(std::get<1>(res1)[0].terms() == 2);
         assert(std::get<1>(res1)[0][0] == Term<PolyRing2>(1, {{1,0}}));
         assert(std::get<1>(res1)[0][1] == Term<PolyRing2>(1, {{0,0}}));
         // a2 = x
         assert(std::get<1>(res1)[1].terms() == 1);
         assert(std::get<1>(res1)[1][0] == Term<PolyRing2>(1, {{1,0}}));
  
  
         auto res2 = divide(f, {f2, f1}); // 80

         // r = x+y+1
         assert(std::get<0>(res2).terms() == 3);
         assert(std::get<0>(res2)[0] == Term<PolyRing2>(1, {{1,0}}));
         assert(std::get<0>(res2)[1] == Term<PolyRing2>(1, {{0,1}}));
         assert(std::get<0>(res2)[2] == Term<PolyRing2>(1, {{0,0}}));
          // a1 = x+y
         assert(std::get<1>(res2)[0].terms() == 2);
         assert(std::get<1>(res2)[0][0] == Term<PolyRing2>(1, {{1,0}}));
         assert(std::get<1>(res2)[0][1] == Term<PolyRing2>(1, {{0,1}}));
         // a2 = 1
         assert(std::get<1>(res2)[1].terms() == 1);
         assert(std::get<1>(res2)[1][0] == Term<PolyRing2>(1, {{0,0}}));
     }
   }

   void testSPoly()
   {
      // f = x^3*y^2 - x^2*y^3 + x
      Polynomial<PolyRing2, GrlexOrder> f;
      f += Term<PolyRing2>(1, {{3, 2}});
      f += Term<PolyRing2>(-1, {{2, 3}});
      f += Term<PolyRing2>(1, {{1, 0}});

      // g = 3*x^4*y + y^2
      Polynomial<PolyRing2, GrlexOrder> g;
      g += Term<PolyRing2>(3, {{4, 1}});
      g += Term<PolyRing2>(1, {{0, 2}});

      // s = -x^3*y^3 - (1/3)*y^3 + x^2
      auto s = makeSPolynomial(f, g);
      assert(s.terms() == 3);
      assert(s[0] == Term<PolyRing2>(-1, {{3,3}}));
      assert(s[1] == Term<PolyRing2>(-1/3.0, {{0,3}}));
      assert(s[2] == Term<PolyRing2>(1, {{2,0}}));
   }

   void testBuchbergers1()
   {
      // f = x^3 - 2xy
      Polynomial<PolyRing2, GrlexOrder> f;
      f += Term<PolyRing2>(1, {{3, 0}});
      f += Term<PolyRing2>(-2, {{1, 1}});
      
      // g = x^2y - 2y^2 + x
      Polynomial<PolyRing2, GrlexOrder> g;
      g += Term<PolyRing2>(1, {{2, 1}});
      g += Term<PolyRing2>(-2, {{0, 2}});
      g += Term<PolyRing2>(1, {{1, 0}});
   
      
      // f1 = f
      Polynomial<PolyRing2, GrlexOrder> f1(f);

      // f2 = g
      Polynomial<PolyRing2, GrlexOrder> f2(g);
 
      // f3 = -x^2
      Polynomial<PolyRing2, GrlexOrder> f3;
      f3 += Term<PolyRing2>(-1, {{2, 0}});
   
      // f4 = -2xy
      Polynomial<PolyRing2, GrlexOrder> f4;
      f4 += Term<PolyRing2>(-2, {{1, 1}});
   
      // f5 = -2y^2+x
      Polynomial<PolyRing2, GrlexOrder> f5;
      f5 += Term<PolyRing2>(-2, {{0, 2}});
      f5 += Term<PolyRing2>(1, {{1, 0}});
   
      auto groebner = runBuchbergers({f, g});
      assert(groebner[0] == f1);
      assert(groebner[1] == f2);
      assert(groebner[2] == f3);
      assert(groebner[3] == f4);
      assert(groebner[4] == f5);

      // f1_reduced = x^2
      Polynomial<PolyRing2, GrlexOrder> f1_reduced;
      f1_reduced += Term<PolyRing2>(1, {{2, 0}});

      // f2_reduced = xy
      Polynomial<PolyRing2, GrlexOrder> f2_reduced;
      f2_reduced += Term<PolyRing2>(1, {{1, 1}});

      // f3_reduced = y^2 - 0.5*x
      Polynomial<PolyRing2, GrlexOrder> f3_reduced;
      f3_reduced += Term<PolyRing2>(1, {{0, 2}});
      f3_reduced += Term<PolyRing2>(-1, {{1, 0}});

      makeMinimalGroebner(groebner);
      makeReducedGroebner(groebner);

      assert(groebner.size() == 3);
      assert(groebner[0] == f1_reduced);
      assert(groebner[1] == f2_reduced);
      assert(groebner[2] == f3_reduced);     
  }

   void testBuchbergers2()
   {
      Polynomial<PolyRing4, LexOrder> l1;
      l1 += Term<PolyRing4>(+3, {{1, 0, 0, 0}});
      l1 += Term<PolyRing4>(-6, {{0, 1, 0, 0}});
      l1 += Term<PolyRing4>(-2, {{0, 0, 1, 0}});
      l1 += Term<PolyRing4>(0,  {{0, 0, 0, 1}});
     
      Polynomial<PolyRing4, LexOrder> l2;
      l2 += Term<PolyRing4>(+2, {{1, 0, 0, 0}});
      l2 += Term<PolyRing4>(-4, {{0, 1, 0, 0}});
      l2 += Term<PolyRing4>(0,  {{0, 0, 1, 0}});
      l2 += Term<PolyRing4>(+4, {{0, 0, 0, 1}});
    
      Polynomial<PolyRing4, LexOrder> l3;
      l3 += Term<PolyRing4>(+1, {{1, 0, 0, 0}});
      l3 += Term<PolyRing4>(-2, {{0, 1, 0, 0}});
      l3 += Term<PolyRing4>(-1, {{0, 0, 1, 0}});
      l3 += Term<PolyRing4>(-1, {{0, 0, 0, 1}});
 
      // f1_reduced = x - 2y + 2w
      Polynomial<PolyRing4, LexOrder> f1_reduced;
      f1_reduced += Term<PolyRing4>(1,  {{1, 0, 0, 0}});
      f1_reduced += Term<PolyRing4>(-2, {{0, 1, 0, 0}});
      f1_reduced += Term<PolyRing4>(+2, {{0, 0, 0, 1}});

      // f2_reduced = z + 3w
      Polynomial<PolyRing4, LexOrder> f2_reduced;
      f2_reduced += Term<PolyRing4>(1, {{0, 0, 1, 0}});
      f2_reduced += Term<PolyRing4>(3, {{0, 0, 0, 1}});

      auto groebner = runBuchbergers({l1, l2, l3});
      makeMinimalGroebner(groebner); 
      makeReducedGroebner(groebner);

      assert(groebner.size() == 2);
      assert(groebner[0] == f1_reduced);
      assert(groebner[1] == f2_reduced);
   }

   void testBuchbergers3()
   {
      Polynomial<PolyRing2, LexOrder> l1 { {1, {{2,0}}}, {-1, {{0,2}}} };
      Polynomial<PolyRing2, LexOrder> l2 { {1, {{2,0}}}, {1, {{0,1}}} };

      // f1_reduced = x^2+y
      Polynomial<PolyRing2, LexOrder> f1_reduced { {1, {{2,0}}}, {1, {{0,1}}} };
  
      // f2_reduced = y^2+y
      Polynomial<PolyRing2, LexOrder> f2_reduced { {1, {{0,2}}}, {1, {{0,1}}} };

      auto groebner = runBuchbergers({l1, l2});
      makeMinimalGroebner(groebner); 
      makeReducedGroebner(groebner);

      assert(groebner.size() == 2);
      assert(groebner[0] == f1_reduced);
      assert(groebner[1] == f2_reduced);
   }


} // namespace Tests




#endif
