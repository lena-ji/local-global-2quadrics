//////////////////////////////////////////////////////////////////////////////////
// The intersection of quadrics
//////////////////////////////////////////////////////////////////////////////////

k:=Rationals();
RP5<u,v,w,x,y,z>:= PolynomialRing(k, 6);

Q1 := u*v + u*w - 4*v*w + 2*v*z + 2*w*z + x^2 - 2*x*z + y^2 - z^2;
Q2 := u*v - u*w + u*y - 2*v^2 + 2*v*x - 2*w*y + 2*w*z + 2*x*z;
P5:=Proj(RP5);
X:= Scheme(P5, Q1) meet Scheme(P5, Q2);

print "The intersection of quadrics is the", X;

//////////////////////////////////////////////////////////////////////////////////
// Check that X is smooth and contains a rational point
//////////////////////////////////////////////////////////////////////////////////

print "X is smooth?", IsNonsingular(X);

point:=P5![1, 0, 0, 0, 0, 0];
point in X;

print "So X contains the point", point;

//////////////////////////////////////////////////////////////////////////////////
// Check that the associated genus 2 curve C is what we say and that the bad  
// primes are 2 and 149743897
//////////////////////////////////////////////////////////////////////////////////


M1 := SymmetricMatrix(Q1);
M2 := SymmetricMatrix(Q2);

AA<T> := PolynomialRing(Rationals());

MT := M1 - T*M2;
GammaBranch:=-Determinant(MT);

print "Check: This recovers f(T) =", GammaBranch;

Gamma := HyperellipticCurve(GammaBranch);

primes:=Factorization(Integers()!Discriminant(Gamma));

print "The primes dividing the discriminant of Gamma are:", [primes[1][1],primes[2][1]];

//////////////////////////////////////////////////////////////////////////////////
// Check that C has two real Weierstrass points
//////////////////////////////////////////////////////////////////////////////////


print "The real branch points of Gamma are";
Roots(GammaBranch, RealField());

//////////////////////////////////////////////////////////////////////////////////
// Making an affine patch of the Fano variety of lines over F_2
//////////////////////////////////////////////////////////////////////////////////

// We use the affine patch for Gr(2,6) via the chart
//
// u  v  w  x  y  z
// [t1 1 0 t3 t5 t7]
// [t2 0 1 t4 t6 t8]

prime:=2;

k:=GF(prime);
basering<u,v,w,x,y,z>:= PolynomialRing(k, 6);

Q1 := u*v + u*w - 4*v*w + 2*v*z + 2*w*z + x^2 - 2*x*z + y^2 - z^2;
Q2 := u*v - u*w + u*y - 2*v^2 + 2*v*x - 2*w*y + 2*w*z + 2*x*z;

R<[t]> := PolynomialRing(k, 8);
RP1<r, s> := PolynomialRing(R, 2);

   
phi := [t[1] * r + t[2] * s, r, s,  
       t[3] * r + t[4] * s, 
       t[5] * r + t[6] * s, 
       t[7] * r + t[8] * s
       ];
    
eqns1 := Coefficients(Evaluate(Q1, phi));
eqns2 := Coefficients(Evaluate(Q2, phi));
ambient := Spec(R);
Fanopatch := Scheme(ambient, eqns1) meet Scheme(ambient, eqns2);
point1:=ambient![1, 1, 0, 0, 1, 1, 0, 0];
IsNonsingular(Fanopatch, point1);

print "When p =", prime, "the Fano variety of lines contains the smooth point", point1;

//////////////////////////////////////////////////////////////////////////////////
// Making an affine patch of the Fano variety of lines over F_149743897
//////////////////////////////////////////////////////////////////////////////////


prime:=149743897;

k:=GF(prime);
basering<u,v,w,x,y,z>:= PolynomialRing(k, 6);

Q1 := u*v + u*w - 4*v*w + 2*v*z + 2*w*z + x^2 - 2*x*z + y^2 - z^2;
Q2 := u*v - u*w + u*y - 2*v^2 + 2*v*x - 2*w*y + 2*w*z + 2*x*z;

R<[t]> := PolynomialRing(k, 8);
RP1<r, s> := PolynomialRing(R, 2);

   
phi := [t[1] * r + t[2] * s, r, s,  
       t[3] * r + t[4] * s, 
       t[5] * r + t[6] * s, 
       t[7] * r + t[8] * s
       ];
    
eqns1 := Coefficients(Evaluate(Q1, phi));
eqns2 := Coefficients(Evaluate(Q2, phi));
ambient := Spec(R);
Fanopatch := Scheme(ambient, eqns1) meet Scheme(ambient, eqns2);
point2:=ambient![10276, 859210, 113976451, 113430900, 122036333, 94785567, 35411179,
25838500];
IsNonsingular(Fanopatch, point2);

print "When p =", prime, "the Fano variety of lines contains the smooth point", point2;

//////////////////////////////////////////////////////////////////////////////////
// The reduction of X modulo p=149743897 is not conical
//////////////////////////////////////////////////////////////////////////////////

k:=GF(149743897);
R<u,v,w,x,y,z>:= PolynomialRing(k, 6);

Q1 := u*v + u*w - 4*v*w + 2*v*z + 2*w*z + x^2 - 2*x*z + y^2 - z^2;
Q2 := u*v - u*w + u*y - 2*v^2 + 2*v*x - 2*w*y + 2*w*z + 2*x*z;

I := ideal<R | Q1, Q2>;
X := Proj(quo<R | I>);

print "Is the reduction of X modulo 149743897 reduced?", IsReduced(X);
print "Is the reduction of X modulo 149743897 irreducible?", IsIrreducible(X);

print "The dimension of Sing(X_149743897) is", Dimension(SingularSubscheme(X));
print "Is Sing(X_149743897) irreducible?", IsIrreducible(SingularSubscheme(X));
print "Sing(X_149743897) is the", PrimeComponents(SingularSubscheme(X));

P:=SingularSubscheme(X)![10925789 , 85737939 , 85378598 , 93099029 , 51694582 , 1];
eqs := [Q1, Q2];
J := JacobianMatrix(eqs);
print "The Jacobian matrix at P := SingularSubscheme(X_149743897) has rank", Rank(Evaluate(J,[10925789 , 85737939 , 85378598 , 93099029 , 51694582 , 1]));
print "so P is non-conical";

//////////////////////////////////////////////////////////////////////////////////
// The reduction of X modulo p=2 is reducible and non-reduced
//////////////////////////////////////////////////////////////////////////////////

k:=GF(2);
R<u,v,w,x,y,z>:= PolynomialRing(k, 6);

Q1 := u*v + u*w - 4*v*w + 2*v*z + 2*w*z + x^2 - 2*x*z + y^2 - z^2;
Q2 := u*v - u*w + u*y - 2*v^2 + 2*v*x - 2*w*y + 2*w*z + 2*x*z;

I := ideal<R | Q1, Q2>;
X := Proj(quo<R | I>);

print "Is the reduction of X modulo 2 reduced?", IsReduced(X);
print "Is the reduction of X modulo 2 irreducible?", IsIrreducible(X);
