//////////////////////////////////////////////////////////////////////////////////
// The intersection of quadrics
//////////////////////////////////////////////////////////////////////////////////

k:=Rationals();
RP5<u,v,w,x,y,z>:= PolynomialRing(k, 6);

Q1 := 1471*u^2 + 13372*u*v + 12461*u*w - 27530*u*x + 22448*u*y + 3650*u*z + 
    30388*v^2 + 56632*v*w - 125128*v*x + 102064*v*y + 16580*v*z + 26386*w^2 -
    116594*w*x + 95092*w*y + 15451*w*z + 128811*x^2 - 210156*x*y - 34132*x*z +
    85760*y^2 + 27832*y*z + 2261*z^2;
Q2 := -1294*u^2 - 11768*u*v - 10966*u*w + 24229*u*x - 19756*u*y - 3214*u*z -
    26744*v^2 - 49844*v*w + 110120*v*x - 89792*v*y - 14600*v*z - 23225*w^2 +
    102616*w*x - 83656*w*y - 13610*w*z - 113357*x^2 + 184884*x*y + 30051*x*z -
    75468*y^2 - 24464*y*z - 1996*z^2;
P5:=Proj(RP5);
X:= Scheme(P5, Q1) meet Scheme(P5, Q2);

print "The intersection of quadrics is the", X;


//////////////////////////////////////////////////////////////////////////////////
// Check that X is smooth and contains a rational point
//////////////////////////////////////////////////////////////////////////////////

print "X is smooth?", IsNonsingular(X);

point:=P5![2, -5/7, -4/7, -1/7, 1/7, 1];
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

// We use the first standard affine patch for Gr(2,6):
//
// u  v  |  w  x  y   z
// [1 0. | t1 t3 t5 t7 ]
// [0 1  | t2 t4 t6 t8 ]

prime:=2;

k:=GF(prime);
basering<u,v,w,x,y,z>:= PolynomialRing(k, 6);

Q1 := 1471*u^2 + 13372*u*v + 12461*u*w - 27530*u*x + 22448*u*y + 3650*u*z + 
    30388*v^2 + 56632*v*w - 125128*v*x + 102064*v*y + 16580*v*z + 26386*w^2 -
    116594*w*x + 95092*w*y + 15451*w*z + 128811*x^2 - 210156*x*y - 34132*x*z +
    85760*y^2 + 27832*y*z + 2261*z^2;
Q2 := -1294*u^2 - 11768*u*v - 10966*u*w + 24229*u*x - 19756*u*y - 3214*u*z -
    26744*v^2 - 49844*v*w + 110120*v*x - 89792*v*y - 14600*v*z - 23225*w^2 +
    102616*w*x - 83656*w*y - 13610*w*z - 113357*x^2 + 184884*x*y + 30051*x*z -
    75468*y^2 - 24464*y*z - 1996*z^2;

R<[t]> := PolynomialRing(k, 8);
RP1<X, Y> := PolynomialRing(R, 2);

   
phi := [X, Y, (t[1] * X + t[2] * Y ),
       (t[3] * X + t[4] * Y ),
        (t[5] * X + t[6] * Y ),
       (t[7] * X + t[8] * Y )
       ];
    
eqns1 := Coefficients(Evaluate(Q1, phi));
eqns2 := Coefficients(Evaluate(Q2, phi));
ambient := Spec(R);
Fanopatch := Scheme(ambient, eqns1) meet Scheme(ambient, eqns2);
point1:=ambient![0, 0, 0, 1, 0, 0, 1, 1];
IsNonsingular(Fanopatch, point1);

print "When p =", prime, "the Fano variety of lines contains the smooth point", point1;


//////////////////////////////////////////////////////////////////////////////////
// Making an affine patch of the Fano variety of lines over F_149743897
//////////////////////////////////////////////////////////////////////////////////


prime:=149743897;

k:=GF(prime);
basering<u,v,w,x,y,z>:= PolynomialRing(k, 6);

Q1 := 1471*u^2 + 13372*u*v + 12461*u*w - 27530*u*x + 22448*u*y + 3650*u*z + 
    30388*v^2 + 56632*v*w - 125128*v*x + 102064*v*y + 16580*v*z + 26386*w^2 -
    116594*w*x + 95092*w*y + 15451*w*z + 128811*x^2 - 210156*x*y - 34132*x*z +
    85760*y^2 + 27832*y*z + 2261*z^2;
Q2 := -1294*u^2 - 11768*u*v - 10966*u*w + 24229*u*x - 19756*u*y - 3214*u*z -
    26744*v^2 - 49844*v*w + 110120*v*x - 89792*v*y - 14600*v*z - 23225*w^2 +
    102616*w*x - 83656*w*y - 13610*w*z - 113357*x^2 + 184884*x*y + 30051*x*z -
    75468*y^2 - 24464*y*z - 1996*z^2;

R<[t]> := PolynomialRing(k, 8);
RP1<X, Y> := PolynomialRing(R, 2);

   
phi := [X, Y, (t[1] * X + t[2] * Y ),
       (t[3] * X + t[4] * Y ),
        (t[5] * X + t[6] * Y ),
       (t[7] * X + t[8] * Y )
       ];
    
eqns1 := Coefficients(Evaluate(Q1, phi));
eqns2 := Coefficients(Evaluate(Q2, phi));
ambient := Spec(R);
Fanopatch := Scheme(ambient, eqns1) meet Scheme(ambient, eqns2);
point2:=ambient![23367516, 38034557, 52902585, 9110600, 36955379, 15930411,78325746,139150875];
IsNonsingular(Fanopatch, point2);

print "When p =", prime, "the Fano variety of lines contains the smooth point", point2;

//////////////////////////////////////////////////////////////////////////////////
// The reduction of X modulo p=149743897
//////////////////////////////////////////////////////////////////////////////////

k:=GF(149743897);
R<u,v,w,x,y,z>:= PolynomialRing(k, 6);

Q1 := 1471*u^2 + 13372*u*v + 12461*u*w - 27530*u*x + 22448*u*y + 3650*u*z + 
    30388*v^2 + 56632*v*w - 125128*v*x + 102064*v*y + 16580*v*z + 26386*w^2 -
    116594*w*x + 95092*w*y + 15451*w*z + 128811*x^2 - 210156*x*y - 34132*x*z +
    85760*y^2 + 27832*y*z + 2261*z^2;
Q2 := -1294*u^2 - 11768*u*v - 10966*u*w + 24229*u*x - 19756*u*y - 3214*u*z -
    26744*v^2 - 49844*v*w + 110120*v*x - 89792*v*y - 14600*v*z - 23225*w^2 +
    102616*w*x - 83656*w*y - 13610*w*z - 113357*x^2 + 184884*x*y + 30051*x*z -
    75468*y^2 - 24464*y*z - 1996*z^2;

I := ideal<R | Q1, Q2>;
X := Proj(quo<R | I>);

print "Is the reduction of X modulo 149743897 reduced?", IsReduced(X);
print "Is the reduction of X modulo 149743897 irreducible?", IsIrreducible(X);

print "The dimension of Sing(X_149743897) is", Dimension(SingularSubscheme(X));
print "Is Sing(X_149743897) irreducible?", IsIrreducible(SingularSubscheme(X));
print "Sing(X_149743897) is the", PrimeComponents(SingularSubscheme(X));

P:=SingularSubscheme(X)![87326150 , 44639405 , 42946036 , 20149960 , 14475939 , 1];
eqs := [Q1, Q2];
J := JacobianMatrix(eqs);
print "The Jacobian matrix at P := SingularSubscheme(X_149743897) has rank", Rank(Evaluate(J,[87326150 , 44639405 , 42946036 , 20149960 , 14475939 , 1]));
print "so P is non-conical";

//////////////////////////////////////////////////////////////////////////////////
// The reduction of X modulo p=2
//////////////////////////////////////////////////////////////////////////////////


k:=GF(2);
R<u,v,w,x,y,z>:= PolynomialRing(k, 6);

Q1 := 1471*u^2 + 13372*u*v + 12461*u*w - 27530*u*x + 22448*u*y + 3650*u*z + 
    30388*v^2 + 56632*v*w - 125128*v*x + 102064*v*y + 16580*v*z + 26386*w^2 -
    116594*w*x + 95092*w*y + 15451*w*z + 128811*x^2 - 210156*x*y - 34132*x*z +
    85760*y^2 + 27832*y*z + 2261*z^2;
Q2 := -1294*u^2 - 11768*u*v - 10966*u*w + 24229*u*x - 19756*u*y - 3214*u*z -
    26744*v^2 - 49844*v*w + 110120*v*x - 89792*v*y - 14600*v*z - 23225*w^2 +
    102616*w*x - 83656*w*y - 13610*w*z - 113357*x^2 + 184884*x*y + 30051*x*z -
    75468*y^2 - 24464*y*z - 1996*z^2;

I := ideal<R | Q1, Q2>;
X := Proj(quo<R | I>);

print "Is the reduction of X modulo 2 reduced?", IsReduced(X);
print "Is the reduction of X modulo 2 irreducible?", IsIrreducible(X);

print "The defining equations of X modulo 2 after the coordinate change u -> u + z are";
1471*(u+z)^2 + 13372*(u+z)*v + 12461*(u+z)*w - 27530*(u+z)*x + 22448*(u+z)*y + 3650*(u+z)*z + 
    30388*v^2 + 56632*v*w - 125128*v*x + 102064*v*y + 16580*v*z + 26386*w^2 -
    116594*w*x + 95092*w*y + 15451*w*z + 128811*x^2 - 210156*x*y - 34132*x*z +
    85760*y^2 + 27832*y*z + 2261*z^2;
-1294*(u+z)^2 - 11768*(u+z)*v - 10966*(u+z)*w + 24229*(u+z)*x - 19756*(u+z)*y - 3214*(u+z)*z -
    26744*v^2 - 49844*v*w + 110120*v*x - 89792*v*y - 14600*v*z - 23225*w^2 +
    102616*w*x - 83656*w*y - 13610*w*z - 113357*x^2 + 184884*x*y + 30051*x*z -
    75468*y^2 - 24464*y*z - 1996*z^2;