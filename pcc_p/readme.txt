////////////////////////////////////////////////// 
// This code is part of the pcc_p MAGMA library //
//                                              //
// copyright (c) 2017 Jan Tuitman               //
//////////////////////////////////////////////////

This is the pcc_p (Point Counting on Curves) MAGMA library v2.16 for computing zeta functions of curves over finite (prime) fields.

The code is based on the algorithm presented in my papers:

"Counting points on curves using a map to P^1", http://arxiv.org/abs/1402.6758
"Counting points on curves using a map to P^1", II, http://arxiv.org/abs/1412.7217

Jan Tuitman
jan.tuitman@wis.kuleuven.be

February 2017

GETTING STARTED
---------------
To use the code:

1) Put all the files into one and the same directory.

2) Open the file "example_p.m" in an editor and uncomment one of the lines specifying 

  Q: a polynomial in Z[x][y], monic in y
  p: a prime number

  or of course add a new line of your own (Q has to be monic in y)!

3) Load "pcc_p.m" and "example_p.m" into MAGMA (in this order).

4) Call chi:=num_zeta(Q,p); 

The code will output the numerator chi of the zeta function of the smooth projective model of the curve defined by Q over the finite field with p elements. 

OPTIONAL PARAMETERS
-------------------
N        : p-adic precision for the Frobenius matrix (by default the code will compute a provably correct precision)

verbose  : prints the time and memory spent on the various steps to the screen

exactcoho: carry out the linear algebra for computing a basis for the cohomology to exact precision, which is a lot slower

An example using all of these options: 

chi:=num_zeta(Q,p: N:=20,verbose:=true,exactcoho:=true); 






