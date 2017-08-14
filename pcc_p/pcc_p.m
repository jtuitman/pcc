////////////////////////////////////////////////// 
// This code is part of the pcc_p MAGMA library //
//                                              //
// copyright (c) 2017 Jan Tuitman               //
//////////////////////////////////////////////////

print "";
print "*****************************";
print "* pcc_p Magma library v2.16 *";
print "*                           *";
print "* by Jan Tuitman            *";
print "*****************************";
print "";

Z:=IntegerRing();
Zx<x>:=PolynomialRing(Z);
Zxy<y>:=PolynomialRing(Zx);

load "auxpolys_p.m";
load "coho_p.m";
load "froblift_p.m";
load "reductions_p.m";
load "zeta_p.m";
print "";

num_zeta:=function(Q,p: N:=0, verbose:=false, exactcoho:=false, W0:=0, Winf:=0);

  if not IsPrime(p) then
    error "p is not prime";
  end if;

  N0:=N;

  ResetMaximumMemoryUsage();
  t0:=Cputime();

  if verbose then
    print "Computing Delta,r,s,W^0,W^inf,G:";
  end if;

  g:=genus(Q,p);
  r,Delta,s:=auxpolys(Q);

  if W0 eq 0 then
    W0:=mat_W0(Q);
  end if;
  if Winf eq 0 then
    Winf:=mat_Winf(Q);
  end if;
  W0inv:=W0^(-1); 
  Winfinv:=Winf^(-1); 

  if (Degree(r) lt 1) or (not smooth(r,p)) or (not (is_integral(W0,p) and is_integral(W0inv,p) and is_integral(Winf,p) and is_integral(Winfinv,p))) then
    error "Bad model for curve";
  end if;

  G:=con_mat(Q,Delta,s);
  G0:=W0*Evaluate(G,Parent(W0[1,1]).1)*W0^(-1)+ddx_mat(W0)*W0^(-1);
  M0:=r*G0;
  Ginf:=Winf*Evaluate(G,Parent(Winf[1,1]).1)*Winf^(-1)+ddx_mat(Winf)*Winf^(-1);
  Jinf,Tinf,Tinfinv:=jordan_inf(Ginf);
  J0,T0,T0inv:=jordan_0(r,G0);
  e0,einf:=ram(J0,Jinf);

  if verbose then
    print "Time (s) :    ", Cputime(t0);
    print "Memory (Mb) : ", GetMaximumMemoryUsage() div (1024^2), "\n";
  end if;

  ResetMaximumMemoryUsage();
  t:=Cputime();

  if N eq 0 then
    prov_prec:=true;
    N:=provable_prec(Q,p,g,W0,Winf,e0,einf);
  else
    prov_prec:=false;
  end if;

  if verbose then
    print "Computing basis H^1(X):"; 
  end if;

  basis,quomap:=basis_coho(Q,p,N,r,W0,Winf,G0,Ginf,J0,Jinf,T0inv,Tinfinv,exactcoho);
  v:=Minimum([0] cat [Valuation(x,p) : x in Eltseq(quomap)]);

  if (v lt 0 and prov_prec) then
    N:=provable_prec(Q,p,g,W0,Winf,e0,einf:val:=v);
  end if;

  if verbose then 
    print "Time (s) :    ", Cputime(t);
    print "Memory (Mb) : ", GetMaximumMemoryUsage() div (1024^2), "\n";
  end if;

  ResetMaximumMemoryUsage();
  t:=Cputime();

  if verbose then
    print "Computing Frobenius lift:";
  end if;

  frobmatb0r:=froblift(Q,p,N-1,r,Delta,s,W0);

  if verbose then 
    print "Time (s) :    ", Cputime(t);
    print "Memory (Mb) : ", GetMaximumMemoryUsage() div (1024^2), "\n";
  end if;

  ResetMaximumMemoryUsage();
  t:=Cputime();

  if verbose then
    print "Computing reduction matrices:";
  end if;

  redlistfin,redlistinf:=red_lists(Q,p,N,r,W0,Winf,G0,Ginf,e0,einf,J0,Jinf,T0,Tinf,T0inv,Tinfinv); 

  if verbose then
    print "Time (s) :    ", Cputime(t);
    print "Memory (Mb) : ", GetMaximumMemoryUsage() div (1024^2), "\n";
  end if;

  ResetMaximumMemoryUsage();
  t:=Cputime();

  if verbose then
    print "Computing Frobenius matrix:";
  end if;

  F:=frobmatrix(Q,p,N,r,W0,Winf,G0,Ginf,frobmatb0r,redlistfin,redlistinf,basis,quomap,verbose);

  if verbose then
    print "";
    print "Time (s) :    ", Cputime(t);
    print "Memory (Mb) : ", GetMaximumMemoryUsage() div (1024^2), "\n";
  end if;

  chi:=revcharpoly(F,p);

  if verbose then
    print "Total time (s):", Cputime(t0), "\n";
    print "The p-adic precision N:", N, "\n";
  end if;


  if not test_num_zeta(chi,p) then
    print "Q =", Q, "p =", p;
    print "Output seems to be wrong. Please report this example to jan.tuitman@kuleuven.be.";
    return 0;
  end if;

  return chi;

end function;
