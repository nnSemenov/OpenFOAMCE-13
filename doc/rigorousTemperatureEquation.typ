= Rigorous Temperature Equation Derived from Energy Equation

== From Specific Enthalpy

Enthalpy equation:
$
  (partial (rho h))/(partial t) + nabla dot (rho bold(U) h) + (partial (rho K))/(partial t) + nabla dot (rho bold(U) K) - (partial p)/(partial t) + nabla dot bold(q) = S_h
$

Full derivative of enthalpy:
$
  dif h = underbrace(((partial h)/(partial p))_p, C_p) dif T + underbrace(((partial h)/(partial p))_T, Ж) dif p
$

Represent $h$ with $T$:
$
  (partial (rho h))/(partial t) + nabla dot (rho bold(U) h)
  & = underbrace(h (partial rho)/(partial t) + h nabla dot (rho bold(U)), =0) + rho (partial h)/(partial t) rho bold(U) dot nabla h \
  & = rho C_p (partial T)/(partial t) + rho Ж (partial p)/(partial t) + rho bold(U) dot C_p nabla T + Ж rho bold(U) dot nabla p \
  & = rho C_p (partial T)/(partial t) + underbrace(nabla dot (rho bold(U) C_p T), "fvm::div") - underbrace(T nabla dot (rho bold(U) C_p), "fvm::SuSp") + Ж(rho (partial p)/(partial t) + rho bold(U) dot nabla p)
$
Diffusion terms becomes:
$
  nabla dot bold(q) & = - nabla dot (kappa nabla T)
$

The whole energy equation becomes:
$
  rho C_p underbrace((partial T)/(partial t), "fvm::ddt") + underbrace(nabla dot (rho bold(U) C_p T), "fvm::div") - underbrace(T nabla dot (rho bold(U) C_p), "fvm::SuSp") + Ж(rho (partial p)/(partial t) + rho bold(U) dot nabla p) + (partial (rho K))/(partial t) + nabla dot (rho bold(U) K) - (partial p)/(partial t) - underbrace(nabla dot (kappa nabla T), "fvm::laplacian") = S_h
$
Other terms should be discretized explicitly.

Isothermal compression enthalpy factor should be calculated from thermodynamics:
$
  ((partial H)/(partial p))_T & = T((partial S)/(partial p))_T + V                       \
                              & = V - T((partial V)/(partial T))_p                       \
                            h & := H/M_w                                                 \
                            Ж & = 1/M_w (V - T((partial V)/(partial T))_p )              \
                              & = 1/rho + limits(R)^~ (z + T((partial z)/(partial T))_p)
$
where $H$ is in J/mol, $h$ is in J/kg.

For ideal gases, $Ж=0$.

For Redlich-Kwong gases, $Ж$ is given as
$
  ((partial V)/(partial T))_p
  = ((R V^2)/p + (a/(2p sqrt(T))+(R b)/p )V + (a b)/(2 p T sqrt(T)) )
  / ( 3V^2 - 2(R T)/p V + (a/(p sqrt(T)) - b^2 - (b R T)/p) )
$

== From Internal Energy

Similar to enthalpy, isothermal energy to pressure derivative is given as
$
  ((partial U)/(partial p))_T = - T ((partial V)/(partial T))_p - p ((partial V)/(partial p))_T \
  ((partial V)/(partial p))_T =
  ( -(R T )/p^2 V^2 -((b R T)/p - a/(p^2 sqrt(T)))V + (a b)/(p^2 sqrt(T)) )
  / ( 3V^2 -2 (R T)/p V + a/(p sqrt(T)) - (b R T)/p - b^2 )
$

== PengRobinson
Similar to Redlich-Kwong EOS, derivative can be calculated from its polynomial form
$
  V^3 + V^2(b-(R T)/p) + V dot (a/p -2 (b R T)/p -3b^2) + b dot ((b R T)/p +b^2 -a/p) =0 \
  ((partial V)/(partial T))_p = 
  (R V^2 + (2 R b - a') V + b(a'-R b))
  / (
  3 p V^2 + 2V (p b- R T) + a - 2 b R T - 3 p b^2
  ) \
  ((partial V)/(partial p))_T =
  (-R T V^2 - (2 R T b - a)V - a b + b^2 R T) /
  (3 p^2V^2 + 2 p V (p b- R T) + a p - 2 b R T p - 3 p^2 b^2)
$
where $a'=(dif a)/(dif T)$