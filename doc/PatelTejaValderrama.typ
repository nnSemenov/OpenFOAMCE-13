= Residual properties for PTV EOS
$
p = (R T)/(V-b) - a/(V(V+b)+c(V-b)) \

V^3 + (c-(R T)/p)V^2 + (a/p - (b+c)(R T)/p - b^2 -2 b c)V + b^2c + (b c R T)/p -(a b)/p = 0 \

$

where
$
a = a_c alpha(T_r) \
alpha(T_r) = [1+F(1-sqrt(T_r))]^2 \
a_c = Omega_a (R T_c)^2 / p_c \
b=Omega_b (R T_c)/p_c \
c=Omega_c (R T_c)/p_c \
Omega_a = 0.66121-0.761057Z_c \
Omega_b = 0.02207+0.20868Z_c \
Omega_c = 0.57765-1.87080Z_c \
F = 0.46283 + 3.58230 omega Z_c + 8.194 17(omega Z_c^2)^2
$

== Derivatives
$
((partial p)/(partial T))_V &= R/(V-b) - a'/(V(V+b)+c(V-b)) \

((partial^2 p)/(partial T^2))_V &= -(a'')/(V(V+b)+c(V-b)) \

((partial p)/(partial V))_T &= - (R T)/(V-b)^2 + (a (2V+b+c))/(V^2 +(b+c)V - b c)^2 \


((partial^2 V)/(partial T^2))_p &= ((2R V +R c+R b-a')((partial V)/(partial T))_p -a''(V-b))/(3p V^2 + 2(p c - R T)V+a-(b+c)R T-(b+2c)b p) \
    & - ((R V^2+(R b+R c-a')V+a'b-R b c)[(6p V+2p c - 2R T)((partial V)/(partial T))_p - 2 R V+a'-R(b+c)])/(3p V^2+2(p c - R T)V+a-(b+c)R T-(b+2c)b p)^2 \

(partial^2p)/(partial V partial T) &= -R/(V-b)^2 + a'(2V+b+c)/(V^2 +(b+c)V - b c)^2
$

Derivatives of $a$:
$
a &= a_c alpha(T_r) \

a'
 &= a_c (dif alpha(T_r))/(dif T_r) dot 1/T_c \
 &= a_c/T_c F (F(sqrt(T_r)-1)-1)/sqrt(T_r) \

a''
 &= a_c/T_c^2 (F(F+1))/(2T_r sqrt(T_r))
$

== Residual internal energy and enthalpy
$
U_R/(R T)
 &= 1/(R T) integral_infinity^V [T((partial p)/(partial T))_V-p] dif V \
 &= 1/(R T) integral_infinity^V (a'T - a)/(V(V+b)+c(V-b)) dif V \
 &= (a'T-a)/(R T sqrt(c^2+6b c+b^2)) ln((2V-sqrt(c^2+6b c+b^2)+c+b)/(2V+sqrt(c^2+6b c+b^2)+c+b)) \
$

$
H^R/(R T) = Z - 1 + U^R/(R T)
$

== Residual specific heat
$
C_V^R/R &= T/R integral_infinity^V ((partial^2 p)/(partial T^2))_V dif V \
    &= - (a''T)/(R sqrt(c^2+6 b c+b^2))ln((2V-sqrt(c^2+6b c+b^2)+c+b)/(2V+sqrt(c^2+6b c+b^2)+c+b)) \
$

$
(C_p^R - C_V^R)/R = -1 + T/R ((partial p)/(partial T))_V ((partial V)/(partial T))_p \
$

== Derivative of residual specific heat
$
((partial C_v^R)/(partial T))_p
  &= (partial /(partial T))_p [T integral_infinity^V ((partial^2p)/(partial T^2))_V dif V] \
  &=  integral_infinity^V ((partial^2p)/(partial T^2))_V dif V + T ((partial V)/(partial T))_p ((partial^2p)/(partial T^2))_V \
  &=- (a'')/sqrt(c^2+6 b c+b^2) ln((2V-sqrt(c^2+6b c+b^2)+c+b)/(2V+sqrt(c^2+6b c+b^2)+c+b)) + T ((partial V)/(partial T))_p ((partial^2p)/(partial T^2))_V\
$

$
((partial (C_p^R-C_V^R))/(partial T))_p
 &= (partial/(partial T))_p dot [T ((partial p)/(partial T))_V ((partial V)/(partial T))_p] \
 &= ((partial V)/(partial T))_p dot [((partial p)/(partial T))_V + T (partial/(partial T))_p dot ((partial p)/(partial T))_V] + T ((partial p)/(partial T))_V ((partial^2 V)/(partial T^2))_p \

(partial/(partial T))_p dot ((partial p)/(partial T))_V
 &= ((partial^2p)/(partial T^2))_V + (partial^2p)/(partial V partial T) ((partial V)/(partial T))_p
$