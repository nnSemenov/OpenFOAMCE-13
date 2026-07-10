= Modeling capillary force in `multiphaseEuler`


In most papers, different fluid phases have their own pressure in momentum equation:
$
    (partial (alpha_i rho_i U_i)) / (partial t) + nabla dot (alpha_i rho_i U_i U_i) = nabla dot tau_i - alpha_i nabla p_i + alpha_i rho_i g...
$

== Capillary pressure model

Capillary pressure is defined to be pressure difference between phases. They are given by specific physical models.

The most common capillary model is Leverett J function, closing J as function of liquid saturation $S_L$:
$
    (P_c)^i_j = p_i - p_j = sqrt(epsilon/K) sigma cos theta J(S_L)
$

$sqrt(K/epsilon)$ is some kind of hydrolic pore diameter for rocks. For packed spheres,
$
    sqrt(epsilon/K) = (1-epsilon) / (epsilon d_p) sqrt(E_1)
$
Where $epsilon$ is porosity, $K$ is permeability, $d_p$ is particle diameter, $E_1$ is first Ergun constant.

== Capillary system
Capillary system organise pressure models in star-like. A phase is assigned to be reference phase (only in this system, no relation to $alpha$ solving), and all other fluid phases need to have a capillary model with it. The term `phaseA_to_phaseB` refers to capillary pressure $p_A - p_B$.

Like most Euler-Euler multiphase solvers, `multiphaseEuler` assumes all phases share common pressure $p$. This is not breakable otherwise we rebuild everything. Here we use an convention to link common pressure with each fluid phase's own pressure:
$
    p = p_"ref"
$

(Volume average is not practicable, because leads to tons of $nabla alpha_i$)

Then we have
$
    -alpha_i nabla p + F_"cap" = -alpha_i nabla p_i \
    F_"cap" = alpha_i nabla(p-p_i) \
    p-p_i = p_"ref" - p_i = (P_c)^"ref"_i \
$

== Relation with `p_rgh`

When with gravity, $p_"rgh"$ is solved instead of $p$. It's defined that
$
    p_"rgh" = p - rho bold(g) dot bold(h) \
    nabla p_"rgh" = nabla p - rho bold(g) - (bold(g) dot bold(h))nabla rho \
$

There may be some details around $(bold(g) dot bold(h))nabla rho$ but that's nothing to do with capillary force. It's still okay that
$
    F_"cap" = alpha_i nabla(p-p_i)
$
/*
== Semi-implicit discretization
For any moving fluid phase $i$ (not reference phase), momentum balance is approximatly good:
$
(D)/(D t)(alpha_i rho_i U_i)=0 \
$
$
0 &= -alpha_i nabla p + alpha_i rho_i g + K_("ref",i)(U_"ref"-U_i)+ sum_(j!=i\ j!="ref")K_(i j)(U_j-U_i) + (dif P_(c,i))/(dif alpha_i) nabla alpha_i \
0 &= -alpha_"ref" nabla p + alpha_"ref" rho_"ref" g + K_("ref",i)(U_i - U_"ref") + sum_(j!=i\ j!="ref")K_(i j)(U_j-U_"ref")
$

Where $j$ refers to all phases including stationary solid. For stationary phases, $U_j=0$. If without drag force, $K_(i j)=0$.

By elimiating pressure gradient:
$
U_i [K_("ref",i)(alpha_i+alpha_"ref")/(alpha_i alpha_"ref")+sum_(j!=i\ j!="ref")K_(i j)/alpha_i] =
    & U_"ref" [(alpha_i+alpha_"ref")/(alpha_i alpha_"ref")K_("ref",i) + sum_(j!=i\ j!="ref")K_("ref",j)/alpha_"ref"] + sum_(j!=i\ j!="ref") + (rho_i - rho_"ref")g + (dif P_(c,i))/(dif alpha_i) nabla alpha_i \
    & + sum_(j!=i\ j!="ref")(K_(i j)/alpha_i + K_("ref",j)/alpha_"ref") U_j
$
Multiplying with $alpha_i$ on both sides:
$
U_i underbrace([K_("ref",i)(1+alpha_i/alpha_"ref") + sum_(j!=i\ j!="ref") K_(i j)],"Deno") =
    & U_"ref" [(1+alpha_i/alpha_"ref")K_("ref",i) + sum_(j!=i\ j!="ref")alpha_i/alpha_"ref"K_("ref",j)] + alpha_i (rho_i-rho_"ref")g \
    & + (dif P_(c,i))/(dif ln alpha_i) nabla alpha_i +  sum_(j!=i\ j!="ref")(K_(i j) + alpha_i/alpha_"ref" K_("ref",j) ) U_j
$

Then, discrete phase flux is expressed as other velocities:
$
alpha_i U_i & = alpha_i ((1+alpha_i/alpha_"ref")K_("ref",i) + sum_(j!=i\ j!="ref")alpha_i/alpha_"ref"K_("ref",j))/("Deno") U_"ref" + (alpha_i^2(rho_i-rho_"ref"))/"Deno" g + (alpha_i)/"Deno" (dif P_(c,i))/(dif ln alpha_i) nabla alpha_i \
            & + alpha_i (sum_(j!=i\ j!="ref")(K_(i j) + alpha_i/alpha_"ref" K_("ref",j))U_j)/"Deno"
$

Taking divergence of each terms:
1. Terms with vectors can be made into flux(`fvc::flux`) and implicitly discretized (`fvm::div`).
2. Gravity term can be discretized with `fvm::SuSp`:
$
nabla dot ((alpha_i^2 (rho_i - rho_"ref"))/"Deno" g) = alpha_i [nabla dot ((alpha_i (rho_i - rho_"ref"))/"Deno" g) + (nabla alpha_i) dot ((rho_i-rho_"ref")/"Deno" g)]
$
3. Capillary term is discretized with `fvm::laplacian`:
$
nabla dot (alpha_i/"Deno" (dif P_(c,i))/(dif ln alpha_i) nabla alpha_i)
$

*/