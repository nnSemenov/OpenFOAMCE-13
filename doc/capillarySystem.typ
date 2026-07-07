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

== Discretization
1. Very hard to be implicit
Although $F_"cap"$ is in propotional with $nabla alpha_L$, it's not very easy to process like turbulence dispersion force. Capillary has no relation with Reynolds' average so there shouldn't be a diffusion term in alpha. This can be done by replacing $U_"dispered"$ with $U_"continuous"$, but too specilizated for 3-phase system and extremely conflicts `multiphaseEuler`. Capillary force is finally discretized explicitly.

2. Explicit on surface
In `multiphaseEuler`, lift, phase pressure (solid's fake pressure) and turbulent dispersion force is added explicitly to momentum equation. But not as a `volVectorField`, but `surfaceScalarField` for stabilization. Pure cell-center often lead to chess-plate ocsillation.
$
F_f = arrow(F) dot arrow(S)_f  \
$
Although `multiphaseEuler` have both cell and face based pressure corrector, forces are all turned into fluxes on surface(`fvc::flux(F)`), just minor difference in order of interpolation and flux. For capillary force, both cell and face are equal:
$
F_f = alpha_i underbrace(nabla (p-p_i) dot arrow(S)_f,"snGrad" dot "magSf")
$
