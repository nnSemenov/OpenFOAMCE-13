# About
Mikeno is a fork of OpenFOAM, it's Frankensteined for chemical engineering usage.

[简体中文](README.md) [English]

## Goal of this fork

1. Better support for modelling chemical engineering process
2. Continually merge upstream updates
3. Build system migration: Replace wmake with CMake
4. Fix unexpected SIGFPE trapping, mainly for single-precision.
   - Compiler optimises float-point computation with SIMD, but they generate unexpected NAN sometimes. These NANs are
     never used, but emit SIGFPE.

## Modifications

### New Utilities
1. `autoCompress`: Compress all heavy fields (non-uniform, including mesh) to gzip, regardless of format(ascii or binary)

### Code and compiler
1. Migrated all subprojects from wmake to modern cmake.
   - Compatible to modern IDEs, high-quality linting from clangd.
   - More convenient to import external dependencies.
   - More efficient parallel build.
   - More convenient for secondary development with `find_package(Mikeno CONFIG)`.
2. Performance optimisation
   - Add `-march=native` for optimisation mode, enabling more vectorisation from compiler.
3. Support `AOCC`
4. Remove redundant reference in solver modules.
   1. Currently cleaned `isothermalFluidSolver` `fluidSolver` `multicomponentFluidSolver` `XiFluidSolver`

### Support arbitrarily high pressure
1. The `pOffset` keyword can be added to `physicalProperties`, allowing gauge pressure for p-V coupling, and absolute pressure for thermophysical. **Small pressure difference will never be flooded by floating point rounding error(especially single precision)**

### Porous media heat transfer
1. `porousMediaFluidSolver` supporting arbitrary number of porous phases, heat transfer between fluid-porous and porous-porous, both supporting thermal equilibrium and non-equilibrium.

### More rigorous thermodynamics
1. Equation of state: add `RedlichKwongGas`, rewrite `PengRobinsonGas`
2. Both real gas EOSs are tested with AspenPlusV14
3. Van der Waals mixing rule
   - Support binary interaction coefficient `k_ij` for `PengRobinsonGas`
4. All thermo models (eq `hConst`) supports `ideal_*` methods that gives ideal gas properties.
5. New mixture model `realGasMulticomponentMixture`
   1. Compute `rho` from mixed real gas EOS
   2. Compute and add residual properties of `Cp` `Cv` `hs` `ha` `es` `ea`. All residual properties are computed from mixed EOS.
6. New mixture model `idealLiquidMulticomponentMixture`
   1. Apply volume weighed mixing for density (no extra volume)
   2. Apply mass weighed mixing for he and Cp Cv (no extra energy)
   3. Arrhenius mixing for mu
   4. Vredeveled mixing for kappa ($n=-2$)
7. Add `phase` keyword to `equationOfState` dictionary, allow using real fluid EOS for liquid.
8. Equation of states implemenets `dCpdT` `dCvdT` (residual specific heat derivative to temperature) for chemical solver usage.
   - Currently only `dCvdT` is rigorous. Derivative of difference of residual specific heat is currently ignored (TOO COMPLEX for cubic EOS)

## Fix unexpected SIGFPE trapping (compile option in brackets)
1. Fix `flowRateInletVelocity` trapped by SIGFPE when writing flow field. This is caused by division in `unitConversion::toUser(const T& t) const`. (`Clang DP Opt`)
2. Fix `chemistryModel` trapped by SIGFPE when calculating reaction rate. This is is trapped in `Foam::Reaction<ThermoType>::C`, probably caused by branching logic. (`Clang DP Release`)
3. Fix unexpected SIGFPE globally by adding `-ffp-exception-behavior=maytrap`. (`Clang`) This won't slow down because CFD performance is restricted to memory access but not float computation.

## Bug fixes
1. Fix conflicts between `AndradeTransport` and chemistry. Fill missing implementation of `operator+` and `operator*`.
2. Fix segfault when `foamPostProcess` runs function object with cellZone.
   - This was caused by `3caf09d88b95092a1f4a6047cf498d47d5e86a7a` which aims to optimize performance.
   - This optimization is reverted, currently no perfect way to fix.

## Pending works
1. More equation of state: Patel-Teja, Martin-Hou
2. Extend porous media heat transfer to multicomponent
3. Stabilize `porousMediaFluid` for non-equilibrium heat transfer with large coefficient or specific area

## Existing Bugs(Up to 2026-04-01):
1. `foamRun` crashes with Lagrangian fields (`test/Lagrangian/boundaries` fails)