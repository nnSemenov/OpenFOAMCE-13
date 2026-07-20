inline scalar ideal_Cv
(
    const scalar p,
    const scalar T
) const
{
    return ideal_Cp(p, T) - this->R();
}

inline scalar Cv
(
    const scalar p,
    const scalar T
) const
{
    return Cp(p, T) - EquationOfState::CpMCv(p, T);
}


// Ideal Sensible internal energy [J/kg]
inline scalar ideal_es(const scalar p, const scalar T) const {
    return ideal_hs(p, T) - p/EquationOfState::rho(p, T);
}

// Ideal Absolute internal energy [J/kg]
inline scalar ideal_ea(const scalar p, const scalar T) const {
    return ideal_ha(p, T) - p/EquationOfState::rho(p, T);
}

inline scalar es
(
    const scalar p,
    const scalar T
) const
{
    return hs(p, T) - p/EquationOfState::rho(p, T);
}

inline scalar ea
(
    const scalar p,
    const scalar T
) const
{
    return ha(p, T) - p/EquationOfState::rho(p, T);
}

inline scalar dedp_T
(
    const scalar p,
    const scalar T
) const {
    return EquationOfState::dedp_T(p, T);
}

inline scalar dhdp_T
(
    const scalar p,
    const scalar T
) const {
    return EquationOfState::dhdp_T(p, T);
}

static_assert(not std::is_same_v<EquationOfState, Foam::specie>,
              "thermo class shouldn't inherit directly from Foam::specie, "
              "don't skip EOS.");