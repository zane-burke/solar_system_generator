"""
A number of physical and mathematical constants
"""
module Constants
    #=
    Physical constant values taken from NIST CODATA recommended values
    Pi taken from https://oeis.org/A000796 (First 64 Digits)
    =#

    # (π²/60)k⁴/ħ³c² 
    const StefanBoltzmannConstant::BigFloat = BigFloat(5.670374419e-8) # W m⁻² K⁻⁴
    export StefanBoltzmannConstant
    const BoltzmannConstant::BigFloat = BigFloat(1.380649e-23) # J K⁻¹
    export BoltzmannConstant
    const GravitationalConstant::BigFloat = BigFloat(6.6743e-11) # m³ kg⁻¹ s⁻²
    export GravitationalConstant
    const PlanckConstant::BigFloat = BigFloat(6.62607015e-34) # J Hz⁻¹
    export PlanckConstant
    # h/2π
    const ReducedPlanckConstant::BigFloat = BigFloat(1.054571817e-34) # J s
    export ReducedPlanckConstant
    const WienDisplacementConstant::BigFloat = BigFloat(2.897771955e-3) # m K
    export WienDisplacementConstant
    const SpeedOfLight::BigFloat = BigFloat(299792458) # m s⁻¹
    export SpeedOfLight
    const ZeroPointLuminosity::BigFloat = BigFloat(3.0128e28) # W
    export ZeroPointLuminosity
    const SolarMass = 1.9885e30 # kg
    export SolarMass
    const SolarRadius = 695700000 # m
    export SolarRadius
end