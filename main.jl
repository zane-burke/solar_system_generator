include("astrovector.jl")
using .Astrovector
include("constants.jl")
using .Constants
include("coordvectors.jl")
using .CoordVector
include("utils.jl")
using .Utils
include("solar.jl")
using .Solar
include("generators.jl")
using .Generators
#= 
    Known Limitations:

    - Only works for Main-Sequence stars
    - Cannot do binary star systems
    - Does not support 3d space

    - Planetary generation is very limited
      - So-called hot Jupiters cannot exist, since the inner planets are always rocky planets
=#

#=

Steps:

    Generate the Star:
        Generate the class of the star
        Generate the temperature of the star
        Generate the mass of the star
        Generate the radius of the star

    Generate the number of planets


    Generate the first planet


    Generate the second planet

    ...

    Generate the last planet

=#

#=
number_of_planets_weights = [0.0438, 0.1350, 0.3238, 0.6049, 0.8802, 0.9974, 0.8802, 0.6049, 0.3238, 0.1350, 0.438]
number_of_planets = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
=#

function main()
    #=
    Predefine:
        Orbital:
            Semi-major Axis
            Eccentricity

        Physical:
            Mass
            Semi-major Axis
    =#
    star = generatestar()
    printstar(star)
end
main()

