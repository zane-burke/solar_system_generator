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
    - The number of planets maxes out at 10, even though a star could theoretically have many more planets than that

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

