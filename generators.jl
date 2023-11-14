module Generators
    using Random
    using StatsBase

    include("utils.jl")
    using .Utils
    include("solar.jl")
    using .Solar
    include("constants.jl")
    using .Constants

    """
    Randomly generates a stellar classification (i.e., G2, K5, B1, etc)

    Returned as follows: Letter, Number

    Currently only generates Main-Sequence stars. 
    """
    function _stellarclass()
        # O=>weights[1], B=>weights[2], ... for A, F, G, K, M
        # Weights based on the true rarity of that type of star
        class_weights = [0.00003, 0.12, 0.61, 3.0, 7.6, 12, 76]
        class_items = ["O", "B", "A", "F", "G", "K", "M"]
        # Heavily weight the subclasses towards lower-subclass, while still allowing for higher subclass stars
        subclass_weights = [0.2656, 0.3438, 0.4311, 0.5299, 0.6438, 0.7781, 0.9418, 1.1514, 1.4435, 1.9285]
        subclass_items = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
        return sample(class_items, Weights(class_weights)), sample(subclass_items, Weights(subclass_weights))
    end

    function _stellar_rotation_velocity(mass::BigFloat)
        mass_constant = 10 ^ 30
        # Based on a logarithmic regression of the Mass vs. Rotation Velocity from McNally 1965, plus a random constant term to provide some randomness
        f(x) = 33.6893424629 + 47.7907633437 * log10(x) + rand(-5:5)
        return f(mass)
    end
    export _stellar_rotation_velocity

    """
    Generates the temperature, mass, and radius of the star based on the class and subclass
    """
    function _stellarpredefs(class::String, subclass::String)
        temp::BigFloat = 0
        solar_masses::BigFloat = 0
        solar_radii::BigFloat = 0

        # Ranges are all taken from wikipedia
        if class == "O"
            # No observational data for O0V or O1V class stars, so there is nothing to really go off of
            if subclass == "0" || "1" || "2" || "3"
                temp = BigFloat(rand(44900:55200))
                solar_masses = BigFloat(rand(120:152))
                solar_radii = BigFloat(rand(15:25))

            elseif subclass == "4"
                temp = BigFloat(rand(42900:44900))
                solar_masses = BigFloat(rand(85.31:120))
                solar_radii = BigFloat(rand(13.43:15))
                
            elseif subclass == "5"
                temp = BigFloat(rand(41400:42900))
                solar_masses = BigFloat(rand(60:85.31))
                solar_radii = BigFloat(rand(12:13.43))

            elseif subclass == "6"
                temp = BigFloat(rand(39500:41400))
                solar_masses = BigFloat(rand(43.71:60))
                solar_radii = BigFloat(rand(10.71:12))

            elseif subclass == "7"
                temp = BigFloat(rand(37100:39500))
                solar_masses = BigFloat(rand(30.85:43.71))
                solar_radii = BigFloat(rand(9.52:10.71))

            elseif subclass == "8"
                temp = BigFloat(rand(35100:37100))
                solar_masses = BigFloat(rand(23:30.85))
                solar_radii = BigFloat(rand(8.5:9.52))

            else
                temp = BigFloat(rand(33300:35100))
                solar_masses = BigFloat(rand(19.63:23))
                solar_radii = BigFloat(rand(7.51:8.5))

            end

        elseif class == "B"
            if subclass == "0"
                temp = BigFloat(rand(31400:33300))
                solar_masses = BigFloat(rand(17.70:19.63))
                solar_radii = BigFloat(rand(7.16:7.51))

            elseif subclass == "1"
                temp = BigFloat(rand(26000:31400))
                solar_masses = BigFloat(rand(11.00:17.70))
                solar_radii = BigFloat(rand(5.71:7.16))

            elseif subclass == "2"
                temp = BigFloat(rand(20600:26000))
                solar_masses = BigFloat(rand(7.30:11.00))
                solar_radii = BigFloat(rand(4.06:5.71))

            elseif subclass == "3"
                temp = BigFloat(rand(17000:20600))
                solar_masses = BigFloat(rand(5.40:7.30))
                solar_radii = BigFloat(rand(3.61:4.06))

            elseif subclass == "4"
                temp = BigFloat(rand(16400:17000))
                solar_masses = BigFloat(rand(5.10:5.40))
                solar_radii = BigFloat(rand(3.46:3.61))

            elseif subclass == "5"
                temp = BigFloat(rand(15700:16400))
                solar_masses = BigFloat(rand(4.70:5.10))
                solar_radii = BigFloat(rand(3.36:3.46))

            elseif subclass == "6"
                temp = BigFloat(rand(14500:15700))
                solar_masses = BigFloat(rand(4.30:4.70))
                solar_radii = BigFloat(rand(3.27:3.36))

            elseif subclass == "7"
                temp = BigFloat(rand(14000:14500))
                solar_masses = BigFloat(rand(3.92:4.30))
                solar_radii = BigFloat(rand(2.94:3.27))

            elseif subclass == "8"
                temp = BigFloat(rand(12300:14000))
                solar_masses = BigFloat(rand(3.38:3.92))
                solar_radii = BigFloat(rand(2.86:2.94))

            else
                temp = BigFloat(rand(10700:12300))
                solar_masses = BigFloat(rand(2.75:3.38))
                solar_radii = BigFloat(rand(2.49:2.86))

            end

        elseif class == "A"
            if subclass == "0"
                temp = BigFloat(rand(9700:10700))
                solar_masses = BigFloat(rand(2.18:2.75))
                solar_radii = BigFloat(rand(2.193:2.49))

            elseif subclass == "1"
                temp = BigFloat(rand(9300:9700))
                solar_masses = BigFloat(rand(2.05:2.18))
                solar_radii = BigFloat(rand(2.136:2.193))

            elseif subclass == "2"
                temp = BigFloat(rand(8800:9300))
                solar_masses = BigFloat(rand(1.98:2.05))
                solar_radii = BigFloat(rand(2.117:2.136))

            elseif subclass == "3"
                temp = BigFloat(rand(8600:9300))
                solar_masses = BigFloat(rand(1.93:1.98))
                solar_radii = BigFloat(rand(1.861:2.117))

            elseif subclass == "4"
                temp = BigFloat(rand(8250:8600))
                solar_masses = BigFloat(rand(1.88:1.93))
                solar_radii = BigFloat(rand(1.794:1.861))

            elseif subclass == "5"
                temp = BigFloat(rand(8100:8250))
                solar_masses = BigFloat(rand(1.86:1.88))
                solar_radii = BigFloat(rand(1.785:1.794))

            elseif subclass == "6"
                temp = BigFloat(rand(7910:8100))
                solar_masses = BigFloat(rand(1.83:1.86))
                solar_radii = BigFloat(rand(1.775:1.785))

            elseif subclass == "7"
                temp = BigFloat(rand(7760:7910))
                solar_masses = BigFloat(rand(1.81:1.83))
                solar_radii = BigFloat(rand(1.750:1.775))

            elseif subclass == "8"
                temp = BigFloat(rand(7590:7760))
                solar_masses = BigFloat(rand(1.77:1.81))
                solar_radii = BigFloat(rand(1.7479:1.750))

            else
                temp = BigFloat(rand(7400:7590))
                solar_masses = BigFloat(rand(1.75:1.77))
                solar_radii = BigFloat(rand(1.747:1.7479))

            end

        elseif class == "F"
            if subclass == "0"
                temp = BigFloat(rand(7220:7400))
                solar_masses = BigFloat(rand(1.61:1.75))
                solar_radii = BigFloat(rand(1.728:1.747))

            elseif subclass == "1"
                temp = BigFloat(rand(7070:7220))
                solar_masses = BigFloat(rand(1.50:1.61))
                solar_radii = BigFloat(rand(1.679:1.728))

            elseif subclass == "2"
                temp = BigFloat(rand(6820:7070))
                solar_masses = BigFloat(rand(1.46:1.50))
                solar_radii = BigFloat(rand(1.622:1.679))

            elseif subclass == "3"
                temp = BigFloat(rand(6750:6820))
                solar_masses = BigFloat(rand(1.44:1.46))
                solar_radii = BigFloat(rand(1.578:1.622))

            elseif subclass == "4"
                temp = BigFloat(rand(6670:6750))
                solar_masses = BigFloat(rand(1.38:1.44))
                solar_radii = BigFloat(rand(1.533:1.578))

            elseif subclass == "5"
                temp = BigFloat(rand(6550:6750))
                solar_masses = BigFloat(rand(1.33:1.38))
                solar_radii = BigFloat(rand(1.473:1.533))

            elseif subclass == "6"
                temp = BigFloat(rand(6350:6550))
                solar_masses = BigFloat(rand(1.25:1.33))
                solar_radii = BigFloat(rand(1.359:1.473))

            elseif subclass == "7"
                temp = BigFloat(rand(6280:6350))
                solar_masses = BigFloat(rand(1.21:1.25))
                solar_radii = BigFloat(rand(1.324:1.359))

            elseif subclass == "8"
                temp = BigFloat(rand(6180:6280))
                solar_masses = BigFloat(rand(1.18:1.21))
                solar_radii = BigFloat(rand(1.221:1.324))

            else
                temp = BigFloat(rand(6050:6180))
                solar_masses = BigFloat(rand(1.13:1.18))
                solar_radii = BigFloat(rand(1.167:1.221))

            end

        elseif class == "G"
            if subclass == "0"
                temp = BigFloat(rand(5930:6050))
                solar_masses = BigFloat(rand(1.06:1.13))
                solar_radii = BigFloat(rand(1.100:1.167))

            elseif subclass == "1"
                temp = BigFloat(rand(5860:5930))
                solar_masses = BigFloat(rand(1.03:1.06))
                solar_radii = BigFloat(rand(1.060:1.100))

            elseif subclass == "2"
                temp = BigFloat(rand(5770:5860))
                solar_masses = BigFloat(rand(1.00:1.03))
                solar_radii = BigFloat(rand(1.012:1.060))

            elseif subclass == "3"
                temp = BigFloat(rand(5720:5770))
                solar_masses = BigFloat(rand(0.99:1.00))
                solar_radii = BigFloat(rand(1.002:1.012))

            elseif subclass == "4"
                temp = BigFloat(rand(5680:5720))
                solar_masses = BigFloat(rand(0.985:0.99))
                solar_radii = BigFloat(rand(0.991:1.002))

            elseif subclass == "5"
                temp = BigFloat(rand(5660:5680))
                solar_masses = BigFloat(rand(0.98:0.985))
                solar_radii = BigFloat(rand(0.977:0.991))

            elseif subclass == "6"
                temp = BigFloat(rand(5600:5660))
                solar_masses = BigFloat(rand(0.97:0.98))
                solar_radii = BigFloat(rand(0.949:0.977))

            elseif subclass == "7"
                temp = BigFloat(rand(5550:5600))
                solar_masses = BigFloat(rand(0.95:0.97))
                solar_radii = BigFloat(rand(0.927:0.949))

            elseif subclass == "8"
                temp = BigFloat(rand(5480:5550))
                solar_masses = BigFloat(rand(0.94:0.95))
                solar_radii = BigFloat(rand(0.914:0.927))

            else
                temp = BigFloat(rand(5380:5480))
                solar_masses = BigFloat(rand(0.90:0.94))
                solar_radii = BigFloat(rand(0.853:0.914))

            end

        elseif class == "K"
            if class == "0"
                temp = BigFloat(rand(5270:5380))
                solar_masses = BigFloat(rand(0.88:0.90))
                solar_radii = BigFloat(rand(0.813:0.853))

            elseif subclass == "1"
                temp = BigFloat(rand(5170:5270))
                solar_masses = BigFloat(rand(0.86:0.88))
                solar_radii = BigFloat(rand(0.797:0.813))

            elseif subclass == "2"
                temp = BigFloat(rand(5100:5170))
                solar_masses = BigFloat(rand(0.82:0.86))
                solar_radii = BigFloat(rand(0.783:0.797))

            elseif subclass == "3"
                temp = BigFloat(rand(4830:5100))
                solar_masses = BigFloat(rand(0.78:0.82))
                solar_radii = BigFloat(rand(0.755:0.783))

            elseif subclass == "4"
                temp = BigFloat(rand(4600:4830))
                solar_masses = BigFloat(rand(0.73:0.78))
                solar_radii = BigFloat(rand(0.713:0.755))

            elseif subclass == "5"
                temp = BigFloat(rand(4440:4600))
                solar_masses = BigFloat(rand(0.70:0.73))
                solar_radii = BigFloat(rand(0.701:0.713))

            elseif subclass == "6"
                temp = BigFloat(rand(4300:4440))
                solar_masses = BigFloat(rand(0.69:0.70))
                solar_radii = BigFloat(rand(0.669:0.701))

            elseif subclass == "7"
                temp = BigFloat(rand(4100:4300))
                solar_masses = BigFloat(rand(0.64:0.69))
                solar_radii = BigFloat(rand(0.630:0.669))

            elseif subclass == "8"
                temp = BigFloat(rand(3990:4100))
                solar_masses = BigFloat(rand(0.62:0.64))
                solar_radii = BigFloat(rand(0.615:0.630))

            else
                temp = BigFloat(rand(3930:3990))
                solar_masses = BigFloat(rand(0.59:0.62))
                solar_radii = BigFloat(rand(0.608:0.615))

            end

        else # M Class Stars
            if class == "0"
                temp = BigFloat(rand(3850:3930))
                solar_masses = BigFloat(rand(0.57:0.59))
                solar_radii = BigFloat(rand(0.588:0.608))

            elseif subclass == "1"
                temp = BigFloat(rand(3660:3850))
                solar_masses = BigFloat(rand(0.50:0.57))
                solar_radii = BigFloat(rand(0.501:0.588))

            elseif subclass == "2"
                temp = BigFloat(rand(3560:3660))
                solar_masses = BigFloat(rand(0.44:0.50))
                solar_radii = BigFloat(rand(0.446:0.501))

            elseif subclass == "3"
                temp = BigFloat(rand(3430:3560))
                solar_masses = BigFloat(rand(0.37:0.44))
                solar_radii = BigFloat(rand(0.361:0.446))

            elseif subclass == "4"
                temp = BigFloat(rand(3210:3430))
                solar_masses = BigFloat(rand(0.23:0.37))
                solar_radii = BigFloat(rand(0.274:0.361))

            elseif subclass == "5"
                temp = BigFloat(rand(3060:3210))
                solar_masses = BigFloat(rand(0.162:0.23))
                solar_radii = BigFloat(rand(0.196:0.274))

            elseif subclass == "6"
                temp = BigFloat(rand(2810:3060))
                solar_masses = BigFloat(rand(0.102:0.162))
                solar_radii = BigFloat(rand(0.137:0.196))

            elseif subclass == "7"
                temp = BigFloat(rand(2680:2810))
                solar_masses = BigFloat(rand(0.090:0.102))
                solar_radii = BigFloat(rand(0.120:0.137))

            elseif subclass == "8"
                temp = BigFloat(rand(2570:2680))
                solar_masses = BigFloat(rand(0.085:0.090))
                solar_radii = BigFloat(rand(0.114:0.120))

            else
                temp = BigFloat(rand(2380:2570))
                solar_masses = BigFloat(rand(0.079:0.085))
                solar_radii = BigFloat(rand(0.102:0.114))
            end
        end
        return temp, solar_masses * Constants.SolarMass, solar_radii * Constants.SolarRadius
    end

    # WIP. Does not work
    function _stellardimensions(radius::BigFloat)
        equatorial_circumference::BigFloat = Solar.circumference(radius)
        flattening::BigFloat = 0
        surface_area::BigFloat = Solar.sphere_surface_area(radius)
        volume::BigFloat = Solar.sphere_volume(radius)

        return equatorial_circumference, flattening, surface_area, volume
    end

    """
    Generates a star and returns a Star
    """
    function generatestar()::Star
        main_class, subclass = _stellarclass()
        stellar_temperature, stellar_mass, stellar_radius = _stellarpredefs(main_class, subclass)
        circumference, flattening, surface_area, volume = _stellardimensions(stellar_radius)
        mean_density = Solar.mean_density(stellar_radius, volume)
        dimensions = Solar.StellarDimensions(stellar_radius, circumference, flattening, surface_area, volume)
        internals = Solar.StellarInternals(stellar_mass, mean_density)
        gravity = Solar.StellarGravity(Solar.surface_gravity(stellar_mass, stellar_radius), Solar.escape_velocity(stellar_mass, stellar_radius), Solar.rochelimit(stellar_radius, mean_density, BigFloat(1.0)))
        output = Solar.StellarOutput(stellar_temperature, Solar.luminosity(surface_area, stellar_temperature), BigFloat(0), Solar.radiance(stellar_temperature))
        rotation = Solar.StellarRotation(0, 0)
        misc = Solar.StellarMiscellaneous()
        bodies = Vector{Solar.Planet}(undef, 0)
        return Star("WIP Star", main_class, subclass, dimensions, internals, gravity, output, rotation, misc, bodies)
    end
    export generatestar

    function _gasgiant()

    end

    function _rockyplanet()
    end

    function generateplanet(M::BigFloat)
        #=
        Generate 
        =#
    end
end