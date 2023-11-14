module Solar
    include("astrovector.jl")
    using .Astrovector
    include("constants.jl")
    using .Constants
    include("utils.jl")
    using .Utils
    include("coordvectors.jl")
    using .CoordVector

    """
    Stores the six elements necessary to uniquely identify an object's orbit

    ...
    # Fields
    - `orbitalEccentricity`: The orbital eccentricity of the object's orbit
    - `semiMajorAxis`: The semi-major axis of the body's orbit
    - `inclination`: The inclination of the body's orbit relative to some reference plane
    - `longitudeOfTheAscendingNode`: The longitude of the ascending node of the orbit
    - `argumentOfPeriapsis`: The argument of periapsis of the orbit
    - `trueAnomaly`: The true anomaly of the orbit
    ...
    """
    struct KeplerianElements
        orbitalEccentricity::BigFloat
        semiMajorAxis::BigFloat
        inclination::BigFloat
        longitudeOfTheAscendingNode::BigFloat
        argumentOfPeriapsis::BigFloat
        trueAnomaly::BigFloat
    end
    export KeplerianElements

    """
    ...
    # Fields
    - `apoapsis`: The apoapsis of the object's orbit
    - `periapsis`: The periapsis of the object's orbit
    - `semiMajorAxis`: The semi-major axis of the object's orbit
    - `semiMinorAxis`: The semi-minor axis of the object's orbit
    - `semiLatusRectum`: The semi-latus rectum of the object's orbit
    - `eccentricity`: The eccentricity of the object's orbit
    ...
    """
    struct OrbitalPath
        apoapsis::BigFloat
        periapsis::BigFloat
        semiMajorAxis::BigFloat
        semiMinorAxis::BigFloat
        semiLatusRectum::BigFloat
        eccentricity::BigFloat
    end
    export OrbitalPath

    """
    ...
    Stores temporal information about the object's orbit

    # Fields
    - `orbitalPeriod`: The orbital period of the object around its central body
    - `averageOrbitalSpeed`: The mean orbital speed of the object
    - `periapsisOrbitalSpeed`: The orbital speed of the object at its periapsis
    - `apoapsisOrbitalSpeed`: The orbital speed of the object at its apoapsis
    ...
    """
    struct OrbitalTime
        orbitalPeriod::BigFloat
        averageOrbitalSpeed::BigFloat
        periapsisOrbitalSpeed::BigFloat
        apoapsisOrbitalSpeed::BigFloat
    end
    export OrbitalTime

    """
    Stores angular information about an object's orbit

    ...
    # Fields
    - `meanAnomaly`: The mean anomaly of the orbit
    - `inclination`: The inclination of the orbit from some reference plane
    - `longitudeOfTheAscendingNode`: The longitude of the ascending node of the planet
    - `argumentOfPeriapsis`: The argument of periapsis of the orbit
    ...
    """
    struct OrbitalAngles
        meanAnomaly::BigFloat
        inclination::BigFloat
        longitudeOfTheAscendingNode::BigFloat
        argumentOfPeriapsis::BigFloat
    end
    export OrbitalAngles

    """
    Stores dimensional information about a planet
    ...
    # Fields
    - `equatorialRadius`: The equatorial radius of the planet
    - `meanRadius`: The mean radius of the planet
    - `polarRadius`: The polar radius of the planet
    - `flattening`: The flattening of the planet
    - `equatorialCircumference`: The equatorial circumference of the planet
    - `meanCircumference`: The mean circumference of the planet
    - `polarCircumference`: The polar circumference of the planet
    - `surfaceArea`: The surface area of the planet
    - `volume`: The volume of the planet
    ...
    """
    struct PlanetaryDimensions
        equatorialRadius::BigFloat
        meanRadius::BigFloat
        polarRadius::BigFloat
        flattening::BigFloat
        equatorialCircumference::BigFloat
        meanCircumference::BigFloat
        polarCircumference::BigFloat
        surfaceArea::BigFloat
        volume::BigFloat
    end
    export PlanetaryDimensions

    """
    ...
    # Fields
    - `mass`: The mass of the planet
    - `meanDensity`: The mean density of the planet
    ...
    """
    struct PlanetaryInternals
        mass::BigFloat
        meanDensity::BigFloat
    end
    export PlanetaryInternals

    """
    Stores gravitational information about a planet

    ...
    # Fields
    - `equatorialSurfaceGravity`: The surface gravity of the planet at the equator
    - `polarSurfaceGravity`: The surface gravity at the planet's poles
    - `meanSurfaceGravity`: The mean surface gravity of the planet
    - `escapeVelocity`: The escape velocity of the planet
    ...
    """
    struct PlanetaryGravity
        equatorialSurfaceGravity::BigFloat
        polarSurfaceGravity::BigFloat
        meanSurfaceGravity::BigFloat
        escapeVelocity::BigFloat
    end
    export PlanetaryGravity

    """
    Stores rotational information about a planet

    ...
    # Fields
    - `siderealRotationPeriod`: The sidereal rotation period of the planet
    - `rotationVelocity`: The rotation velocity of the planet
    - `axialTilt`: The axial tilt of the planet
    ...
    """
    struct PlanetaryRotation
        siderealRotationPeriod::BigFloat
        rotationVelocity::BigFloat
        axialTilt::BigFloat
    end
    export PlanetaryRotation

    """
    Stores output related information about the planet

    ...
    # Fields
    - `albedo`: The albedo of the planet's surface
    - `blackBodyTemperature`: The idealized blackbody temperature of the planet
    ...
    """
    struct PlanetaryOutput
        albedo::BigFloat
        blackBodyTemperature::BigFloat
    end
    export PlanetaryOutput

    """
    Stores uncategorized information about a planet
    ...
    # Fields
    - ``:
    ...
    """
    struct PlanetaryMiscellaneous

    end

    """
    ...
    # Fields
    - `path`: The path of the object's orbit
    - `time`: Any temporal elements of the object's orbit
    - `angles`: Any angular elements of the object's orbit
    ...
    """
    struct Orbit
        path::OrbitalPath
        time::OrbitalTime
        angles::OrbitalAngles
    end
    export Orbit
    """
    Stores all non-orbital information about a planet
    ...
    # Fields
    - `dimensions`: Dimensional information about a planet``
    - `internals`: Internal information about a planet
    - `gravity`: Gravitational information about a planet
    - `rotation`: Rotational information about a planet
    - `output`: Output information about a planet
    - `miscellaneous`: Uncategorized information about a planet
    ...
    """
    struct Planetary
        dimensions::PlanetaryDimensions
        internals::PlanetaryInternals
        gravity::PlanetaryGravity
        rotation::PlanetaryRotation
        output::PlanetaryOutput
        miscellaneous::PlanetaryMiscellaneous
    end
    export Planetary

    """
    Stores all information about a planet

    ...
    # Fields
    - `name`: The name of the planet
    - `orbit`: Stores the orbital information about the planet
    - `physical`: Stores information about the planet itself
    ...
    """
    struct Planet
        name::String
        orbit::Orbit
        physical::Planetary
    end
    export Planet

    """
    Stores the dimensions of a star

    ...
    # Fields
    - `equatoralRadius`: The equatoral radius of the star
    - `equatorialCircumference`: The equatorial circumference of the star
    - `flattening`: The flattening of the star
    - `surfaceArea`: The surface area of the star
    - `volume`: The volume of the star
    ...
    """
    struct StellarDimensions
        equatorialRadius::BigFloat
        equatorialCircumference::BigFloat
        flattening::BigFloat
        surfaceArea::BigFloat
        volume::BigFloat
    end
    export StellarDimensions

    """
    Stores the internal information about the star

    ...
    # Fields
    - `mass`: The mass of the star
    - `averageDensity`: The average density of the star
    ...
    """
    struct StellarInternals
        mass::BigFloat
        averageDensity::BigFloat
    end
    export StellarInternals

    """
    Stores the gravitational information about the star

    ...
    # Fields
    - `equatorialSurfaceGravity`: The surface gravity of the star at the equator
    - `escapeVelocity`: The escape velocity of the star, specifically at the equator
    - `rocheLimit`: The roche limit of the star
    ...
    """
    struct StellarGravity
        equatorialSurfaceGravity::BigFloat
        escapeVelocity::BigFloat
        rocheLimit::BigFloat
    end
    export StellarGravity

    """
    Stores the output related information about the star

    ...
    # Fields
    - `temperature`: The temperature of the star
    - `luminosity`: The luminosity of the star
    - `colorBV`: The color of the star on the B-V color index
    - `meanRadiance`: The mean radiance of the star
    ...
    """
    struct StellarOutput
        temperature::BigFloat
        luminosity::BigFloat
        colorBV::BigFloat
        meanRadiance::BigFloat
    end
    export StellarOutput

    """
    Stores rotational information about a star

    ...
    # Fields
    - `siderealRotationPeriod`: The sidereal rotation period of the star
    - `rotationVelocity`: The rotation velocity of the star
    ...
    """
    struct StellarRotation
        siderealRotationPeriod::BigFloat
        rotationVelocity::BigFloat
    end
    export StellarRotation

    """
    Stores uncategorized information about the star
    ...
    # Fields
    - ``:
    ...
    """
    struct StellarMiscellaneous
        
    end

    """
    Contains all information about the star

    ...
    # Fields
    - `name`: The name of the star
    - `dimensions`: The dimensional information about the star
    - `internals`: The internal information about the star
    - `gravity`: The gravitational information about the star
    - `output`: The output information about the star
    - `rotation`: The rotational information about the star
    ...
    """
    struct Star
        name::String
        class::String
        subclass::String
        dimensions::StellarDimensions
        internals::StellarInternals
        gravity::StellarGravity
        output::StellarOutput
        rotation::StellarRotation
        miscellaneous::StellarMiscellaneous
        bodies::Vector{Planet}
    end
    export Star

    function printstar(star)
        name = star.name
        class = string(star.class, star.subclass)
        dimensions = star.dimensions
        internals = star.internals
        gravity = star.gravity
        output = star.output
        rotation = star.rotation
        # Deprecated for now
        #miscellaneous = Star.miscellaneous
        #bodies = Star.bodies

        println("Name of star: ", name)
        println("Spectral classification: ", class)
        println("Equatorial Radius: ", round(dimensions.equatorialRadius, sigdigits = 5))
        println("Equatorial Circumference: ", round(dimensions.equatorialCircumference, sigdigits = 5))
        println("Flattening: ", round(dimensions.flattening, sigdigits = 5))
        println("Surface area: ", round(dimensions.surfaceArea, sigdigits = 5))
        println("Volume: ", round(dimensions.volume, sigdigits = 5))
        println("Mass: ", round(internals.mass, sigdigits = 5))
        println("Average density: ", round(internals.averageDensity, sigdigits = 5))
        println("Surface gravity: ", round(gravity.equatorialSurfaceGravity, sigdigits = 5))
        println("Escape velocity: ", round(gravity.escapeVelocity, sigdigits = 5))
        println("Surface Temperature: ", round(output.temperature, sigdigits = 5))
        println("Luminosity: ", round(output.luminosity, sigdigits = 5))
        println("Color (B-V): ", round(output.colorBV, sigdigits = 5))
        println("Radiance: ", round(output.meanRadiance, sigdigits = 5))
        println("Sidereal rotation period: ", round(rotation.siderealRotationPeriod, sigdigits = 5))
        println("Rotation velocity: ", round(rotation.rotationVelocity, sigdigits = 5))
    end
    export printstar
    """
    Returns the semi-minor axis of an ellipse with semi-major axis a and eccentricity e

    ...
    # Arguments
    - `a`: Semi-major Axis
    - `e`: Eccentricity
    ...
    """
    function semi_minor_axis(a::BigFloat, e::BigFloat)
        return a * sqrt(1 - e*e)
    end

    """
    Returns the semi-latus rectum of an ellipse with semi-major axis a and eccentricity e

    ...
    # Arguments
    - `a`: Semi-major axis
    - `e`: Eccentricity
    ...
    """
    function semi_latus_rectum(a::BigFloat, e::BigFloat)
        return a * (1 - e*e)
    end

    """
    Returns the apoapsis of an ellipse

    ...
    # Arguments
    - `a`: Semi-major axis
    - `e`: Eccentricity
    ...
    """
    function apoapsis(a::BigFloat, e::BigFloat)
        return a * (1 + e)
    end

    """
    Returns the periapsis of an ellipse

    ...
    # Arguments
    - `a`: Semi-major axis
    - `e`: Eccentricity
    ...
    """
    function periapsis(a::BigFloat, e::BigFloat)
        return a * (1 - e)
    end

    """
    Returns the orbital period of a body

    ...
    # Arguments
    - `a`: Semi-major axis of the orbit
    - `μ`: The standard gravitational parameter
    ...
    """
    function orbital_period(a::BigFloat, μ::BigFloat)
        return 2 * π * sqrt(/(a*a*a, μ))
    end

    """
    Returns the mean orbital speed of an object orbiting a central body
    
    ...
    # Arguments
    - `semi`: Semi-major Axis
    - `apo`: apoapsis
    - `peri`: periapsis
    - `μ`: Standard Gravitational Parameter
    - `N`: The number of sections to divide the orbit into
    ...
    """
    function mean_orbital_speed(semi::BigFloat, apo::BigFloat, peri::BigFloat, μ::BigFloat, N=64)
        incrementation = (apo - peri) / N
        summation = 0
        for i = 1:N
            r = peri + i * incrementation
            summation += instantaneous_orbital_speed(semi, r, μ)
        end

        return summation / N
    end

    """
    Returns the instantaneous orbital speed of an object with semi major axis a, at a current distance of r from the central body
    
    ...
    # Arguments
    - `a`: Semi-major Axis of orbit
    - `r`: Current distance from central body
    - `μ`: Standard gravitational parameter
    ...
    """
    function instantaneous_orbital_speed(a::BigFloat, r::BigFloat, μ::BigFloat)
        return sqrt(μ * ((2 / r) - (1 / a)))
    end

    """
    Finds an approximation for the flattening of a body
    
    ...
    # Arguments
    - `ω`: The rotation velocity of the body
    - `a`: The semi-major axis of the body
    - `M`: The mass of the body
    ...
    """
    function planetary_flattening(ω::BigFloat, a::BigFloat, m::BigFloat)
        return /(5 * ω*ω * a*a*a, 4 * GravitationalConstant * m)
    end

    """
    Converts the flattening of an ellipse to the eccentricity

    ...
    # Arguments
    - `f`: The flattening of the ellipse
    ...
    """
    function flattening_to_eccentricity(f::BigFloat)
        return sqrt(2*f - f^2)
    end

    """
    Returns the mean radius of an ellipse

    ...
    # Arguments
    - `a`: Semi-major axis
    - `b`: Semi-minor axis
    ...
    """
    function mean_radius(a::BigFloat, b::BigFloat)
        return (a*b) / agm(a, b)
    end

    """
    Returns the surface area of an oblate spheroid (e.g., a planet)

    ...
    # Arguments
    - `a`: The semi-major axis of the object
    - `b`: The semi-minor axis of the object
    - `e`: The eccentricity of the object 
    ...
    """
    function spheroid_surface_area(a::BigFloat, b::BigFloat, e::BigFloat)
        constant_term = 2 * π * a * a
        log_factor = /(π * b*b, e)
        log_term = log(/(1 + e, 1 - e))
        return constant_term + log_factor * log_term
    end

    """
    Returns the volume of a spheroid

    ...
    # Arguments
    - `a`: Semi-major axis
    - `b`: Semi-minor axis
    ...
    """
    function spheroid_volume(a::BigFloat, b::BigFloat)
        return (4 * π * a * a * b) / 3
    end

    """
    Returns the surface area of a sphere

    ...
    # Arguments
    - `r`: The radius of the sphere
    ...
    """
    function sphere_surface_area(r::BigFloat)
        return 4 * π * r * r
    end

    """
    Returns the volume of a sphere

    ...
    # Arguments
    - `r`: The radius of the sphere
    ...
    """
    function sphere_volume(r::BigFloat)
        return /(4 * π * r * r * r, 3)
    end

    """
    Returns the circumference of a sphere

    ...
    # Arguments
    - `r`: The radius of the object
    ...
    """
    function circumference(r::BigFloat)
        return 2 * π * r
    end

    """
    Returns the mean density of an object

    ...
    # Arguments
    - `m`: The mass of the object
    - `V`: The volume of the object
    ...
    """
    function mean_density(m::BigFloat, V::BigFloat)
        return m / V
    end

    """
    Returns the surface gravity at distance r from the center of a planet

    ...
    # Arguments
    - `m`: The mass of the body
    - `r`: The distance from the center of the body
    ...
    """
    function surface_gravity(m::BigFloat, r::BigFloat)
        return /(GravitationalConstant * m, r*r)
    end

    """
    Returns the required escape velocity of an object from a central body with mass m at a distance r from the center of the body.
    
    ...
    # Arguments
    - `m`: The mass of the body
    - `r`: The distance of the object from the center of the body
    ...
    """
    function escape_velocity(m::BigFloat, r::BigFloat)
        return sqrt(/(2*GravitationalConstant*m, r))
    end

    """
    Returns the sidereal rotation period of a body

    ...
    # Arguments
    - `C`: The circumference of the body
    - `v`: The rotation velocity of the body
    ...
    """
    function sidereal_rotation_period(C::BigFloat, v::BigFloat)
        return C / v
    end

    """
    Returns the irradiance/flux density of a star

    ...
    # Arguments
    - `L`: The luminosity/radiant flux of the star
    - `A`: The surface area of the star
    ...
    """
    function irradiance(L::BigFloat, A::BigFloat)
        return L / A
    end

    """
    Returns the radiance of a star

    ...
    # Arguments
    - `T`: The temperature of the star
    ...
    """
    function radiance(T::BigFloat)
        return /(StefanBoltzmannConstant * T * T * T * T, π)
    end
    
    """
    Returns the luminosity of a star

    ...
    # Arguments
    - `A`: Surface area of the star
    - `T`: Temperature of the star
    ...
    """
    function luminosity(A::BigFloat, T::BigFloat)
        return StefanBoltzmannConstant * A * T * T * T * T
    end
    
    """
    Returns the specific relative angular momentum of an object

    ...
    # Arguments
    - `e`: The eccentricity of the object's orbit
    - `μ`: The standard gravitational parameter of the object
    - `a`: The semi-major axis of the object's orbit
    ...
    """
    function specific_relative_angular_momentum(e::BigFloat, μ::BigFloat, a::BigFloat)
        return sqrt(μ*a * (1 - e)^2)
    end
    """
    Returns the specific kinetic energy of an object

    ...
    # Arguments
    - `v`: The velocity of the object
    ...
    """
    function specific_kinetic_energy(v::BigFloat)
        return (v*v) / 2
    end
    
    """
    Returns the specific potential energy of an object

    ...
    # Arguments
    - `μ`: The standard gravitational parameter
    - `r`: The distance of the object from the central body
    ...
    """
    function specific_potential_energy(μ::BigFloat, r::BigFloat)
        return -μ / r
    end
    
    """
    Returns the specific orbital energy of an object

    ...
    # Arguments
    - `μ`: The standard gravitational parameter
    - `a`: The semi-major axis of the object's orbit
    ...
    """
    function specific_orbital_energy(μ::BigFloat, a::BigFloat)
        return -μ / (2 * a)
    end
    
    """
    Returns the mean motion of an object

    ...
    # Arguments
    - `T`: The orbital period of the object
    ...
    """
    function mean_motion()
        return /(2*π, T)
    end

    """
    Returns the Roche limit of a body

    ...
    # Arguments
    - `R`: The radius of the primary body
    - `Ρ`: The density of the primary body 
    - `ρ`: The density of the secondary body
    ...
    """
    function rochelimit(R::BigFloat, Ρ::BigFloat, ρ::BigFloat)
        return R * cbrt(2 * /(Ρ, ρ))
    end

    """
    Returns the sphere of influence of a body

    ...
    # Arguments
    - `a`: The semi-major axis of the orbit of the smaller body
    - `M`: The mass of the larger object
    - `m`: The mass of the smaller object
    ...
    """
    function sphereofinfluence(a::BigFloat, M::BigFloat, m::BigFloat)
        return a * (m / M)^(2/5)
    end

    """
    Returns the hill sphere of a body

    ...
    # Arguments
    - `a`: The semi-major axis of the smaller body's orbit
    - `e`: The eccentricity of the smaller body's orbit
    - `M`: The mass of the larger body
    - `m`: The mass of the smaller body
    ...
    """
    function hillsphere(a::BigFloat, e::BigFloat, M::BigFloat, m::BigFloat)
        return periapsis(a, e) * cbrt(m / (3 * (M + m)))
    end

    """
    Defines the path that an orbiting body takes around a central body in polar coordinates

    ...
    # Arguments
    - `ℓ`: Angular momentum of the orbiting body → mr²θ′
    - `m`: The mass of the body
    - `μ`: The standard gravitational parameter
    - `e`: The eccentricity of the body's orbit
    - `θ`: The true anomaly of the body's orbit
    ...
    """
    function orbit(ℓ, m, μ, e, θ)
        return /(ℓ*ℓ, μ*m*m) * /(1, 1 + e*cos(θ))
    end

    """
    Returns the coefficient of the inverse square law central force

    ...
    # Arguments
    - `M`: The mass of the central body
    - `m`: The mass of the smaller, orbiting body
    ...
    """
    function centralforcecoefficient(M::BigFloat, m::BigFloat)
        return m * GravitationalConstant * sqrt(M * (M + m))
    end
end # End Module