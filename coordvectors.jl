module CoordVector
    #=
    https://en.wikipedia.org/wiki/Coordinate_vector
    https://en.wikipedia.org/wiki/Norm_(mathematics)
    https://physics.nist.gov/cuu/Constants/
    =#

    struct Vec2
        x
        y
    end

    struct Vec3
        x
        y
        z
    end

    """
    Returns the norm of a 2D vector v
    """
    function normalize(v::Vec2)
        return sqrt((v.x)^2 + (v.y)^2)
    end

    """
    Returns the norm of a 3D vector v
    """
    function normalize(v::Vec3)
        return sqrt((v.x)^2 + (v.y)^2 + (v.z)^2)
    end

    """
    Returns the unit vector of a 2D vector v
    """
    function unit(v::Vec2)
        norm = normalize(v)
        x, y = v.x / norm, v.y / norm
        return Vec2(x, y)
    end

    """
    Returns the unit vector of a 3D vector v
    """
    function unit(v::Vec3)
        norm = normalize(v)
        x, y, z = v.x / norm, v.y / norm, v.z / norm
        return Vec3(x, y, z)
    end

    """
    Returns the magnitude of the cross product of two vectors v := [a, b, 0], u := [c, d, 0]
    """
    function cross(v::Vec2, u::Vec2)
        return u.y*v.x - u.x*v.y
    end

    """
    Returns the cross product of two vectors v, u
    """
    function cross(v::Vec3, u::Vec3)
        x = *(u.z, v.y) - *(u.y, v.z)
        y = *(u.x, v.z) - *(u.z, v.x)
        z = *(u.y, v.x) - *(u.x, v.y)
        return Vec3(x, y, z)
    end

    """
    Returns the dot product of two vectors v, u
    """
    function dot(v::Vec2, u::Vec2)
        return *(v.x, u.x) + *(v.y, u.y)
    end

    """
    Returns the dot product of two vectors v, u
    """
    function dot(v::Vec3, u::Vec3)
        return *(v.x, u.x) + *(v.y, u.y) + *(v.z, u.z)
    end

    """
    Returns the magnitude of a vector v
    """
    function magnitude(v::Vec2)
        return sqrt((v.x)^2 + (v.y)^2)
    end

    """
    Returns the magnitude of a vector v
    """
    function magnitude(v::Vec3)
        return sqrt((v.x)^2 + (v.y)^2 + (v.z)^2)
    end
end # End module