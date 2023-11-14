module Utils
    """
    Arithmetic mean two of numbers x, y
    """
    function arithmetic(x, y)
        return (x + y) / 2
    end
    export arithmetic

    """
    Arithmetic mean of N numbers a, b, c, ..., n
    """
    function arithmeticN(numbers)
        return sum(numbers) / length(numbers)
    end
    export arithmeticN

    """
    Geometric mean of two numbers x, y
    """
    function geometric(x, y)
        return sqrt(x * y)
    end
    export geometric

    """
    Geometric mean of N numbers a, b, c, ..., n
    """
    function geometricN(numbers)
        return (prod(numbers))^(1 / length(numbers))
    end
    export geometricN

    """
    Returns arith(x, y), geo(x, y)
    """
    function next_agm(x, y)
        return arithmetic(x, y), geometric(x, y)
    end

    """
    Arithmetic-geometric mean of numbers x, y

    N is the number of iterations to perform, defaulting to 10
    """
    function agm(x, y, N=10)
        a_n = arithmetic(x, y)
        g_n = geometric(x, y)
        for i = 1:N
            a_n, g_n = next_agm(a_n, g_n)
        end

        return a_n
    end
    export agm

    """
    Returns the harmonic mean of 2 numbers x and y
    """
    function harmonic(x, y)
        xinv, yinv = 1 / x, 1 / y
        return 2 / (xinv  + yinv)
    end
    export harmonic

    """
    Returns the harmonic mean of N numbers
    """
    function harmonicN(numbers)
        inverses = sum([(1 / i) for i in numbers])
        elements = length(numbers)
        return elements / inverses
    end
    export harmonicN

    """
    Implements Ramanujan's 1914 approximation of the circumference of an ellipse.
    
    a: Semi-major Axis

    b: Semi-minor Axis

    Does not provide good approximations of the Elliptic integral
    """
    function ellipsecircumference(a, b)
        apb = a+b
        amb = a-b
        lambda = 3 * /(amb^2, apb^2)
        return π*apb * /(1 + lambda, 10 + sqrt(4 - lambda))
    end
    export ellipsecircumference
end # End module