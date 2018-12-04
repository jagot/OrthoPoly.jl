module OrthoPoly
using Polynomials

macro recur_poly_gen(name, rule, p₁, extra_args...)
    generator = gensym()
    @eval begin
        if isempty($extra_args)
            ($name)(n::Integer, x = Poly([0, 1])) = ($generator)(n::Integer, x)[1]
        else
            ($name)(n::Integer, $(extra_args...), x = Poly([0, 1])) = ($generator)(n::Integer, x, $(extra_args...))[1]
        end
        function ($generator)(n::Integer, x::T, $(extra_args...)) where T
            if n == 0
                one(x), zero(x)
            elseif n == 1
                $p₁, one(x)
            else
                p₋₁, p₋₂ = ($generator)(n-1, x, $(extra_args...))
                $rule, p₋₁
            end
        end
    end
end

# Jacobi polynomials are orthogonal on the interval [-1,1] with respect
# to the weight (1-x)^α(1+x)^β. Legendre and Chebyshev polynomials are
# special cases.
# TODO: General Jacobi, Gegenbauer and associated Legendre polynomials
# (the latter are not polynomials for odd values of m).
@recur_poly_gen(chebyshevT, 2x*p₋₁ - p₋₂, x) # Chebyshev polynomials, first kind
@recur_poly_gen(chebyshevU, 2x*p₋₁ - p₋₂, 2x) # Chebyshev polynomials, second kind
@recur_poly_gen(legendre, ((2n - 1)x*p₋₁ - (n - 1)p₋₂)/n, x) # Legendre polynomials

# Hermite polynomials are orthogonal on the interval (-∞,∞) with
# weight exp(-x²/a) where a = 2 for probabilists and a = 1 for
# physicists.
@recur_poly_gen(hermite_prob, x*p₋₁ - (n - 1)p₋₂, x) # Probabilists' Hermite polynomials
@recur_poly_gen(hermite, 2x*p₋₁ - 2(n - 1)p₋₂, 2x) # Physicists' Hermite polynomials

# Associated Laguerre polynomials are orthogonal on the interval [0,∞)
# with respect to the weight x^α exp(-x), ∀ α ∈ ℝ.
@recur_poly_gen(laguerre, ((2n - 1 - x)p₋₁ - (n - 1)p₋₂)/n, 1 - x) # Laguerre polynomials
@recur_poly_gen(laguerre_assoc, ((2n - 1 + α - x)p₋₁ - (n - 1 + α)p₋₂)/n, 1 + α - x, α)

export chebyshevT, chebyshevU, legendre, laguerre, laguerre_assoc, hermite, hermite_prob

end # module
