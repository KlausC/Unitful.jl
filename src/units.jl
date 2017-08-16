@generated function *(a0::FreeUnits, a::FreeUnits...)

    # Sort the units uniquely. This is a generated function so that we
    # don't have to figure out the units each time.
    linunits = Vector{Unit}()

    for x in (a0, a...)
        xp = x.parameters[1]
        append!(linunits, xp[1:end])
    end

    # linunits is an Array containing all of the Unit objects that were
    # found in the type parameters of the FreeUnits objects (a0, a...)
    sort!(linunits, by=x->power(x))
    sort!(linunits, by=x->tens(x))
    sort!(linunits, by=x->name(x))

    # [m,m,cm,cm^2,cm^3,nm,m^4,µs,µs^2,s]
    # reordered as:
    # [nm,cm,cm^2,cm^3,m,m,m^4,µs,µs^2,s]

    # Collect powers of a given unit into `c`
    c = Vector{Unit}()
    if !isempty(linunits)
        i = start(linunits)
        oldstate = linunits[i]
        p = 0//1
        while !done(linunits, i)
            (state, i) = next(linunits, i)
            if tens(state) == tens(oldstate) && name(state) == name(oldstate)
                p += power(state)
            else
                if p != 0
                    push!(c, Unit{name(oldstate),dimtype(oldstate)}(tens(oldstate), p))
                end
                p = power(state)
            end
            oldstate = state
        end
        if p != 0
            push!(c, Unit{name(oldstate),dimtype(oldstate)}(tens(oldstate), p))
        end
    end
    # results in:
    # [nm,cm^6,m^6,µs^3,s]

    d = (c...)
    f = typeof(mapreduce(dimension, *, NoDims, d))
    :(FreeUnits{$d,$f}())
end
*(a0::ContextUnits, a::ContextUnits...) =
    ContextUnits(*(FreeUnits(a0), FreeUnits.(a)...),
                    *(FreeUnits(upreferred(a0)), FreeUnits.((upreferred).(a))...))
FreeOrContextUnits = Union{FreeUnits, ContextUnits}
*(a0::FreeOrContextUnits, a::FreeOrContextUnits...) =
    *(ContextUnits(a0), ContextUnits.(a)...)
*(a0::FixedUnits, a::FixedUnits...) =
    FixedUnits(*(FreeUnits(a0), FreeUnits.(a)...))

"""
```
*(a0::Units, a::Units...)
```

Given however many units, multiply them together. This is actually handled by
a few different methods, since we have `FreeUnits`, `ContextUnits`, and `FixedUnits`.

Collect [`Unitful.Unit`](@ref) objects from the type parameter of the
[`Unitful.Units`](@ref) objects. For identical units including SI prefixes
(i.e. cm ≠ m), collect powers and sort uniquely by the name of the `Unit`.
The unique sorting permits easy unit comparisons.

Examples:

```jldoctest
julia> u"kg*m/s^2"
kg m s^-2

julia> u"m/s*kg/s"
kg m s^-2

julia> typeof(u"m/s*kg/s") == typeof(u"kg*m/s^2")
true
```
"""
*(a0::Units, a::Units...) = FixedUnits(*(FreeUnits(a0), FreeUnits.(a)...))
# Logic above is that if we're not using FreeOrContextUnits, at least one is FixedUnits.

/(x::Units, y::Units) = *(x,inv(y))
//(x::Units, y::Units)  = x/y

# Both methods needed for ambiguity resolution
^(x::Unit{U,D}, y::Integer) where {U,D} = Unit{U,D}(tens(x), power(x)*y)
^(x::Unit{U,D}, y::Number) where {U,D} = Unit{U,D}(tens(x), power(x)*y)

# A word of caution:
# Exponentiation is not type-stable for `Units` objects.
# Dimensions get reconstructed anyway so we pass () for the D type parameter...
^(x::FreeUnits{N}, y::Integer) where {N} = *(FreeUnits{map(a->a^y, N), ()}())
^(x::FreeUnits{N}, y::Number) where {N} = *(FreeUnits{map(a->a^y, N), ()}())

^(x::ContextUnits{N,D,P}, y::Integer) where {N,D,P} =
    *(ContextUnits{map(a->a^y, N), (), typeof(P()^y)}())
^(x::ContextUnits{N,D,P}, y::Number) where {N,D,P} =
    *(ContextUnits{map(a->a^y, N), (), typeof(P()^y)}())

^(x::FixedUnits{N}, y::Integer) where {N} = *(FixedUnits{map(a->a^y, N), ()}())
^(x::FixedUnits{N}, y::Number) where {N} = *(FixedUnits{map(a->a^y, N), ()}())

@generated function Base.literal_pow(::typeof(^), x::FreeUnits{N}, ::Type{Val{p}}) where {N,p}
    y = *(FreeUnits{map(a->a^p, N), ()}())
    :($y)
end
@generated function Base.literal_pow(::typeof(^), x::ContextUnits{N,D,P}, ::Type{Val{p}}) where {N,D,P,p}
    y = *(ContextUnits{map(a->a^p, N), (), typeof(P()^p)}())
    :($y)
end
@generated function Base.literal_pow(::typeof(^), x::FixedUnits{N}, ::Type{Val{p}}) where {N,p}
    y = *(FixedUnits{map(a->a^p, N), ()}())
    :($y)
end

# Since exponentiation is not type stable, we define a special `inv` method to enable fast
# division. For julia 0.6.0, the appropriate methods for ^ and * need to be defined before
# this one!
for (fun,pow) in ((:inv, -1//1), (:sqrt, 1//2), (:cbrt, 1//3))
    # The following are generated functions to ensure type stability.
    @eval @generated function ($fun)(x::FreeUnits)
        unittuple = map(x->x^($pow), x.parameters[1])
        y = *(FreeUnits{unittuple,()}())    # sort appropriately
        :($y)
    end

    @eval @generated function ($fun)(x::ContextUnits)
        unittuple = map(x->x^($pow), x.parameters[1])
        promounit = ($fun)(x.parameters[3]())
        y = *(ContextUnits{unittuple,(),typeof(promounit)}())   # sort appropriately
        :($y)
    end

    @eval @generated function ($fun)(x::FixedUnits)
        unittuple = map(x->x^($pow), x.parameters[1])
        y = *(FixedUnits{unittuple,()}())   # sort appropriately
        :($y)
    end
end

function tensfactor(x::Unit)
    p = power(x)
    if isinteger(p)
        p = Integer(p)
    end
    tens(x)*p
end

@generated function tensfactor(x::Units)
    tunits = x.parameters[1]
    a = mapreduce(tensfactor, +, 0, tunits)
    :($a)
end

# This is type unstable but
# a) this method is not called by the user
# b) ultimately the instability will only be present at compile time as it is
# hidden behind a "generated function barrier"
function basefactor(inex, ex, eq, tens, p)
    # Sometimes (x::Rational)^1 can fail for large rationals because the result
    # is of type x*x so we do a hack here
    function dpow(x,p)
        if p == 0
            1
        elseif p == 1
            x
        elseif p == -1
            1//x
        else
            x^p
        end
    end

    if isinteger(p)
        p = Integer(p)
    end

    eq_is_exact = false
    output_ex_float = (10.0^tens * float(ex))^p
    eq_raised = float(eq)^p
    if isa(eq, Integer) || isa(eq, Rational)
        output_ex_float *= eq_raised
        eq_is_exact = true
    end

    can_exact = (output_ex_float < typemax(Int))
    can_exact &= (1/output_ex_float < typemax(Int))
    can_exact &= isinteger(p)

    can_exact2 = (eq_raised < typemax(Int))
    can_exact2 &= (1/eq_raised < typemax(Int))
    can_exact2 &= isinteger(p)

    if can_exact
        if eq_is_exact
            # If we got here then p is an integer.
            # Note that sometimes x^1 can cause an overflow error if x is large because
            # of how power_by_squaring is implemented for Rationals, so we use dpow.
            x = dpow(eq*ex*(10//1)^tens, p)
            return (inex^p, isinteger(x) ? Int(x) : x)
        else
            x = dpow(ex*(10//1)^tens, p)
            return ((inex*eq)^p, isinteger(x) ? Int(x) : x)
        end
    else
        if eq_is_exact && can_exact2
            x = dpow(eq,p)
            return ((inex * ex * 10.0^tens)^p, isinteger(x) ? Int(x) : x)
        else
            return ((inex * ex * 10.0^tens * eq)^p, 1)
        end
    end
end

@inline basefactor(x::Unit{U}) where {U} = basefactor(basefactors[U]..., 1, 0, power(x))

function basefactor(x::Units{U}) where {U}
    fact1 = map(basefactor, U)
    inex1 = mapreduce(x->getfield(x,1), *, 1.0, fact1)
    float_ex1 = mapreduce(x->float(getfield(x,2)), *, 1, fact1)
    can_exact = (float_ex1 < typemax(Int))
    can_exact &= (1/float_ex1 < typemax(Int))
    if can_exact
        inex1, mapreduce(x->getfield(x,2), *, 1, fact1)
    else
        inex1*float_ex1, 1
    end
end
