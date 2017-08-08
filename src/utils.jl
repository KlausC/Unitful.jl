@inline isunitless(::Units) = false
@inline isunitless(::Units{()}) = true

@inline numtype(::Quantity{T}) where {T} = T
@inline numtype(::Type{Quantity{T,D,U}}) where {T,D,U} = T
@inline dimtype(u::Unit{U,D}) where {U,D} = D

"""
    ustrip(x::Number)
Returns the number out in front of any units. This may be different from the value
in the case of dimensionless quantities. See [`uconvert`](@ref) and the example
below. Because the units are removed, information may be lost and this should
be used with some care.

This function is mainly intended for compatibility with packages that don't know
how to handle quantities. This function may be deprecated in the future.

```jldoctest
julia> ustrip(2u"μm/m") == 2
true

julia> uconvert(NoUnits, 2u"μm/m") == 2//1000000
true
```
"""
@inline ustrip(x::Number) = x/unit(x)

"""
    ustrip(x::Array{Q}) where {Q <: Quantity}
Strip units from an `Array` by reinterpreting to type `T`. The resulting
`Array` is a "unit free view" into array `x`. Because the units are
removed, information may be lost and this should be used with some care.

This function is provided primarily for compatibility purposes; you could pass
the result to PyPlot, for example. This function may be deprecated in the future.

```jldoctest
julia> a = [1u"m", 2u"m"]
2-element Array{Quantity{Int64, Dimensions:{𝐋}, Units:{m}},1}:
 1 m
 2 m

julia> b = ustrip(a)
2-element Array{Int64,1}:
 1
 2

julia> a[1] = 3u"m"; b
2-element Array{Int64,1}:
 3
 2
```
"""
@inline ustrip(x::Array{Q}) where {Q <: Quantity} = reinterpret(numtype(Q), x)

"""
    ustrip(A::AbstractArray{Q}) where {Q <: Quantity}
Strip units from an `AbstractArray` by making a new array without units using
array comprehensions.

This function is provided primarily for compatibility purposes; you could pass
the result to PyPlot, for example. This function may be deprecated in the future.
"""
ustrip(A::AbstractArray{Q}) where {Q <: Quantity} = (numtype(Q))[ustrip(x) for x in A]

"""
    ustrip(x::AbstractArray{T}) where {T <: Number}
Fall-back that returns `x`.
"""
@inline ustrip(A::AbstractArray{T}) where {T <: Number} = A


ustrip(A::Diagonal{T}) where {T <: Quantity} = Diagonal(ustrip(A.diag))
ustrip(A::Bidiagonal{T}) where {T <: Quantity} =
    Bidiagonal(ustrip(A.dv), ustrip(A.ev), A.isupper)
ustrip(A::Tridiagonal{T}) where {T <: Quantity} =
    Tridiagonal(ustrip(A.dl), ustrip(A.d), ustrip(A.du))
ustrip(A::SymTridiagonal{T}) where {T <: Quantity} =
    SymTridiagonal(ustrip(A.dv), ustrip(A.ev))

"""
    unit{T,D,U}(x::Quantity{T,D,U})
Returns the units associated with a quantity.

Examples:

```jldoctest
julia> unit(1.0u"m") == u"m"
true

julia> typeof(u"m")
Unitful.FreeUnits{(Unitful.Unit{:Meter,Unitful.Dimensions{(Unitful.Dimension{:Length}(1//1),)}}(0, 1//1),),Unitful.Dimensions{(Unitful.Dimension{:Length}(1//1),)}}
```
"""
@inline unit(x::Quantity{T,D,U}) where {T,D,U} = U()

"""
    unit{T,D,U}(x::Type{Quantity{T,D,U}})
Returns the units associated with a quantity type, `ContextUnits(U(),P())`.

Examples:

```jldoctest
julia> unit(typeof(1.0u"m")) == u"m"
true
```
"""
@inline unit(::Type{Quantity{T,D,U}}) where {T,D,U} = U()


"""
    unit(x::Number)
Returns a `Unitful.Units{(), Dimensions{()}}` object to indicate that ordinary
numbers have no units. This is a singleton, which we export as `NoUnits`.
The unit is displayed as an empty string.

Examples:

```jldoctest
julia> typeof(unit(1.0))
Unitful.FreeUnits{(),Unitful.Dimensions{()}}
julia> typeof(unit(Float64))
Unitful.FreeUnits{(),Unitful.Dimensions{()}}
julia> unit(1.0) == NoUnits
true
```
"""
@inline unit(x::Number) = NoUnits
@inline unit(x::Type{T}) where {T <: Number} = NoUnits

"""
    dimension(x::Number)
    dimension{T<:Number}(x::Type{T})
Returns a `Unitful.Dimensions{()}` object to indicate that ordinary
numbers are dimensionless. This is a singleton, which we export as `NoDims`.
The dimension is displayed as an empty string.

Examples:

```jldoctest
julia> typeof(dimension(1.0))
Unitful.Dimensions{()}
julia> typeof(dimension(Float64))
Unitful.Dimensions{()}
julia> dimension(1.0) == NoDims
true
```
"""
@inline dimension(x::Number) = NoDims
@inline dimension(x::Type{T}) where {T <: Number} = NoDims

"""
    dimension{U,D}(u::Units{U,D})
Returns a [`Unitful.Dimensions`](@ref) object corresponding to the dimensions
of the units, `D()`. For a dimensionless combination of units, a
`Unitful.Dimensions{()}` object is returned.

Examples:

```jldoctest
julia> dimension(u"m")
𝐋

julia> typeof(dimension(u"m"))
Unitful.Dimensions{(Unitful.Dimension{:Length}(1//1),)}

julia> typeof(dimension(u"m/km"))
Unitful.Dimensions{()}
```
"""
@inline dimension(u::Units{U,D}) where {U,D} = D()

"""
    dimension{T,D}(x::Quantity{T,D})
Returns a [`Unitful.Dimensions`](@ref) object `D()` corresponding to the
dimensions of quantity `x`. For a dimensionless [`Unitful.Quantity`](@ref), a
`Unitful.Dimensions{()}` object is returned.

Examples:

```jldoctest
julia> dimension(1.0u"m")
𝐋

julia> typeof(dimension(1.0u"m/μm"))
Unitful.Dimensions{()}
```
"""
@inline dimension(x::Quantity{T,D}) where {T,D} = D()*dimension(T)
@inline dimension(::Type{Quantity{T,D,U}}) where {T,D,U} = D()*dimension(T)

@deprecate(dimension(x::AbstractArray{T}) where {T<:Number}, dimension.(x))
@deprecate(dimension(x::AbstractArray{T}) where {T<:Units}, dimension.(x))

"""
    mutable struct DimensionError{T,S} <: Exception
      x::T
      y::S
    end
Thrown when dimensions don't match in an operation that demands they do.
Display `x` and `y` in error message.
"""
mutable struct DimensionError{T,S} <: Exception
    x::T
    y::S
end
Base.showerror(io::IO, e::DimensionError) =
    print(io, "DimensionError: $(e.x) and $(e.y) are not dimensionally compatible.");
