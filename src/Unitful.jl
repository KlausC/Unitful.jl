__precompile__(true)
module Unitful

import Base: ==, <, <=, +, -, *, /, //, ^
import Base: show, convert
import Base: abs, abs2, float, fma, muladd, inv, sqrt, cbrt
import Base: min, max, floor, ceil, real, imag, conj
import Base: exp, exp10, exp2, expm1, log, log10, log1p, log2
import Base: sin, cos, tan, cot, sec, csc, atan2, cis, vecnorm

import Base: mod, rem, div, fld, cld, trunc, round, sign, signbit
import Base: isless, isapprox, isinteger, isreal, isinf, isfinite, isnan
import Base: copysign, flipsign
import Base: prevfloat, nextfloat, maxintfloat, rat, step #, linspace
import Base: length, float, start, done, next, last, one, zero, colon#, range
import Base: getindex, eltype, step, last, first, frexp
import Base: Integer, Rational, typemin, typemax
import Base: steprange_last, unsigned

import Base.LinAlg: istril, istriu

export unit, dimension, uconvert, ustrip, upreferred
export @dimension, @derived_dimension, @refunit, @unit, @u_str
export Quantity
export DimensionlessQuantity
export NoUnits, NoDims

const unitmodules = Vector{Module}()
const basefactors = Dict{Symbol,Tuple{Float64,Rational{Int}}}()

include("types.jl")
const promotion = Dict{Symbol,Unit}()

include("user.jl")
include("utils.jl")
include("dimensions.jl")
include("units.jl")
include("quantities.jl")
include("display.jl")
include("promotion.jl")
include("conversion.jl")
include("range.jl")
include("fastmath.jl")
include("pkgdefaults.jl")
include("level.jl")

function __init__()
    # @u_str should be aware of units defined in module Unitful
    Unitful.register(Unitful)
end

end
