export FixedReferenceLevel,
    ArbitraryReferenceLevel

"""
    @logunit(symb,abbr,name,base,ref,pre)
Define a log scale-based unit. Rather than specifying a dimension like in
[`@refunit`](@ref), `ref` should be a number equal to zero of the
log unit being defined. The `base` of the logarithm should be specified.

Returns the [`Unitful.FreeUnits`](@ref) object to which `symb` is bound.

Usage example: `@logunit dBm "dBm" dBm 10 1mW 10`
"""
macro logunit(symb,abbr,name,base,ref,pre)
    x = Expr(:quote, name)
    quote
        if $base isa Number
            d = Unitful.dimension($(esc(ref)))
            Unitful.abbr(::Unitful.LogUnit{$(esc(x)), typeof(d), $base, $ref, $pre}) = $abbr
            Unitful.@log_unit_symbols($(esc(symb)), $(esc(name)), d, $base, $ref, $pre)
            $(esc(symb))
        else
            error("logarithm base must be a number.")
        end
    end
end

"""
    @log_unit_symbols(symb,name,dimension,base,ref,pre)
Not called directly by the user. Given a unit symbol and a unit's name,
will define units without SI power-of-ten prefixes.

Example: `@unit_symbols ft Foot 𝐋` results in `ft` getting defined but not `kft`.
"""
macro log_unit_symbols(symb,name,dimension,base,ref,pre)
    s = Symbol(symb)
    z = Expr(:quote, name)
    u = :(Unitful.LogUnit{$z, typeof($dimension), $base, $ref, $pre}())
    esc(quote
        const $s = Unitful.FreeUnits{($u,), typeof(Unitful.dimension($u))}()
    end)
end

"""
    struct LogInfo{N,B,P}
Describes a logarithmic unit. Type parameters include:
- `N`: The name of the logarithmic unit, e.g. `:Bel`, `:Neper`.
- `B`: The base of the logarithm
- `P`: A prefactor to multiply the logarithm
"""
struct LogInfo{N,B,P} end
abbr(::LogInfo{:Bel})     = "B"
abbr(::LogInfo{:Decibel}) = "dB"
abbr(::LogInfo{:Neper})   = "Np"
base(::LogInfo{N,B}) where {N,B} = B
prefactor(::LogInfo{N,B,P}) where {N,B,P} = P

abstract type Level{L,S,T<:Number} <: Number end
dimension(x::Level) = dimension(reflevel(x))
dimension(x::Type{T}) where {L,S,T<:Level{L,S}} = dimension(S)
function abbr(x::Level{L,S}) where {L,S}
    if dimension(S) == NoDims
        return abbr(L())
    else
        return join([abbr(L()), " (", reflevel(x), ")"])
    end
end

function uconvert(a::Units, x::Level)
    if dimension(a) == dimension(x)
        uconvert(a, x.val)
    else
        throw(DimensionError(a,x))
    end
end

"""
    struct ArbitraryReferenceLevel{L, S<:Number, T<:Number} <: Number
A logarithmic scale-based level with arbitrary reference. Details about the logarithmic
scale are encoded in `L`, which is an object of type `T <: LogInfo`.

The type parameter `D` contains dimension information, e.g. `typeof(NoUnits)` or
`typeof(𝐌*𝐋^2/𝐓^3)`. Powers of a log-based unit are not allowed. `R` is a reference
quantity, as done e.g. with the `dBm` unit, which is referred to 1 mW.
"""
struct ArbitraryReferenceLevel{L<:LogInfo, S<:Number, T<:Number} <: Level{L,S,T}
    ref::S
    val::T
    function ArbitraryReferenceLevel{L,S,T}(x,y) where {L,S,T}
        dimension(x) != dimension(y) && throw(DimensionError(x,y))
        return new{L,S,T}(x,y)
    end
end
function ArbitraryReferenceLevel{L}(ref::Number, val::Number) where L <: LogInfo
    dimension(ref) != dimension(val) && throw(DimensionError(ref, val))
    return ArbitraryReferenceLevel{L, typeof(ref), typeof(val)}(ref, val)
end
Base.convert(T::Type{<:ArbitraryReferenceLevel}, x::Level) = T(reflevel(x), x.val)
reflevel(x::ArbitraryReferenceLevel) = x.ref

"""
    struct FixedReferenceLevel{L<:LogInfo, S, T<:Number} <: Number
A logarithmic scale-based level with arbitrary reference. Details about the logarithmic
scale are encoded in `L`, which is a type `<: LogInfo`.
"""
struct FixedReferenceLevel{L<:LogInfo, S, T<:Number} <: Level{L,S,T}
    val::T
    function FixedReferenceLevel{L,S,T}(x) where {L,S,T}
        dimension(S) != dimension(x) && throw(DimensionError(S,x))
        return new{L,S,T}(x)
    end
    function FixedReferenceLevel{L,S,T}(r,x) where {L,S,T}
        r != S && throw(ArgumentError("invalid FixedReferenceLevel constructor."))
        return new{L,S,T}(x)
    end
end
function FixedReferenceLevel{L,S}(val::Number) where {L,S}
    dimension(S) != dimension(val) && throw(DimensionError(S, val))
    return FixedReferenceLevel{L,S,typeof(val)}(val)
end
function FixedReferenceLevel{L}(ref::Number, val::Number) where {L}
    dimension(ref) != dimension(val) && throw(DimensionError(ref,val))
    return FixedReferenceLevel{L,ref,typeof(val)}(val)
end
Base.convert(T::Type{<:FixedReferenceLevel}, x::Level) = T(x.val)
reflevel(x::FixedReferenceLevel{L,S}) where {L,S} = S

Level{L,S}(ref, val) where {L,S<:Number} = ArbitraryReferenceLevel{L,S,typeof(val)}(ref, val)
Level{L,S}(ref, val) where {L,S} = FixedReferenceLevel{L,S,typeof(val)}(ref, val)

for (_short,_long,_base,_pre) in (  (:B, :Bel, 10, 1),
                                (:dB, :Decibel, 10, 10),
                                (:Np, :Neper, e, 1//2))
    li = Symbol("li_",_short)
    @eval begin
        const $li = LogInfo{$(QuoteNode(_long)),$_base,$_pre}
        const $_short = FixedReferenceLevel{$li, 1}
    end
end
*(x::Real, y::Type{T}) where {L<:LogInfo, S, T<:FixedReferenceLevel{L,S}} =
    (T)(S*expfn(L())(x/((1+isrootpower(L,dimension(S)))*prefactor(L()))))
*(x::Type{T}, y::Real) where {L<:LogInfo, S, T<:FixedReferenceLevel{L,S}} = *(y,x)

const dBV  = FixedReferenceLevel{li_dB, 1V}
const dBu  = FixedReferenceLevel{li_dB, sqrt(0.6)V}
const dBμV = FixedReferenceLevel{li_dB, 1μV}
const dBµV = dBμV # different encoding of mu
const dBm  = FixedReferenceLevel{li_dB, 1mW}

ustrip(x::Level{T,S}) where {T<:LogInfo, S} =
    (1+isrootpower(T,dimension(S)))*prefactor(T())*(logfn(T()))(x.val/reflevel(x))

# TODO: some more dimensions?
isrootpower(x,y) = isrootpower_warn(x,y)
isrootpower(::Type{<:LogInfo}, ::typeof(dimension(W))) = false
isrootpower(::Type{<:LogInfo}, ::typeof(dimension(V))) = true
isrootpower(::Type{<:LogInfo}, ::typeof(dimension(A))) = true

# Default to power or root-power as appropriate for the given logarithmic unit
function isrootpower_warn(x,y)
    irp = isrootpower(x)
    str = ifelse(irp, "root-power", "power")
    warn("assuming ratios of quantities of dimension ", y, " are ", str, " ratios. Define ",
        "`Unitful.isrootpower(::Type{<:Unitful.LogInfo}, ::typeof($y))` to fix.")
    return irp
end

isrootpower(t::Type{<:LogInfo}, ::typeof(NoDims)) = isrootpower(t)
isrootpower(::Type{li_B}) = false
isrootpower(::Type{li_dB}) = false
isrootpower(::Type{li_Np}) = true

for _T in (:Number, :Quantity) #ambiguity resolution
    # Addition
    @eval +(x::$_T, y::Level) = +(y,x)
    @eval +(x::Level, y::$_T) = +(promote(x,y)...)

    # @eval +(x::FixedReferenceLevel{L,S}, y::$_T) where {L,S} =
    #     FixedReferenceLevel{L,S}(x.val+y)
    # @eval +(x::ArbitraryReferenceLevel{L}, y::$_T) where {L} =
    #     ArbitraryReferenceLevel{L}(reflevel(x), x.val+y)

    # Subtraction
    @eval -(x::$_T, y::FixedReferenceLevel{L,S}) where {L,S} =
        FixedReferenceLevel{L,S}(x-y.val)
    @eval -(x::$_T, y::ArbitraryReferenceLevel{L}) where {L} =
        ArbitraryReferenceLevel{L}(reflevel(y), x-y.val)
    @eval -(x::FixedReferenceLevel{L,S}, y::$_T) where {L,S} =
        FixedReferenceLevel{L,S}(x.val-y)
    @eval -(x::ArbitraryReferenceLevel{L}, y::$_T) where {L} =
        ArbitraryReferenceLevel{L}(reflevel(x), x.val-y)
end

+(x::FixedReferenceLevel{L,S}, y::FixedReferenceLevel{L,S}) where {L,S} =
    FixedReferenceLevel{L,S}(x.val + y.val)
+(x::FixedReferenceLevel{L,S1}, y::FixedReferenceLevel{L,S2}) where {L,S1,S2} =
    +(promote(x,y)...)
+(x::ArbitraryReferenceLevel{L}, y::ArbitraryReferenceLevel{L}) where {L} = x.val+y.val
+(x::Level, y::Level) = +(promote(x,y)...)

# Multiplication
*(y::Number, x::Level) = *(x,y)
*(x::Level, y::Number) = (typeof(x))(reflevel(x), x.val*y)

# Division
/(x::Level{L,S}, y::Number) where {L,S} = Level{L,S}(reflevel(x), x.val / y)
/(x::Number, y::Level) = error("cannot divide a number by a level.")

valtype(x::Type{X}) where {L,S,T,X<:Level{L,S,T}} = T
function (Base.promote_rule(::Type{S}, ::Type{T})
        where {A1,A2,B1,B2, C1<:Number, C2<:Number, S<:Level{A1,B1,C1}, T<:Level{A2,B2,C2}})
    if A1 == A2
        if B1 isa Number && B2 isa Number
            if B1 == B2
                # Use convert(promote_type(typeof(B1), typeof(B2)), B1) instead of B1?
                return FixedReferenceLevel{A1, B1, promote_type(C1,C2)}
            else
                return ArbitraryReferenceLevel{A1, promote_type(typeof(B1), typeof(B2)),
                    promote_type(C1,C2)}
            end
        elseif B1 isa Number
            return ArbitraryReferenceLevel{A1, promote_type(typeof(B1), B2),
                promote_type(C1,C2)}
        elseif B2 isa Number
            return ArbitraryReferenceLevel{A1, promote_type(B1, typeof(B2)),
                promote_type(C1,C2)}
        else
            return ArbitraryReferenceLevel{A1, promote_type(B1,B2), promote_type(C1,C2)}
        end
    else
        return promote_type(C1,C2)
    end
end

function Base.promote_rule(::Type{X}, ::Type{N}) where {X<:Level, N<:Number}
    return promote_type(valtype(X), N)
end
function Base.promote_rule(::Type{N}, ::Type{X}) where {X<:Level, N<:Number}
    return promote_type(valtype(X), N)
end
function Base.promote_rule(::Type{X}, ::Type{Quantity{T,D,U}}) where {X<:Level,T<:Number,D,U}
    return promote_type(valtype(X), Quantity{T,D,U})
end
function Base.promote_rule(::Type{Quantity{T,D,U}}, ::Type{X}) where {X<:Level,T<:Number,D,U}
    return promote_type(valtype(X), Quantity{T,D,U})
end

Base.convert(::Type{Quantity{T,D,U}}, x::Level) where {T,D,U} =
    convert(Quantity{T,D,U}, x.val)
Base.convert(::Type{Quantity{T}}, x::Level) where {T<:Number} =
    convert(Quantity{T}, x.val)
Base.convert(::Type{S}, x::Level) where {S<:Number} = convert(S, x.val)

function Base.show(io::IO, x::FixedReferenceLevel)
    print(io, ustrip(x), " ", abbr(x))
    nothing
end
abbr(::FixedReferenceLevel{li_dB, 1mW}) = "dBm"
abbr(::FixedReferenceLevel{li_dB, 1V}) = "dBV"
abbr(::FixedReferenceLevel{li_dB, sqrt(0.6)V}) = "dBu"
abbr(::FixedReferenceLevel{li_dB, 1μV}) = "dBμV"
abbr(x::FixedReferenceLevel{L,1}) where {L} = abbr(L())

function Base.show(io::IO, x::ArbitraryReferenceLevel)
    print(io, ustrip(x), " ", abbr(x))
    nothing
end

function Base.show(io::IO, x::Quantity{<:Level})
    print(io, "[")
    show(io,x.val)
    print(io, "]")
    if !isunitless(unit(x))
        print(io," ")
        show(io, unit(x))
    end
    nothing
end

for li in (:B, :dB, :Np)
    @eval begin
        $(Expr(:export, Symbol("@",li)))

        macro ($li)(r::Union{Real,Symbol})
            s = $(Symbol("_", li))
            :(($s)($r, 1))
        end

        macro ($li)(expr::Expr)
            s = $(Symbol("_", li))
            expr.args[1] != :/ &&
                throw(ArgumentError(join(["usage: `@", $(String(li)), " (a)/(b)`"])))
            length(expr.args) != 3 &&
                throw(ArgumentError(join(["usage: `@", $(String(li)), " (a)/(b)`"])))
            :(($s)($(expr.args[2]), $(expr.args[3])))
        end

        function $(Symbol("_", li))(num::Number, den::Number)
            dimension(num) != dimension(den) && throw(DimensionError(num,den))
            if den in fixedlevels
                return FixedReferenceLevel{$(Symbol("li_", li)), den}(num)
            else
                return ArbitraryReferenceLevel{$(Symbol("li_", li))}(den, num)
            end
        end

        function $(Symbol("_", li))(num::Number, den::Units)
            $(Symbol("_", li))(num, 1*den)
        end
    end
end

"""
When these are encountered in the denominator of the @dB macro, a FixedReferenceLevel
is used instead of an ArbitraryReferenceLevel.
"""
const fixedlevels = Any[1mW, 1V, sqrt(0.6)V, 1μV, 1]

# # powerratio(x::Decibels) = exp10(ustrip(x)/10)
# # rootpowerratio(x::Decibels) = exp10(ustrip(x)/20)
# # amplituderatio = rootpowerratio
# # fieldratio = rootpowerratio

"""
    logfn(x::LogInfo)
=Returns the appropriate logarithm function to use in calculations involving the
logarithmic unit / quantity. For example, decibel-based units yield `log10`,
Neper-based yield `ln`, and so on. Returns `x->log(base, x)` as a fallback.
"""
function logfn end
logfn(x::LogInfo{N,10}) where {N} = log10
logfn(x::LogInfo{N,2}) where {N} = log2
logfn(x::LogInfo{N,e}) where {N} = log
logfn(x::LogInfo{N,B}) where {N,B} = x->log(B,x)

"""
    expfn(x::LogInfo)
Returns the appropriate exponential function to use in calculations involving the
logarithmic unit / quantity. For example, decibel-based units yield `exp10`,
Neper-based yield `exp`, and so on. Returns `x->(base)^x` as a fallback.
"""
function expfn end
expfn(x::LogInfo{N,10}) where {N} = exp10
expfn(x::LogInfo{N,2}) where {N} = exp2
expfn(x::LogInfo{N,e}) where {N} = exp
expfn(x::LogInfo{N,B}) where {N,B} = x->B^x
