using StaticArrays

struct Domain{S, T, G, N, L} <: AbstractArray{T, N}
    data::MArray{S, T, N, L}

    function Domain{T}(size::Int...; stencil_size::Int=1) where T
        all(size .> 0) || throw(ArgumentError("size needs to be positive!"))
        stencil_size > 0 || throw(ArgumentError("stencil_size needs to be positive!"))
        stencil_size <= minimum(size) || throw(ArgumentError("stencil_size needs to be smaller than the data size!"))

        G = stencil_size >> 1
        size = size .+ 2G
        S = Tuple{size...}
        N = length(size)
        L = prod(size)

        data = zeros(MArray{S, T, N, L})
        return new{S, T, G, N, L}(data)
    end
end

Base.size(::Domain{S, T, G}) where {S, T, G} = tuple(S.parameters...) .- 2G
Base.IndexStyle(::Type{<:Domain}) = IndexCartesian()
Base.@propagate_inbounds function Base.getindex(domain::Domain{S, T, G, N}, I::Vararg{Int, N}) where {S, T, G, N}
    return domain.data[(I .+ G)...]
end
Base.@propagate_inbounds function Base.setindex!(domain::Domain{S, T, G, N}, value, I::Vararg{Int, N}) where {S, T, G, N}
    return domain.data[(I .+ G)...] = value
end

function Base.similar(domain::Domain{S, T, G}, ::Type{X}, dims::Dims) where {S, T, G, X}
    return Domain{X}(dims...; stencil_size = 2G + 1)
end

function Base.view(domain::Domain{S, T, G}) where {S, T, G}
    range(n) = (1 + G):(n + G)
    @view domain.data[range.(size(domain))...]
end

@inline _boundary_pre(s, g) = 1:(s + 2g)
@inline _boundary_post(s, g) = (1 + g):(s + g)
@inline function _boundary_pivot(s, g)
    a = 1:g
    b = a .+ g
    d = b .+ s
    c = d .- g
    return (a, b, c, d)
end

struct Boundary{A, B}
    pre::NTuple{A, UnitRange{Int}}
    post::NTuple{B, UnitRange{Int}}
    pivot::NTuple{4, UnitRange{Int}}
end
@inline Boundary(s, g, i) = Boundary(
    _boundary_pre.(s[1:(i - 1)], g),
    _boundary_post.(s[(i + 1):length(s)], g),
    _boundary_pivot(s[i], g)
)
@inline Boundary(domain::Domain, i) = Boundary(size(domain), stencil_radius(domain), i)
@inline Base.getindex(I::Boundary, i::Int) = (I.pre..., I.pivot[i], I.post...)
@inline Base.getindex(I::Boundary, i::Int, di::Int) = (I.pre..., I.pivot[i] .+ di, I.post...)
@inline Base.getindex(I::Boundary, ::typeof(!), i::Int) = (I.pre..., reverse(I.pivot[i]), I.post...)
@inline Base.getindex(I::Boundary, ::typeof(!), i::Int, di::Int) = (I.pre..., reverse(I.pivot[i] .+ di), I.post...)

function fixed!(domain::Domain, i::Int, x)
    I = Boundary(domain, i)

    @inbounds begin
        domain.data[I[1]...] = x
        domain.data[I[4]...] = x
    end
end
function fixed!(domain::Domain, x)
    for i in 1:ndims(domain)
        fixed!(domain, i)
    end
end

function periodic!(domain::Domain, i::Int)
    I = Boundary(domain, i)

    @inbounds begin
        domain.data[I[1]...] = domain.data[I[3]...]
        domain.data[I[4]...] = domain.data[I[2]...]
    end

    return nothing
end
function periodic!(domain::Domain)
    for i in 1:ndims(domain)
        periodic!(domain, i)
    end
end

function mirror!(domain::Domain, i::Int)
    I = Boundary(size(domain), stencil_radius(domain), i)

    @inbounds begin
        domain.data[I[1]...] = domain.data[I[!, 2, 1]...]
        domain.data[I[4]...] = domain.data[I[!, 3, -1]...]
    end
end
function mirror!(domain::Domain)
    for i in 1:ndims(domain)
        mirror!(domain, i)
    end
end

abstract type BoundaryType end
function bound!(domain::Domain, type::BoundaryType)
    for i in 1:ndims(domain)
        bound!(domain, type, i)
    end
end
function bound!(domain::Domain{S, T, G, N}, types::Vararg{<:BoundaryType, N}) where {S, T, G, N}
    for (i, type) in enumerate(types)
        bound!(domain, type, i)
    end
end

struct Fixed{T} <: BoundaryType
    value::T
end
function bound!(domain::Domain, fixed::Fixed, i::Int)
    fixed!(domain, i, fixed.value)
end

struct Periodic <: BoundaryType end
function bound!(domain::Domain, ::Periodic, i::Int)
    periodic!(domain, i)
end

struct Mirror <: BoundaryType end
function bound!(domain::Domain, ::Mirror, i::Int)
    mirror!(domain, i)
end

function update!(f, domain::Domain)
    update!(f, domain, domain)
end

function update!(f, input::Domain{S, T, G, N}, output::Domain{S, T, G, N}) where {S, T, G, N}
    dI = CartesianIndex(Iterators.repeated(G, N)...)

    for I in keys(output)
        stencil = @view input.data[I:(I + 2dI)]
        output.data[I + dI] = f(stencil)
    end

    return output
end

stencil_radius(domain::Domain{S, T, G}) where {S, T, G} = G
stencil_size(domain::Domain{S, T, G}) where {S, T, G} = 2G + 1