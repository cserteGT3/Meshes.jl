# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Frustum(bot, top)

A frustum (truncated cone) with `bot` and `top` disks.
See <https://en.wikipedia.org/wiki/Frustum>.

See also [`FrustumSurface`](@ref).
"""
struct Frustum{T} <: Primitive{3,T}
  bot::Disk{T}
  top::Disk{T}

  function Frustum{T}(bot, top) where {T}
    bn = normal(plane(bot))
    tn = normal(plane(top))
    @assert bn ⋅ tn ≈ 1 "Bottom and top plane must be parallel"
    @assert center(bot) ≉ center(top) "Bottom and top centers need to be distinct"
    new(bot, top)
  end
end

Frustum(bot::Disk{T}, top::Disk{T}) where {T} = Frustum{T}(bot, top)

paramdim(::Type{<:Frustum}) = 3

bottom(f::Frustum) = f.bot

top(f::Frustum) = f.top

height(f::Frustum) = height(boundary(f))

axis(f::Frustum) = axis(boundary(f))

paramdim(::Type{<:Frustum}) = 3

axis(f::Frustum) = Line(center(bottom(f)), center(top(f)))

# around the axis, along the radius, height
function (f::Frustum{T})(φ, r, z) where {T}
  if (φ < 0 || φ > 1) || (z < 0 || z > 1) || (r < 0 || r > 1)
    throw(DomainError((φ, r, z), "f(φ, r, z) is not defined for φ, r, z outside [0, 1]³."))
  end
  rb = radius(bottom(f))
  rt = radius(top(f))
  a = axis(f)
  d = a(1) - a(0)
  l = norm(d)

  # rotation to align z axis with cylinder axis
  Q = rotation_between(d, Vec{3,T}(0, 0, 1))

  # scale coordinates
  φₛ = 2T(π) * φ
  zₛ = z * l
  rₛ = r * (rb * (l - zₛ) + rt * zₛ)

  x_local = cos(φₛ) * rₛ / l
  y_local = sin(φₛ) * rₛ / l
  z_local = zₛ
  p_local = Vec{3,T}(x_local, y_local, z_local)

  center(bottom(f)) + Q' * p_local
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Frustum{T}}) where {T}
  bottom = rand(rng, Disk{T})
  ax = normal(plane(bottom))
  topplane = Plane{T}(center(bottom) + rand(T) * ax, ax)
  top = Disk{T}(topplane, rand(T))
  Frustum(bottom, top)
end
