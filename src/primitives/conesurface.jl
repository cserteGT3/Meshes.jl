# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ConeSurface(base, apex)

A cone surface with `base` disk and `apex`.
See https://en.wikipedia.org/wiki/Cone.

See also [`Cone`](@ref).
"""
struct ConeSurface{T} <: Primitive{3,T}
  base::Disk{T}
  apex::Point{3,T}
end

ConeSurface(base::Disk, apex::Tuple) = ConeSurface(base, Point(apex))

paramdim(::Type{<:ConeSurface}) = 2

base(c::ConeSurface) = c.base

apex(c::ConeSurface) = c.apex

height(c::ConeSurface) = norm(center(base(c)) - apex(c))

halfangle(c::ConeSurface) = atan(radius(base(c)), height(c))

axis(c::ConeSurface) = Line(apex(c), center(base(c)))

# around the axis, height
function (c::ConeSurface{T})(φ, z) where {T}
  if (φ < 0 || φ > 1) || (z < 0 || z > 1)
    throw(DomainError((φ, z), "c(φ, z) is not defined for φ, z outside [0, 1]²."))
  end
  a = axis(c)
  d = a(1) - a(0)
  l = norm(d)
  rb = radius(base(c))
  
  # rotation to align z axis with cylinder axis
  Q = rotation_between(d, Vec{3,T}(0, 0, 1))

  # scale coordinates
  φₛ = 2T(π) * φ
  zₛ = z * l

  x_local = cos(φₛ) * (zₛ * rb) / l
  y_local = sin(φₛ) * (zₛ * rb) / l
  z_local = zₛ
  p_local = Vec{3,T}(x_local, y_local, z_local)

  apex(c) + Q' * p_local
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{ConeSurface{T}}) where {T} =
  ConeSurface(rand(rng, Disk{T}), rand(rng, Point{3,T}))
