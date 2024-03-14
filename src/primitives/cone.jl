# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Cone(base, apex)

A cone with `base` disk and `apex`.
See <https://en.wikipedia.org/wiki/Cone>.

See also [`ConeSurface`](@ref).
"""
struct Cone{T} <: Primitive{3,T}
  base::Disk{T}
  apex::Point{3,T}
end

Cone(base::Disk, apex::Tuple) = Cone(base, Point(apex))

paramdim(::Type{<:Cone}) = 3

base(c::Cone) = c.base

apex(c::Cone) = c.apex

height(c::Cone) = norm(center(base(c)) - apex(c))

halfangle(c::Cone) = atan(radius(base(c)), height(c))

axis(c::Cone) = Line(apex(c), center(base(c)))

# around the axis, opening angle, height
function (c::Cone{T})(φ, Ψ, z) where {T}
  if (φ < 0 || φ > 1) || (z < 0 || z > 1) || (Ψ < 0 || Ψ > 1)
    throw(DomainError((φ, Ψ, z), "f(φ, Ψ, z) is not defined for φ, Ψ, z outside [0, 1]³."))
  end
  a = axis(c)
  d = a(1) - a(0)
  l = norm(d)
  
  # rotation to align z axis with cylinder axis
  Q = rotation_between(d, Vec{3,T}(0, 0, 1))

  # scale coordinates
  φₛ = 2T(π) * φ
  zₛ = z * l
  Ψₛ = halfangle(c) * Ψ

  x_local = cos(φₛ) * (zₛ * tan(Ψₛ)) / l
  y_local = sin(φₛ) * (zₛ * tan(Ψₛ)) / l
  z_local = zₛ
  p_local = Vec{3,T}(x_local, y_local, z_local)

  apex(c) + Q' * p_local
end

Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Cone{T}}) where {T} =
  Cone(rand(rng, Disk{T}), rand(rng, Point{3,T}))
