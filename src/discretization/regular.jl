# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    RegularDiscretization(n1, n2, ..., np)

A method to discretize primitive geometries with
`n1×n2×...×np` elements sampled regularly along
each parametric dimensions. The adequate number
of points is calculated for each type of geometry
and passed to [`RegularSampling`](@ref).
"""
struct RegularDiscretization{N} <: DiscretizationMethod
  sizes::Dims{N}
end

RegularDiscretization(sizes::Vararg{Int,N}) where {N} = RegularDiscretization(sizes)

function discretize(geometry::Geometry, method::RegularDiscretization)
  if isparametrized(geometry)
    sz = fitdims(method.sizes, paramdim(geometry))
    ip = isperiodic(geometry)
    np = @. sz + !ip

    points = sample(geometry, RegularSampling(np))
    topo = GridTopology(sz, ip)

    SimpleMesh(collect(points), topo)
  else
    box = boundingbox(geometry)
    grid = discretize(box, method)
    view(grid, geometry)
  end
end

# --------------
# SPECIAL CASES
# --------------

function discretize(box::Box, method::RegularDiscretization)
  sz = fitdims(method.sizes, paramdim(box))
  CartesianGrid(extrema(box)..., dims=sz)
end

function discretize(sphere::Sphere{3}, method::RegularDiscretization)
  nx, ny = fitdims(method.sizes, paramdim(sphere))

  # sample points regularly
  sampler = RegularSampling(nx, ny)
  points = collect(sample(sphere, sampler))

  # connect regular samples with quadrangles
  topo = GridTopology((nx - 1, ny - 1))
  middle = collect(elements(topo))
  offset = nx * ny - nx
  for i in 1:(nx - 1)
    u = offset + i
    v = offset + i + 1
    w = i + 1
    z = i
    quad = connect((u, v, w, z))
    push!(middle, quad)
  end

  # connect north pole with triangles
  north = map(1:(ny - 1)) do j
    u = nx * ny + 1
    v = 1 + (j - 1) * nx
    w = 1 + (j) * nx
    connect((u, v, w))
  end
  u = nx * ny + 1
  v = 1 + (ny - 1) * nx
  w = 1
  push!(north, connect((u, v, w)))

  # connect south pole with triangles
  south = map(1:(ny - 1)) do j
    u = nx * ny + 2
    v = (j) * nx
    w = (j + 1) * nx
    connect((u, w, v))
  end
  u = nx * ny + 2
  v = ny * nx
  w = nx
  push!(south, connect((u, w, v)))

  connec = [middle; north; south]

  SimpleMesh(points, connec)
end

discretize(ball::Ball{2}, method::RegularDiscretization) = _rball(ball, method)

discretize(disk::Disk, method::RegularDiscretization) = _rball(disk, method)

function _rball(ball, method::RegularDiscretization)
  nx, ny = fitdims(method.sizes, paramdim(ball))

  # sample points regularly
  sampler = RegularSampling(nx, ny)
  points = collect(sample(ball, sampler))

  # connect regular samples with quadrangles
  topo = GridTopology((nx - 1, ny), (false, true))
  quads = collect(elements(topo))

  # connect center with triangles
  tris = map(1:(ny - 1)) do j
    u = nx * ny + 1
    v = 1 + (j - 1) * nx
    w = 1 + (j) * nx
    connect((u, v, w))
  end
  u = nx * ny + 1
  v = 1 + (ny - 1) * nx
  w = 1
  push!(tris, connect((u, v, w)))

  connec = [quads; tris]

  SimpleMesh(points, connec)
end

function discretize(cylsurf::CylinderSurface, method::RegularDiscretization)
  nx, ny = fitdims(method.sizes, paramdim(cylsurf))

  # sample points regularly
  sampler = RegularSampling(nx, ny)
  points = collect(sample(cylsurf, sampler))

  # connect regular samples with quadrangles
  topo = GridTopology((nx - 1, ny - 1))
  middle = collect(elements(topo))
  for j in 1:(ny - 1)
    u = (j) * nx
    v = (j - 1) * nx + 1
    w = (j) * nx + 1
    z = (j + 1) * nx
    quad = connect((u, v, w, z))
    push!(middle, quad)
  end

  # connect south pole with triangles
  south = map(1:(nx - 1)) do i
    u = nx * ny + 1
    v = i + 1
    w = i
    connect((u, v, w))
  end
  u = nx * ny + 1
  v = 1
  w = nx
  push!(south, connect((u, v, w)))

  # connect north pole with triangles
  offset = nx * ny - nx
  north = map(1:(nx - 1)) do i
    u = nx * ny + 2
    v = offset + i + 1
    w = offset + i
    connect((u, w, v))
  end
  u = nx * ny + 2
  v = nx * ny - nx + 1
  w = nx * ny
  push!(north, connect((u, w, v)))

  connec = [middle; south; north]

  SimpleMesh(points, connec)
end

function discretize(frustsurf::FrustumSurface, method::RegularDiscretization)
  nx, ny = fitdims(method.sizes, paramdim(frustsurf))

  # sample points regularly
  sampler = RegularSampling(nx, ny)
  points = collect(sample(frustsurf, sampler))

  # connect regular samples with quadrangles
  topo = GridTopology((nx - 1, ny - 1))
  middle = collect(elements(topo))
  for j in 1:(ny - 1)
    u = (j) * nx
    v = (j - 1) * nx + 1
    w = (j) * nx + 1
    z = (j + 1) * nx
    quad = connect((u, v, w, z))
    push!(middle, quad)
  end

  # connect south pole with triangles
  south = map(1:(nx - 1)) do i
    u = nx * ny + 1
    v = i + 1
    w = i
    connect((u, v, w))
  end
  u = nx * ny + 1
  v = 1
  w = nx
  push!(south, connect((u, v, w)))

  # connect north pole with triangles
  offset = nx * ny - nx
  north = map(1:(nx - 1)) do i
    u = nx * ny + 2
    v = offset + i + 1
    w = offset + i
    connect((u, w, v))
  end
  u = nx * ny + 2
  v = nx * ny - nx + 1
  w = nx * ny
  push!(north, connect((u, w, v)))

  # currently: only discretize to conical surface and skip bottom-top planes
  #connec = [middle; south; north]
  connec = middle

  SimpleMesh(points, connec)
end
