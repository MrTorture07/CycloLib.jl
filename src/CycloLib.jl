module CycloLib
    
include("gear.jl")
include("materials.jl")
include("forcecalc.jl")
include("sector.jl")

export CircularGearing
export DefaultCircularGearing
export cycloid_point, cycloid_points
export cycloid_normal, cycloid_normals
export satellite_point, satellite_points
export cycloid_curvature, satellite_curvature, unit_curvature

export SectorGearing
export excenter_coords, excenter_angle

end
