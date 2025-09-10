include("materials.jl")
include("gear.jl")

"""Возвращает приведённый модуль упругости"""
function unit_young(E1::Number, E2::Number, μ1::Number, μ2::Number)::Float64
    inv_E = (1-μ1^2)/E1 + (1-μ2^2)/E2
    return 1 / inv_E
end

"""Возвращает приведённый модуль упругости"""
function unit_young(mat1::Material, mat2::Material)
    return unit_young(mat1.E, mat2.E, mat1.mu, mat2.mu)
end

"""Возвращает приведённый модуль упругости"""
function unit_young(gear::CircularGearing)
    return unit_young(gear.pin_material, gear.satellite_material)
end

"""Считает линейную контактную жёсткость"""
function linear_contact_stiffness(gear::CircularGearing)::Float64
    return unit_young(gear) * gear.b_w * π/4
end

function herz_stress(force, curvature, E_unit, b)
    @assert force >= 0
    @assert curvature >= 0
    return sqrt(force * curvature * E_unit / π / b)
end