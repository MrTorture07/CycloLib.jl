include("gear.jl")
include("utils.jl")

"""
Сателлиты идут от верхнего по часовой стрелке
   
Параметры

    `gear`

    `z_s`

    `sector_angle`
"""
mutable struct SectorGearing
    gear::CircularGearing
    satellite_quantity::Int32
    sector_angle::Float64
    # excenter_center::Vector{Float64}
    
    """
    Параметры
    `gear`
    `z_s`
    `sector_angle`
    """
    function SectorGearing(;gear, z_s, sector_angle)
        # gear.z_p % z_s ==0 ? 1 : throw(ArgumentError("Число цевок должно быть кратно числу сателлитов"))
        return new(gear, z_s, sector_angle)
    end
end

function excenter_coords(sgear::SectorGearing, current_angle::Number, satellite_number::Int)
    angle = excenter_angle(sgear, current_angle, satellite_number)
    e = sgear.gear.e
    x_e = e * sin(angle)
    y_e = e * cos(angle)
    return [x_e, y_e]
end

function excenter_angle(sgear::SectorGearing, current_angle, satellite_number::Int)
    return current_angle + (satellite_number - 1) * 2π / sgear.satellite_quantity
end

function sector_satellite_points(sgear::SectorGearing, satellite_number::Int; N=1000)
    θ = sgear.sector_angle
    θ_1, θ_2 = -θ/2, θ/2
    
    t = collect(range(θ_1, θ_2, N))

    exc_phase = excenter_angle(sgear, 0, satellite_number)

    # Создаём вектор с координатами сателлита
    s_points = satellite_points(t, sgear.gear, t_02=exc_phase)

    return s_points
end

"""Возвращает массив с точками контакта сателлита  №`satellite_number`"""
function contact_points(sgear::SectorGearing, satellite_number::Int, current_angle::Number)

    gear = sgear.gear
    s = gear.s
    z_p = gear.z_p
    c_p = zeros(z_p, 2)

    exc_phase = excenter_angle(sgear, 0, satellite_number)

    for i in 1:z_p
        t_contact = 2π * (i - 1) / z_p + current_angle / z_p
        c_p[i, :] .= satellite_point(t_contact, gear, t_02=exc_phase)
    end
    return c_p
end

"""Возвращает массив с нормалями точек контакта сателлита  №`satellite_number`"""
function normal_contact_points(sgear::SectorGearing, satellite_number::Int, current_angle::Number)

    gear = sgear.gear
    s = gear.s
    z_p = gear.z_p
    n_p = zeros(z_p, 2)

    exc_phase = excenter_angle(sgear, 0, satellite_number)

    for i in 1:z_p
        t_contact = 2π * (i - 1) / z_p + current_angle / z_p
        n_p[i, :] .= cycloid_normal(t_contact, gear, t_02=exc_phase)
    end
    return n_p
end

function is_pin_in_contact(sgear::SectorGearing, pin_number::Int, cur_angle::Number)::Bool
    z_p = sgear.gear.z_p
    s = sgear.gear.s
    current_t = (pin_number-1) * 2π / z_p + cur_angle / z_p
    clamp_angle!(current_t)
    θ = sgear.sector_angle
    θ_1, θ_2 = -θ/2, θ/2
    return angle_is_between(current_t, θ_1, θ_2)
end

"""Возвращает массив с граничными точками сателлита. Отсчёт идёт от левой точки профиля сателлита и кончается правой точкой профиля, чтобы
замыкаться вместе с массивом sector_satellite_points"""
function satellite_border(sgear::SectorGearing, satellite_number::Int, height::Number, N::Int=50)
    θ = sgear.sector_angle
    s = sgear.gear.s
    a_p = sgear.gear.a_p
    d_p = sgear.gear.d_p
    e = sgear.gear.e
    θ_1 = -θ/2
    θ_2 = θ/2
    t_02 = 2π *(satellite_number - 1) / sgear.satellite_quantity

    p4 = satellite_point(θ_1, sgear.gear, t_02 = t_02)
    p1 = satellite_point(θ_2, sgear.gear, t_02 = t_02)

    dir1 = p1 .- [0; 0]
    dir2 = p4 .- [0; 0]

    dir1 = dir1 / hypot(dir1[1], dir1[2])
    dir2 = dir2 / hypot(dir2[1], dir2[2])

    r = a_p / 2 - s * (d_p/2 + height + e)
    r_p1 = hypot(p1...)
    r_p4 = hypot(p4...)
    side_height_1 = s * (r_p1 - r)
    side_height_2 = s * (r_p4 - r)

    # side_height_1 = abs(a_p / 2 - s * (d_p / 2 + e) - hypot(p1[1], p1[2]) + s * height)
    # true_height_2 = a_p / 2 - s * (d_p / 2 + e) - hypot(p1[1], p1[2])

    p2 = p1 .- s * dir1 * side_height_1
    p3 = p4 .- s * dir2 * side_height_2

    β1 = atan(p1[1], p1[2])
    β2 = atan(p4[1], p4[2])
    t = range(β1,  β2, N)
    
    # r_0 = hypot(p1[1], p1[2])
    # r = r_0 - s * height
    x = r .* sin.(t)
    y = r .* cos.(t)
    
    return([p1 p2 [x y]' p3 p4]')
    # return([x y])
    # return(β1, β2)

end

function full_satellite_border(sgear::SectorGearing, satellite_number::Int, height::Number, N_bor::Int=50, N_sat::Int=1000)
    sb = satellite_border(sgear, satellite_number, height, N_bor)
    sp = sector_satellite_points(sgear, satellite_number, N=N_sat)
    return [sp;sb]
end

    

"""Переводит массив с координатами точек, относящимися к разным сателлитам в систему координат редуктора"""
function translate_to_reducer_coords!(A::AbstractArray, sgear::SectorGearing; current_angle::Number=0)
    arr_size = size(A)
    for i in 1:arr_size[1]
        for j in 1:arr_size[2]
            ec = excenter_coords(sgear, current_angle, i)
            A[i, j, :] .+= ec
        end
    end
end