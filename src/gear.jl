include("materials.jl")

"""Основная структура для работы с библиотекой.

Параметры

'z_p'

'a_p'

'e'

'd_p'

's'

'z_c'

'λ'

'satellite_material'

'pin_material'

'c'

'b_w'
"""
mutable struct CircularGearing
    z_p::Int32
    a_p::Float64
    e::Float64
    d_p::Float64
    s::Int8
    z_c::Int32
    λ::Float64
    satellite_material::Union{Nothing, Material}
    pin_material::Union{Nothing, Material}
    c::Int8
    b_w::Float64

    """Конструктор для передачи с параметрами `z_p, a_p, λ, d_p`"""
    
    function CircularGearing(;z_p, a_p, λ, d_p, s=1, c=1, b_w=1)
        e = λ * a_p / (2 * z_p)
        z_c = z_p - s
        return new(z_p, a_p, e, d_p, s, z_c, λ, nothing, nothing, c, b_w)
    end

end

"""
`z_p` - число цевок
`a_p` - делительный диаметр
`e` - эксцентриситет
`d_p` - диаметр цевки
`s` - параметр эпи/гипо
`c` - параметр обычная/обратная
"""
function CircularGearingByExcentricity(;z_p, a_p, e, d_p, s=1, c=1, b_w=1)
    λ = 2 * z_p * e / a_p
    z_c = z_p - s
    return CircularGearing(z_p = z_p, a_p=a_p, λ=λ, d_p=d_p, b_w=b_w, s=s, c=c)
end

"""Возвращает передачу с параметрами z_p = 24, a_p = 100, λ = 0.7, d_p=6, b_w=6"""
function DefaultCircularGearing()
    return CircularGearing(;z_p = 24, a_p = 100, λ = 0.7, d_p=6, b_w=6)
end

function Base.display(x::CircularGearing)
    type = x.s > 0 ? "epicycloid" : "hypocycloid"
    text = "Gearing(z_p=$(x.z_p), a_p=$(x.a_p), e=$(x.e), d_p=$(x.d_p), λ=$(x.λ), type=$(type))"
    return text
end

"""Возвращает координаты точки циклоиды передачи `gear`, определённой параметром `t`"""
function cycloid_point(t::Number, gear::CircularGearing; t_01=0, t_02=0)
    a_p = gear.a_p
    e = gear.e
    z_p = gear.z_p
    s = gear.s
    Cx = s * a_p / 2 * sin(t + t_01) - e * sin(z_p * t + t_02)
    Cy = a_p / 2 * cos(t + t_01) - e * cos(z_p * t + t_02)
    return [Cx, Cy]
end

"""Возвращает матрицу с координатами точкек циклоиды передачи `gear`, определённых вектором `t`"""
function cycloid_points(t::AbstractVector, gear::CircularGearing; t_01=0, t_02=0)::Matrix{Float64}
    a_p = gear.a_p
    e = gear.e
    z_p = gear.z_p
    s = gear.s
    n = length(t)
    res = Matrix{Float64}(undef, n, 2)
    for i in 1:n
        Cx = s * a_p / 2 * sin(t[i] + t_01) - e * sin(z_p * t[i] + t_02)
        Cy = a_p / 2 * cos(t[i] + t_01) - e * cos(z_p * t[i] + t_02)
        res[i, 1] = Cx
        res[i, 2] = Cy
    end
    return res
end

"""Возвращает нормаль точки циклоиды передачи `gear`, определённой параметром `t`"""
function cycloid_normal(t::Number, gear::CircularGearing; t_01=0, t_02=0)
    λ = gear.λ
    s = gear.s
    z_p = gear.z_p
    sqr = sqrt(1 + λ^2 - 2 * λ * s * cos((z_p - s)*t - s * t_01 + t_02))
    Nx = (-sin(t + t_01) + λ * sin(z_p * t + t_02)) / sqr
    Ny = (-s * cos(t + t_01) + λ * cos(z_p * t + t_02)) / sqr
    return [Nx, Ny]
end

"""Возвращает нормаль точки циклоиды передачи `gear```, определённой вектором `t`"""
function cycloid_normals(t::AbstractVector, gear::CircularGearing; t_01=0, t_02=0)
    λ = gear.λ
    s = gear.s
    z_p = gear.z_p
    n = length(t)
    res = Matrix{Float64}(undef, n, 2)
    for i in 1:n
        sqr = sqrt(1 + λ^2 - 2 * λ * s * cos((z_p - s)*t[i] - s * t_01 + t_02))
        Nx = (-sin(t[i] + t_01) + λ * sin(z_p * t[i] + t_02)) / sqr
        Ny = (-s * cos(t[i] + t_01) + λ * cos(z_p * t[i] + t_02)) / sqr
        res[i, :] .= Nx, Ny
    end
    return res
end

"""Возвращает координаты точки сателлита передачи `gear`, определённой параметром `t`"""
function satellite_point(t::Number, gear::CircularGearing; t_01=0, t_02=0)
    Cx, Cy = cycloid_point(t, gear, t_01=t_01, t_02=t_02)
    Nx, Ny = cycloid_normal(t, gear, t_01=t_01, t_02=t_02)
    d_p = gear.d_p
    Px = Cx + d_p / 2 * Nx * gear.c
    Py = Cy + d_p / 2 * Ny * gear.c
    return [Px, Py]
end

"""Возвращает координаты точки сателлита передачи `gear`, определённой параметром `t`"""
function satellite_points(t::AbstractVector, gear::CircularGearing; t_01=0, t_02=0)
    cp = cycloid_points(t, gear, t_01=t_01, t_02=t_02)
    cn = cycloid_normals(t, gear, t_01=t_01, t_02=t_02)
    d_p = gear.d_p
    n = length(t)
    res = Matrix{Float64}(undef, n, 2)
    for i in 1:n
        Px = cp[i, 1] + d_p / 2 * cn[i, 1] * gear.c
        Py = cp[i, 2] + d_p / 2 * cn[i, 2] * gear.c
        res[i, :] .= Px, Py
    end
    return res
end


"""Возвращает кривизну циклоиды передачи `gear` в точке, заданной параметром `t`"""
function cycloid_curvature(t::Number, gear::CircularGearing; t_01=0, t_02=0)
    s = gear.s
    λ = gear.λ
    z_p = gear.z_p
    z_c = gear.z_c
    a_p = gear.a_p
    k_c = - 2 / a_p * 
    (
        s + λ^2 * z_p - λ * cos(z_c*t - s * t_01 + t_02)*(1 + s*z_p)
    )/
    (
        sqrt(1 + λ^2 - 2 * s * λ * cos(z_c * t - s * t_01 + t_02))^3
    )
    return k_c
end

"""Возвращает кривизну сателлита передачи `gear` в точке, заданной параметром `t`"""
function satellite_curvature(t::Number, gear::CircularGearing; t_01=0, t_02=0)
    d_p = gear.d_p
    s = gear.s
    c = gear.c
    k_c = cycloid_curvature(t, gear, t_01=t_01, t_02=t_02)
    k_p_inv = 1 / k_c + d_p / 2
    k_p = 1 / k_p_inv
    return k_p
end

"""Возвращает приведённую кривизну сателлита передачи `gear`` в точке, заданной параметром `t`"""
function unit_curvature(t, gear::CircularGearing; t_01=0, t_02=0)
    s = gear.s
    d_p = gear.d_p
    c = gear.c
    k_p = satellite_curvature(t, gear, t_01=t_01, t_02=t_02)
    k_u = 2 / d_p - k_p
    return k_u
end

