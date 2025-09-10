struct Material
    E::Float64
    mu::Float64
end

const mSteel40 = Material(2.1e5, 0.3);
const mKaprolonPA16 = Material(2e3, 0.44);
const mPolyethylene = Material(1.1e3, 0.42)