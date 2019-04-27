module Parameter

    struct System_para
        ns::Int64
        Nx::Int64
        Ny::Int64
        Nz::Int64
        Ns::Int64
        Ne::Int64
        nup::Int64
        ndown::Int64
        unit_vec::Array{Array{Float64, 1}}
        reciprocal_lattice_vec::Array{Array{Float64, 1}}
        link_mat::Array{Array{Any,1},1}
        link_list::Array{Array{Any}, 1}
        pos::Array{Array{Float64}}
    end

    mutable struct Hubbard_para
        t1::Float64
        t2::Float64
        t3::Float64
        U::Float64
        μ::Float64
        V::Float64
    end

    mutable struct Spectral_func_para
        η::Float64
        n_lanczos_vec::Int64
        Nk::Int64
        NΩ::Int64
        Ω::Array{Float64, 1}
        G::Array{Complex, 2}
    end

    mutable struct Dynamical_structure_factor_para
        η::Float64
        n_lanczos_vec::Int64
        Nk::Int64
        NΩ::Int64
        Ω::Array{Float64, 1}
        G::Array{Complex, 2}
    end

    mutable struct RIXS_para
        Vd::Float64
        η::Float64
        Γ::Float64
        n_lanczos_vec::Int64
        # 入射光の振動数
        ωin::Float64
        NΩ::Int64
        # ωin - ωout
        Ω::Array{Float64, 1}
        # x-ray の運動量変化
        NQ::Int64
        Q::Array{Float64, 1}
        G::Array{Complex, 2}
    end
end