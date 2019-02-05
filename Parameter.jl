module Parameter

    struct System_para
        ns::Int64
        Nx::Int64
        Ny::Int64
        Ns::Int64
        Ne::Int64
        unit_vec::Array{Array{Float64, 1}}
        reciprocal_lattice_vec::Array{Array{Float64, 1}}
        link_mat::Array{Int64, 2}
        link_list::Array{Array{Any}, 1}
        pos::Array{Array{Float64}}
    end

    mutable struct Hubbard_para
        t::Float64
        U::Float64
        μ::Float64
        Vd::Float64
    end

    mutable struct Spectral_func_para
        η::Float64
        n_lanczos_vec::Int64
        Nb::Int64
        NΩ::Int64
        Ω::Array{Float64, 1}
        G::Array{Complex, 2}
    end

    mutable struct RIXS_para
        Vd::Float64
        # 入射光の振動数
        NΩ::Int64
        Ω::Array{Float64, 1}
        # x-ray の運動量変化
        NQ::Int64
        Q::Array{Float64, 1}
        G::Array{Complex, 2}
    end
end