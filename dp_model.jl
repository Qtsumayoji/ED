using LinearAlgebra
using PyPlot
using PyCall
#import Arpack

sns = pyimport("seaborn")

include("./Model.jl")
include("./Krylov.jl")
include("./Fermion.jl")
include("./Spectrum.jl")
include("./Parameter.jl")
include("./Make_singleband_Hubbard.jl")

PyPlot.rc("font",family ="Times New Roman")

BLAS.set_num_threads(4)

Udd = 8.0
Upp = Udd/2.0
Upc = Upp/0.8

function plot_spectra(Nk, K, NΩ, Ω, G)
    X = zeros(Nk + 1, NΩ + 1)
    Y = zeros(Nk + 1, NΩ + 1)
    ΔX = K[2] - K[1]
    ΔY = Ω[2] - Ω[1]
    for i in 1:NΩ + 1
        for m in 1:Nk + 1
            X[m, i] = K[1] + (m - 1.5)*ΔX
            Y[m, i] = Ω[1] + (i - 1.5)*ΔY
        end
    end

    PyPlot.figure(figsize=(8,6))
    PyPlot.pcolormesh(X, Y, -1.0/pi*imag(G))
    PyPlot.colorbar()
    PyPlot.xlabel("k", size=20)
    PyPlot.ylabel("ω", size=20)
    PyPlot.tight_layout()


    Iω = zeros(NΩ)
    for i in 1:NΩ
        # B.Zの両端の波数をカラープロットで見ているので片側を省くための-1
        for m in 1:Nk - 1
            Iω[i] += -1.0/pi*imag(G[m, i])
        end
    end
    PyPlot.figure(figsize=(8,6))
    PyPlot.xlabel("Energy", size=20)
    PyPlot.ylabel("Intensity", size=20)
    PyPlot.tight_layout()
    PyPlot.plot(Ω, Iω, c="black")
end

function main_RIXS(Egs, φgs, system_para, basis)
    g = system_para.reciprocal_lattice_vec
    g1 = g[1]
    η = 0.5
    Γ = 0.5
    n_lanczos_vec = 200
    # 入射光の振動数
    ωin = 4.85
    
    Ns = system_para.Ns
    NΩ = 1000
    # ωin - ωout
    Ω = range(0.0, stop=-15.0, length=NΩ)
    # x-ray の運動量変化
    NQ = 1
    Q = zeros(NQ)
    G = zeros(Complex, NQ, NΩ)
    RIXS_para = Parameter.RIXS_para(Upc, η, Γ, n_lanczos_vec, ωin, NΩ, -Ω, NQ, Q, G)
    
    @time Spectrum.calc_O2p_RIXS_spectrum(1, Egs, φgs, system_para, RIXS_para, basis)

    #PyPlot.subplot()
    PyPlot.figure(figsize=(10, 6))
    Intensity = 1.0/π*imag(G[1, :])
    PyPlot.plot(Ω, Intensity)
    PyPlot.tight_layout()
end

function main_XAS(Egs, φgs, system_para, basis)
    g1 = system_para.reciprocal_lattice_vec[1]
    η = 0.5
    Γ = 0.1
    ωin = 0.0
    n_lanczos_vec = 200
    Ns = system_para.Ns
    Nk = 1
    K = zeros(Nk)
    NΩin = 10000
    Ωin = range(0.0, stop=20.0, length=NΩin)
    G = zeros(Complex, Nk, NΩin)
    RIXS_para = Parameter.RIXS_para(Upc, η, Γ, n_lanczos_vec, ωin, NΩin, Ωin, Nk, K, G)

    Spectrum.calc_O2p_XAS_spectrum(Egs, φgs, system_para, RIXS_para, basis)

    PyPlot.figure(figsize=(10, 6))
    PyPlot.xlabel("Energy", size=20)
    PyPlot.ylabel("Intensity", size=20)
    PyPlot.tight_layout()
    PyPlot.plot(Ωin, 1.0/pi*imag(G[1, :]), c="black") 
end

function main_dynamical_structure_factor(Egs, φgs, system_para, basis)
    g1 = system_para.reciprocal_lattice_vec[1]
    η = 0.2
    n_lanczos_vec = 200
    Ns = system_para.Ns
    Nk = 1
    K = zeros(Nk)
    NΩ = 1000
    Ω = range(0.0, stop=15.0, length=NΩ)
    G = zeros(Complex, Nk, NΩ)
    DSF_para = Parameter.Dynamical_structure_factor_para(η, n_lanczos_vec, Nk, NΩ, Ω, G)

    H = Fermion.calc_3band_dp_model(system_para, basis)
    q = [0.0;0.0;0.0]
    Spectrum.calc_dynamical_structure_factor(1, q, Egs, φgs, H, system_para, DSF_para, basis)

    PyPlot.figure(figsize=(8,6))
    PyPlot.plot(Ω, -1.0/pi*imag(G[1,:]))
end

function dp_model()
    # 格子モデルの作成
    #link_mat, link_list, pos = Model.square_lattice(Nx, Ny)
    input_file = "corner_sharing_dp.txt"
    link_mat, link_list, pos, system_size = Model.read_model(input_file)
    unit_vec = Model.get_square_lattice_unit_vec()

    ns = system_size[1]
    Nx = system_size[2]
    Ny = system_size[3]
    Nz = system_size[4]
    Ns = system_size[5]
    # 電子数
    Ne = 16

    # 逆格子ベクトル
    V = dot(unit_vec[1], cross(unit_vec[2], unit_vec[3]))
    g1 = 2.0*pi*cross(unit_vec[2], unit_vec[3])/V
    g2 = 2.0*pi*cross(unit_vec[3], unit_vec[1])/V
    g3 = 2.0*pi*cross(unit_vec[1], unit_vec[2])/V
    g = [g1, g2, g3]

    system_para = Parameter.System_para(ns, Nx, Ny, Nz, Ns, Ne, unit_vec, g, link_mat, link_list, pos)
    Model.show_links(link_mat)

    # 電子数=Ne,totSz=0の基底を作成
    basis = Fermion.make_n_basis(Ns, Ne)
    nupspin = div(Ne, 2)
    ndownspin = Ne - nupspin
    basis = Fermion.make_s_basis(Ns, nupspin, ndownspin, basis)
    dim = length(basis)

    println("calc_dp_model...")
    #H = @time Fermion.calc_dp_model(H_para, system_para, basis)
    H = Fermion.calc_3band_dp_model(system_para, basis)

    println("Diagonalization...")
    E, tri = @time Krylov.lanczos(H, minite=200, maxite=3000, ϵ=1E-10, nev=1)
    Egs = E[1]
    println("Egs=",Egs)
    println("Egs/Ns=",Egs/Ns)

    println("inverse_iteration...")
    φgs = @time Krylov.inverse_iteration(H, Egs, maxite=2000, ϵ_inv=1E-10, ϵ_cg=1E-7)
    φgs /= norm(φgs)

    #main_XAS(Egs, φgs, system_para, basis)
    main_RIXS(Egs, φgs, system_para, basis)
    #main_dynamical_structure_factor(Egs, φgs, system_para, basis)
    
    PyPlot.show()
end

dp_model()
