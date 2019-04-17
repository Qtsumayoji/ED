using LinearAlgebra
using SparseArrays
using Distributions
using PyCall

@pyimport matplotlib.pyplot as plt
@pyimport seaborn as sns

include("./Model.jl")
include("./Krylov.jl")
include("./Fermion.jl")
include("./Spectrum.jl")
include("./Parameter.jl")
include("./Make_singleband_Hubbard.jl")
#BLAS.set_num_threads(2)

t = 1.0
t2 = 0.0#-0.35
U = 10*t
V = 1.5
μ = 0.0

Make_singleband_Hubbard.main(t, t2, U, V, μ)

function plot_spectra(system_para, Nk, NΩ, Ω, G)
    g = system_para.reciprocal_lattice_vec
    g1 = g[1]

    X = zeros(Nk + 1, NΩ)
    Y = zeros(Nk + 1, NΩ)
    for i in 1:NΩ
        for m in 1:Nk + 1
            X[m, i] = m - 1.5#(m-1)/Ns*g1[1] - 1/Ns*g1[1]/2
            Y[m, i] = Ω[i]
        end
    end

    plt.figure(figsize=(10,8))
    plt.pcolormesh(X, Y, -1.0/pi*imag(G))
    plt.colorbar()

    plt.figure(figsize=(10,8))
    Iω = zeros(NΩ)
    for i in 1:NΩ
        for m in 1:Nk
            Iω[i] += -1.0/pi*imag(G[m, i])
        end
    end
    plt.plot(Ω, Iω)
end

function main_spectral_func(Egs, φgs, system_para, basis)
    println("calc_spectral_func...")
    g1 = system_para.reciprocal_lattice_vec[1]
    η = 0.1*t
    n_lanczos_vec = 200
    Nk = system_para.Ns
    NΩ = 1000
    Ω = range(-12, stop=12, length=NΩ)
    G = zeros(Complex, Nk, NΩ)
    spectral_func_para = Parameter.Spectral_func_para(η, n_lanczos_vec, Nk, NΩ, Ω, G)

    for m in 1:Nk
        q = (m-1)*g1/Nk
        Spectrum.calc_spectral_func(m, q, Egs, φgs, system_para, spectral_func_para, basis)
    end

    plot_spectra(system_para, Nk, NΩ, Ω, G)
end

function main_RIXS(Egs, φgs, system_para, basis)
    g = system_para.reciprocal_lattice_vec
    g1 = g[1]
    Vd = 15.0*t
    η = 0.2*t
    Γ = t
    n_lanczos_vec = 200
    # 入射光の振動数
    ωin = -14.0*t
    
    Ns = system_para.Ns
    NΩ = 1000
    # ωin - ωout
    Ω = range(0.0, stop=15.0, length=NΩ)
    # x-ray の運動量変化
    NQ = system_para.Ns + 1
    Q = zeros(NQ)
    G = zeros(Complex, NQ, NΩ)
    RIXS_para = Parameter.RIXS_para(Vd, η, Γ, n_lanczos_vec, ωin, NΩ, Ω, NQ, Q, G)
    
    H = @time Fermion.calc_extHubbard_model(system_para, basis)
    @time for m in 1:NQ
        # 外場の波数
        q = (m-1)/Ns*g1 - g1/2
       @time Spectrum.calc_RIXS_spectrum(m, q, Egs, φgs, H, system_para, RIXS_para, basis)
    end

    #plt.subplot()
    plt.figure(figsize=(10,8))
    sns.set_palette("hsv",NQ)
    for i in 1:NQ
        b = [i*0.5 for j in 1:NΩ]
        Intensity = 1.0/π*imag(G[i, :]) + b
        plt.plot(Ω, Intensity,label=string(i))
        plt.legend(loc="best")
    end
    plot_spectra(system_para, NQ, NΩ, Ω, -G)
end

function main_dynamical_structure_factor(Egs, φgs, system_para, basis)
    println("calc_dynamical_structure_factor")
    g1 = system_para.reciprocal_lattice_vec[1]
    η = 0.2*t
    n_lanczos_vec = 200
    Nk = 80
    NΩ = 1000
    Ω = range(0.0, stop=15.0, length=NΩ)
    G = zeros(Complex, Nk, NΩ)
    DSF_para = Parameter.Dynamical_structure_factor_para(η, n_lanczos_vec, Nk, NΩ, Ω, G)

    H = Fermion.calc_extHubbard_model(system_para, basis)
    for m in 1:Nk
        q = (m-1)/(Nk-1)*[1pi;0;0]
        Spectrum.calc_dynamical_structure_factor(m, q, Egs, φgs, H, system_para, DSF_para, basis)
    end

    sns.set_palette("hsv",Nk)
    for i in 1:Nk
        b = [i*0.5 for j in 1:NΩ]
        plt.plot(Ω, -1.0/pi*imag(G[i,:]),label=string(i))
        plt.legend(loc="best")
    end

    plot_spectra(system_para, Nk, NΩ, Ω, G)
end

function Hubbard_model()
    # 格子モデルの作成
    #link_mat, link_list, pos = Model.square_lattice(Nx, Ny)
    input_file = "1d_extHubbard.txt"
    link_mat, link_list, pos, system_size = Model.read_model(input_file)
    unit_vec = Model.get_square_lattice_unit_vec()

    ns = system_size[1]
    Nx = system_size[2]
    Ny = system_size[3]
    Nz = system_size[4]
    Ns = system_size[5]
    # 電子数
    Ne = Ns

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

    println("calc_Hubbard_model...")
    #H = @time Fermion.calc_Hubbard_model(H_para, system_para, basis)
    H = Fermion.calc_extHubbard_model(system_para, basis)

    println("Diagonalization...")
    E, tri = @time Krylov.lanczos(H, minite=200, maxite=3000, ϵ=1E-10, nev=1)
    #A = zeros(dim, dim)
    #A .= H
    #H = []
    #E ,U = eigen(A)
    Egs = E[1]
    println("Egs/Ns=",Egs/Ns)

    println("inverse_iteration...")
    φgs = @time Krylov.inverse_iteration(H, Egs, maxite=2000, ϵ_inv=1E-10, ϵ_cg=1E-7)
    φgs /= norm(φgs)

    # <niσ>
    #for i in 1:Ns
    #    φ = zeros(dim)
    #    for j in 1:dim
    #        state = basis[j]
    #        φ[j] = Fermion.number_op(i, Ns, state)*φgs[j]
    #    end
    #    println("<n",i,">=",φgs'*φ)
    #end

    #println("calc_charge_correlation_func")
    #@time Fermion.calc_charge_correlation_func(φgs, system_para, basis)

    main_spectral_func(Egs, φgs, system_para, basis)
    plt.show()

    #main_RIXS(Egs, φgs, system_para, basis)
    #plt.show()

    #main_dynamical_structure_factor(Egs, φgs, system_para, basis)
    #plt.show()

    #sns.set_style("white")
    #plt.figure(figsize=(10, 8))
    #plt.subplot(211)
    #@time plt.plot(1:length(E), E)
    #plt.subplot(212)
    #@time plt.hist(E, bins = 200)
    #plt.show()
end

Hubbard_model()
