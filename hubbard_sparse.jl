using LinearAlgebra
using SparseArrays
using Distributions
using PyCall

@pyimport pylab as plt
@pyimport seaborn as sns

include("./model.jl")
include("./Krylov.jl")
include("./Fermion.jl")
include("./Spectrum.jl")
include("./Parameter.jl")

#unit cell内のサイト数
const ns = 1
const Nx = 8
const Ny = 1
const Ns = Nx*Ny*ns
# 電子数
const Ne = Ns 

function plot_spectral_func(system_para, spectral_func_para)
    g = system_para.reciprocal_lattice_vec
    g1 = g[1]
    Nk = spectral_func_para.Nk
    NΩ = spectral_func_para.NΩ
    Ω = spectral_func_para.Ω
    G = spectral_func_para.G

    X = zeros(Nk + 1, NΩ)
    Y = zeros(Nk + 1, NΩ)
    for i in 1:NΩ
        for m in 1:Nk + 1
            X[m, i] = m - 1.5#(m-1)/Ns*g1[1] - 1/Ns*g1[1]/2
            Y[m, i] = Ω[i]
        end
    end

    plt.figure(figsize=(10,8))
    plt.pcolormesh(X, Y, -1.0/pi*imag(G), cmap="magma")
    plt.colorbar()
    plt.show()

    Iω = zeros(NΩ)
    for i in 1:NΩ
        for m in 1:Nk
            Iω[i] += -1.0/pi*imag(G[m, i])
        end
    end
    plt.plot(Ω, Iω)
    plt.show()
end

function plot_dynamical_structure_factor(system_para, dynamical_structure_factor_para)
    g = system_para.reciprocal_lattice_vec
    g1 = g[1]
    Nk = dynamical_structure_factor_para.Nk
    NΩ = dynamical_structure_factor_para.NΩ
    Ω = dynamical_structure_factor_para.Ω
    G = dynamical_structure_factor_para.G

    X = zeros(Nk + 1, NΩ)
    Y = zeros(Nk + 1, NΩ)
    for i in 1:NΩ
        for m in 1:Nk + 1
            X[m, i] = m - 1.5#(m-1)/Ns*g1[1] - 1/Ns*g1[1]/2
            Y[m, i] = Ω[i]
        end
    end

    plt.figure(figsize=(10,8))
    plt.pcolormesh(X, Y, -1.0/pi*imag(G), cmap="magma")
    plt.colorbar()
    plt.show()

    Iω = zeros(NΩ)
    for i in 1:NΩ
        for m in 1:Ns
            Iω[i] += -1.0/pi*imag(G[m, i])
        end
    end
    plt.plot(Ω, Iω)
    plt.show()
end

function test()
    # 格子モデルの作成
    link_mat, link_list, pos = Model.square_lattice(Nx, Ny)
    unit_vec = Model.get_square_lattice_unit_vec()

    # 逆格子ベクトル
    V = dot(unit_vec[1], cross(unit_vec[2], unit_vec[3]))
    g1 = 2.0*pi*cross(unit_vec[2], unit_vec[3])/V
    g2 = 2.0*pi*cross(unit_vec[3], unit_vec[1])/V
    g3 = 2.0*pi*cross(unit_vec[1], unit_vec[2])/V
    g = [zeros(3) for i in 1:3]
    g[1] = g1; g[2] = g2; g[3] = g3

    system_para = Parameter.System_para(ns, Nx, Ny, Ns, Ne, unit_vec, g, link_mat, link_list, pos)
    Model.show_links(link_mat)

    # 電子数とtotSzの保存した基底を作成
    basis = Fermion.make_n_basis(Ns, Ne)
    nupspin = div(Ne, 2)
    ndownspin = Ne - nupspin
    basis = Fermion.make_s_basis(Ns, nupspin, ndownspin, basis)

    # Hubbard modelのパラメータ
    t = 0.35; U = 10.0*t; μ = U/2
    H_para = Parameter.Hubbard_para(t, U, μ, 0.0)

    println("calc_Hubbard_model...")
    H = @time Fermion.calc_Hubbard_model(H_para, system_para, basis)

    println("Diagonalization...")
    E, tri = @time Krylov.lanczos(H, minite=200, maxite=3000, ϵ=1E-10, nev=1)

    Egs = E[1]
    println("Egs/Ns=",Egs/Ns)

    println("inverse_iteration...")
    φgs = @time Krylov.inverse_iteration(H, Egs, maxite=2000, ϵ_inv=1E-8, ϵ_cg=1E-5)
    φgs /= norm(φgs)

    #println("calc_charge_correlation_func")
    #@time Fermion.calc_charge_correlation_func(φgs, system_para, basis)

    println("calc_spectral_func...")
    η = 1e-2
    n_lanczos_vec = 200
    Nk = 50
    NΩ = 500
    Ω = range(-U, stop=U, length=NΩ)
    G = zeros(Complex, Nk, NΩ)
    spectral_func_para = Parameter.Spectral_func_para(η, n_lanczos_vec, Nk, NΩ, Ω, G)
    
    for m in 1:Nk
        # Γ-M
        q = (m-1)*g1/Nk
        # Γ-K
        #q = (m-1)*g1/Nx + (m-1)*g2/Ny
        Spectrum.calc_spectral_func(m, q, Egs, φgs, H_para, system_para, spectral_func_para, basis)
    end
    plot_spectral_func(system_para, spectral_func_para)

    println("calc_dynamical_structure_factor")
    η = 1e-2
    n_lanczos_vec = 200
    Nk = 100
    NΩ = 1000
    Ω = range(-U/100, stop=U/100, length=NΩ)
    G = zeros(Complex, Nk, NΩ)
    DSF_para = Parameter.Dynamical_structure_factor_para(η, n_lanczos_vec, Nk, NΩ, Ω, G)

    H = Fermion.calc_Hubbard_model(H_para, system_para, basis)
    for m in 1:Nk
        q = (m - 1)/Nk*g1
        Spectrum.calc_dynamical_structure_factor(m, q, Egs, φgs, H, system_para, DSF_para, basis)
    end

    plot_dynamical_structure_factor(system_para, DSF_para)

    #sns.set_style("white")
    #plt.figure(figsize=(10, 8))
    #plt.subplot(211)
    #@time plt.plot(1:length(E), E)
    #plt.subplot(212)
    #@time plt.hist(E, bins = 200)
    #plt.show()
end

test()