using LinearAlgebra
using SparseArrays
using Distributions
using PyCall

@pyimport pylab as plt
@pyimport seaborn as sns

include("./Model.jl")
include("./Krylov.jl")
include("./Fermion.jl")
include("./Spectrum.jl")
include("./Parameter.jl")

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

function plot_RIXS(system_para, RIXS_para)
    g = system_para.reciprocal_lattice_vec
    g1 = g[1]
    NQ = RIXS_para.NQ
    NΩ = RIXS_para.NΩ
    Ω = RIXS_para.Ω
    G = RIXS_para.G

    X = zeros(NQ + 1, NΩ)
    Y = zeros(NQ + 1, NΩ)
    for i in 1:NΩ
        for m in 1:NQ + 1
            X[m, i] = m - 1.5#(m-1)/Ns*g1[1] - 1/Ns*g1[1]/2
            Y[m, i] = Ω[i]
        end
    end

    plt.figure(figsize=(10,8))
    plt.pcolormesh(X, Y, 1.0/pi*imag(G), cmap="magma")
    plt.colorbar()
    plt.show()

    Iω = zeros(NΩ)
    for i in 1:NΩ
        for m in 1:NQ
            Iω[i] += 1.0/pi*imag(G[m, i])
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
        for m in 1:Nk
            Iω[i] += -1.0/pi*imag(G[m, i])
        end
    end
    plt.plot(Ω, Iω)
    plt.show()
end

function test()
    # 格子モデルの作成
    #link_mat, link_list, pos = Model.square_lattice(Nx, Ny)
    link_mat, link_list, pos, Nxyz = Model.read_model("1d_extHubbard.txt")
    unit_vec = Model.get_square_lattice_unit_vec()

    ns = Nxyz[1]
    Nx = Nxyz[2]
    Ny = Nxyz[3]
    Nz = Nxyz[4]
    Ns = ns*Nx*Ny*Nz
    # 電子数
    Ne = Ns

    # 逆格子ベクトル
    V = dot(unit_vec[1], cross(unit_vec[2], unit_vec[3]))
    g1 = 2.0*pi*cross(unit_vec[2], unit_vec[3])/V
    g2 = 2.0*pi*cross(unit_vec[3], unit_vec[1])/V
    g3 = 2.0*pi*cross(unit_vec[1], unit_vec[2])/V
    g = [zeros(3) for i in 1:3]
    g[1] = g1; g[2] = g2; g[3] = g3

    system_para = Parameter.Model_para(ns, Nx, Ny, Ns, Ne, unit_vec, g, link_mat, link_list, pos)
    Model.show_links(link_mat)

    # 電子数とtotSzの保存した基底を作成
    basis = Fermion.make_n_basis(Ns, Ne)
    nupspin = div(Ne, 2)
    ndownspin = Ne - nupspin
    basis = Fermion.make_s_basis(Ns, nupspin, ndownspin, basis)

    println("calc_Hubbard_model...")
    #H = @time Fermion.calc_Hubbard_model(H_para, system_para, basis)
    H = Fermion.calc_extHubbard_model(system_para, basis)

    println("Diagonalization...")
    E, tri = @time Krylov.lanczos(H, minite=200, maxite=3000, ϵ=1E-10, nev=1)

    Egs = E[1]
    println("Egs/Ns=",Egs/Ns)

    println("inverse_iteration...")
    φgs = @time Krylov.inverse_iteration(H, Egs, maxite=2000, ϵ_inv=1E-10, ϵ_cg=1E-7)
    φgs /= norm(φgs)

    #println("calc_charge_correlation_func")
    #@time Fermion.calc_charge_correlation_func(φgs, system_para, basis)

    println("calc_spectral_func...")
    t = 1.0; U = 10.0*t
    η = 0.1*t
    n_lanczos_vec = 200
    Nk = Nx
    NΩ = 1000
    Ω = range(-2U, stop=2U, length=NΩ)
    G = zeros(Complex, Nk, NΩ)
    spectral_func_para = Parameter.Spectral_func_para(η, n_lanczos_vec, Nk, NΩ, Ω, G)
    
    for m in 1:Nk
        # Γ-M
        q = (m-1)*g1/Nk
        # Γ-K
        #q = (m-1)*g1/Nx + (m-1)*g2/Ny
        #Spectrum.calc_spectral_func(m, q, Egs, φgs, system_para, spectral_func_para, basis)
    end
    #plot_spectral_func(system_para, spectral_func_para)

    Vd = 15.0*t
    η = 0.1*t
    Γ = t
    n_lanczos_vec = 200
    # 入射光の振動数
    ωin = -21.2*t
    NΩ = 1000
    # ωin - ωout
    Ω = range(-24, stop=-10, length=NΩ)
    # x-ray の運動量変化
    NQ = 2
    Q = range(0.0, stop=pi, length=NQ)
    G_RIXS = zeros(Complex, NQ, NΩ)
    RIXS_para = Parameter.RIXS_para(Vd, η, Γ, n_lanczos_vec, ωin, NΩ, Ω, NQ, Q, G_RIXS)
    #G_REXS = zeros(Complex, NQ, NΩ)
    #REXS_para = Parameter.RIXS_para(Vd, η, 0.0, n_lanczos_vec, ωin, NΩ, Ω, NQ, Q, G_REXS)
    H = @time Fermion.calc_extHubbard_model(system_para, basis)
    
    @time for m in 1:NQ
        # 外場の波数
        q = (m-1)/(NQ-1)*[pi;pi;0]
        Spectrum.calc_XAS_spectrum(m, q, Egs, φgs, H, system_para, RIXS_para, basis)
    end

    sns.set_palette("hsv",NQ)
    for i in 1:NQ
        plt.plot(Ω, 1.0/pi*imag(G_RIXS[i,:]),label=string(i))
        plt.legend(loc="best")
    end

    plot_RIXS(system_para, RIXS_para)
   
#    println("calc_dynamical_structure_factor")
#    η = 0.1*t
#    n_lanczos_vec = 200
#    Nk = 8
#    NΩ = 1000
#    Ω = range(2.0, stop=15.0, length=NΩ)
#    G = zeros(Complex, Nk, NΩ)
#    DSF_para = Parameter.Dynamical_structure_factor_para(η, n_lanczos_vec, Nk, NΩ, Ω, G)
#
#    H = Fermion.calc_extHubbard_model(system_para, basis)
#    for m in 1:Nk
#        q = m/Nk*[pi;pi;0]#(m-1)/Nk*g1
#        Spectrum.calc_dynamical_structure_factor(m, q, Egs, φgs, H, system_para, DSF_para, basis)
#    end
#    sns.set_palette("hsv",Nk)
#    for i in 1:Nk
#        plt.plot(Ω, -1.0/pi*imag(G[i,:]),label=string(i))
#        plt.legend(loc="best")
#    end
#
#    plot_dynamical_structure_factor(system_para, DSF_para)

    #sns.set_style("white")
    #plt.figure(figsize=(10, 8))
    #plt.subplot(211)
    #@time plt.plot(1:length(E), E)
    #plt.subplot(212)
    #@time plt.hist(E, bins = 200)
    #plt.show()
end

test()