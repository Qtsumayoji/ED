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
include("./bit_operations.jl")

PyPlot.rc("font",family ="Times New Roman")

BLAS.set_num_threads(4)

Udd = 8.0
Upp = Udd/2.0
#Upc = Upp/0.8
Upc = 5.0

function main_RIXS(Egs, φgs, system_para, basis)
    g = system_para.reciprocal_lattice_vec
    g1 = g[1]
    ηs = [0.25]
    for η in ηs
    #η = 0.5
    Γ = 0.5
    n_lanczos_vec = 200
    # 入射光の振動数
    ωin = 4.2
    
    Ns = system_para.Ns
    NΩ = 2000
    # ωin - ωout
    Ω = range(0.0, stop=-10.0, length=NΩ)
    # x-ray の運動量変化
    NQ = 1
    Q = zeros(NQ)
    G = zeros(Complex, NQ, NΩ)

    #RIXS_para = Parameter.RIXS_para(-Upc, η, Γ, n_lanczos_vec, ωin, NΩ, -Ω, NQ, Q, G)
    #H = Fermion.calc_3band_dp_model(system_para, basis)
    #@time Spectrum.calc_1s4p_RIXS_spectrum(1, [0.0;0.0;0.0], Egs, φgs, H, system_para, RIXS_para, basis)

    RIXS_para = Parameter.RIXS_para(Upc, η, Γ, n_lanczos_vec, ωin, NΩ, Ω, NQ, Q, G)
    @time Spectrum.calc_O2p_RIXS_spectrum(1, Egs, φgs, system_para, RIXS_para, basis)

    #PyPlot.subplot()
    Intensity = 1.0/π*imag(G[1, :])
    PyPlot.plot(Ω, Intensity)
    PyPlot.tight_layout()
    end
end

function main_XAS(Egs, φgs, system_para, basis)
    g1 = system_para.reciprocal_lattice_vec[1]
    Γ = 0.25
    n_lanczos_vec = 200
    Ns = system_para.Ns
    NQ = 1
    Q = zeros(NQ)
    NΩ = 10000
    Ω = range(0.0, stop=10.0, length=NΩ)
    G = zeros(Complex, NQ, NΩ)
    XAS_para = Parameter.XAS_para(Upc, Γ, n_lanczos_vec, NΩ, Ω, NQ, Q, G)

    Spectrum.calc_O2p_XAS_spectrum(1, Egs, φgs, system_para, XAS_para, basis)

    PyPlot.figure(figsize=(10, 6))
    PyPlot.xlabel("Energy", size=20)
    PyPlot.ylabel("Intensity", size=20)
    PyPlot.tight_layout()
    PyPlot.plot(Ω, 1.0/pi*imag(G[1, :]), c="black") 
end

function main_dynamical_structure_factor(Egs, φgs, system_para, basis)
    g1 = system_para.reciprocal_lattice_vec[1]
    η = 0.2
    n_lanczos_vec = 200
    Ns = system_para.Ns
    Nk = 1
    K = zeros(Nk)
    NΩ = 1000
    Ω = range(0.0, stop=-15.0, length=NΩ)
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
    input_file = "tanioka.txt"
    #input_file = "corner_sharing_dp_deg.txt"
    link_mat, link_list, pos, system_size = Model.read_model(input_file)
    unit_vec = Model.get_square_lattice_unit_vec()

    ns = system_size[1]
    Nx = system_size[2]
    Ny = system_size[3]
    Nz = system_size[4]
    Ns = system_size[5]

    # 逆格子ベクトル
    V = dot(unit_vec[1], cross(unit_vec[2], unit_vec[3]))
    g1 = 2.0*pi*cross(unit_vec[2], unit_vec[3])/V
    g2 = 2.0*pi*cross(unit_vec[3], unit_vec[1])/V
    g3 = 2.0*pi*cross(unit_vec[1], unit_vec[2])/V
    g = [g1, g2, g3]

    # ホール数
    Ne = 2
    # 電子数=Ne=2, totSz=nup - ndown=0 の基底を作成
    nup   = 1
    ndown = 1
    #basis = bit_operations.tow_particle_basis(Ns, 0)
    basis = bit_operations.two_particle_Sz_0_basis(Ns)
    #basis = bit_operations.one_particle_basis(Ns,0)
    #for i in basis println(string(i, base=2)) end

    system_para = Parameter.System_para(ns, Nx, Ny, Nz, Ns, Ne, nup, ndown, unit_vec, g, link_mat, link_list, pos)
    Model.show_links(link_mat)

    println("calc_dp_model...")
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
end

function dp_model2()
    # 格子モデルの作成
    #link_mat, link_list, pos = Model.square_lattice(Nx, Ny)
    input_file = "tanioka.txt"
    #input_file = "corner_sharing_dp_deg.txt"
    link_mat, link_list, pos, system_size = Model.read_model(input_file)
    unit_vec = Model.get_square_lattice_unit_vec()

    ns = system_size[1]
    Nx = system_size[2]
    Ny = system_size[3]
    Nz = system_size[4]
    Ns = system_size[5]

    # 逆格子ベクトル
    V = dot(unit_vec[1], cross(unit_vec[2], unit_vec[3]))
    g1 = 2.0*pi*cross(unit_vec[2], unit_vec[3])/V
    g2 = 2.0*pi*cross(unit_vec[3], unit_vec[1])/V
    g3 = 2.0*pi*cross(unit_vec[1], unit_vec[2])/V
    g = [g1, g2, g3]

    # ホール数
    Ne = 2
    # 電子数=Ne=2, totSz=nup - ndown=0 の基底を作成
    nup   = 1
    ndown = 1
    basis = bit_operations.tow_particle_basis(Ns, 0)

    system_para = Parameter.System_para(ns, Nx, Ny, Nz, Ns, Ne, nup, ndown, unit_vec, g, link_mat, link_list, pos)
    Model.show_links(link_mat)

    println("calc_dp_model...")
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
end


function main()
    PyPlot.figure(figsize=(12, 6))
    dp_model()
    dp_model2()
    PyPlot.show()
end
main()