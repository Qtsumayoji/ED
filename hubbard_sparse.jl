using LinearAlgebra
using SparseArrays
using Distributions
using PyCall

@pyimport pylab as plt
@pyimport seaborn as sns

include("./model.jl")
include("./Krylov.jl")
include("./Fermion.jl")

#unit cell内のサイト数
const ns = 1
const Nx = 8
const Ny = 1
const Ns = Nx*Ny*ns
# 電子数
const Ne = Ns

function test()
    link_mat, link_list, pos = Model.square_lattice(Nx, Ny)
    unit_vec = Model.get_square_lattice_unit_vec()
    Model.show_links(link_mat)

    basis = Fermion.make_n_basis(Ns, Ne)
    nupspin = div(Ne, 2)
    ndownspin = Ne - nupspin
    basis = Fermion.make_s_basis(Ns, nupspin, ndownspin, basis)

    t = 1.0; U = 5.0; μ = U/2
    H_parameter = [t; U; μ]
    println("calc_Hubbard_model...")
    H = Fermion.calc_Hubbard_model(H_parameter, Ns, basis, link_list)

    println("Diagonalization...")
    E, tri = @time Krylov.lanczos(H, minite=200, maxite=3000, ϵ=1E-10, nev=10)

    Egs = E[1]
    println("Egs=",Egs)

    println("inverse_iteration...")
    φgs = @time Krylov.inverse_iteration(H, Egs, maxite=2000, ϵ_inv=1E-8, ϵ_cg=1E-5)
    φgs /= norm(φgs)

    H = Fermion.calc_Hubbard_model(H_parameter, Ns, basis, link_list)
    println("compere ",Egs-φgs'*H*φgs)

    println("calc_1d_spectral_func...")
    X, Y, G = Fermion.calc_1d_spectral_func(Egs, φgs, H_parameter, Ns, Ne, pos, unit_vec, basis, link_list)

    plt.figure(figsize=(10,8))
    plt.pcolormesh(X, Y, -1.0/pi*imag(G), cmap="magma")
    plt.colorbar()
    plt.show()

    NΩ = length(G[1,:])
    Ω = Y[1,:]
    Iω = zeros(NΩ)
    for i in 1:NΩ
        for m in 1:Ns
            Iω[i] += -1.0/pi*imag(G[m, i])
        end
    end
    plt.plot(Ω, Iω)
    plt.show()
    
    sns.set_style("white")
    plt.figure(figsize=(10, 8))
    plt.subplot(211)
    @time plt.plot(1:length(E), E)
    plt.subplot(212)
    @time plt.hist(E, bins = 200)

    plt.savefig("dos.png")
    plt.show()
end

test()