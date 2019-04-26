using LinearAlgebra
using SparseArrays

include("./model.jl")
include("./Krylov.jl")
include("./Spin.jl")

#unit cell内のサイト数
const ns = 1
const Nx = 10
const Ny = 2
const Ns = Nx*Ny*ns
const nstate = 2^Ns
const totSz = 0
const Nb = 1

function test()
    link_mat, link_list, pos = Model.square_lattice(Nx, Ny)
    unit_vec = Model.get_square_lattice_unit_vec()
    Model.show_links(link_mat)

    basis = Spin.make_sz_basis(Ns, totSz)

    println("dim=",length(basis))

    # 無理やり開放端境界条件
    #link = [[] for i in 1:Ns]
    #push!(link[1], link_list[1][1])
    #for i in 2:Ns-1
    #    push!(link[i], link_list[i][1])
    #    push!(link[i], link_list[i][2])
    #end 
    #push!(link[Ns], link_list[Ns][2])
    #link_list = link
    
    J = 1.0
    println("make hamiltonian")
    H_mat = @time Spin.calc_Heisenberg_model(J, Ns, basis, link_list)

    println("Diagonalization...")
    E, U = @time Krylov.lanczos(H_mat, maxite=3000, ϵ=1E-15, nev=1)

    Egs = E[1]
    println("Egs/Ns = ",Egs/Ns)
end

test()