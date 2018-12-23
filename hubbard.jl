using LinearAlgebra
using SparseArrays
using Distributions
using PyCall
@pyimport pylab as plt
@pyimport seaborn as sns

include("./model.jl")
include("./Krylov.jl")

#unit cell内のサイト数
const ns = 2
const Nx = 2
const Ny = 2
const Ns = Nx*Ny*ns
const nstate = 4^Ns
const filling = 1
const Nb = 1


function make_n_basis()
    basis = [[] for i in 1:32]
    for i in 1:nstate
        push!(basis[count_ones(i)], i)
    end

    return basis
end

function make_s_basis(ndown, nup, basis)
    sbasis = []
    for base in basis
        if count_ones(base & bitfrag[Ns]) == ndown
            if count_ones(base & ~bitfrag[Ns]) == nup
                push!(sbasis, base)
            end
        end
    end

    return sbasis
end

function make_reverse_basis(basis)
    reverse_basis = Dict()
    dim = length(basis)
    for i in 1:dim
        reverse_basis[basis[i]] = i
        #println(bits(basis[i]))
    end

    return reverse_basis
end

#fermion factorの計算に使用
bitfrag = Array{Int64}(undef, 32)
for i in 1:32
    bitfrag[i] = ~((1 << i) - 1)
    #println(bits(bitfrag[i]))
end

function calc_fermion_factor(i::Int64, state::Int64)
    #return 1.0
    return (-1)^(count_ones(state & bitfrag[i]))
end

function test_calc_fermion_factor()
    test1 = 2134
    println(bitstring(test1))
    i = 4
    println(calc_fermion_factor(i, test1)," ",count_ones(test1 & bitfrag[i]))
    test2 =1223
    println(bitstring(test2))
    i = 2
    println(calc_fermion_factor(i, test2)," ",count_ones(test2 & bitfrag[i]))
end
#test_calc_fermion_factor()

#c_i
function ci(i::Int64, state::Int64)
    if state & (1 << (i - 1)) == 0
        sign = calc_fermion_factor(i, state)
        state = xor(state , (1 << (i - 1)))
        return sign, state
    else
        return 0, 0
    end
end

#c_i*a_j
function ciaj(i::Int64, j::Int64, state::Int64)
    if (state & (1 << (j - 1))) == 0
        return 0,0
    else
        fermion_factor = calc_fermion_factor(j, state)
        state = xor(state , (1 << (j - 1)))
        if (state & (1 << (i - 1))) == 0
            fermion_factor *= calc_fermion_factor(i, state)
            state = xor(state , (1 << (i - 1)))
            return fermion_factor, state
        else
            return 0, 0
        end
    end
end

function test_ciaj()
    println(bin(10))
    sig,state = ciaj(1,2,10)
    println(bin(state))
    println(sig)
    sig,state = ciaj(1+2,2+2,10)
    println(bin(state))
    println(sig)
end
#test_ciaj()

function calc_Hkij(i::Int64, j::Int64, state::Int64, t::float, H_mat::Array{Float64}, reverse_basis::Dict)
    if i == j
        return
    end

    #up spin
    sig, ket = ciaj(i, j, state)
    if ket != 0
        id1 = reverse_basis[state]
        id2 = reverse_basis[ket]
        H_mat[id1, id2] += -t*sig
    end

    #down spin
    sig, ket = ciaj(i + Ns, j + Ns, state)
    if ket != 0
        id1 = reverse_basis[state]
        id2 = reverse_basis[ket]
        H_mat[id1, id2] += -t*sig
    end
end

#n↓n↑
function coulmb_repulsion(i::Int64, state::Int64)
    if state & (1 << (i - 1)) == 0
        return 0
    else
        if state & (1 << (i - 1 + Ns)) == 0
            return 0
        else
            return 1
        end
    end
end

#n
function chemical_potential(i::Int64, state::Int64)
    if state & (1 << (i - 1)) == 0
        return 0
    else
        return 1
    end
end

function calc_Hv(state::Int64, U, μ, H_mat, reverse_basis::Dict)
    id1 = reverse_basis[state]
    for i in 1:Ns
        #n↑n↓
        ket = coulmb_repulsion(i, state)
        if ket != 0
            H_mat[id1, id1] += U
        end
        #n↑
        ket = chemical_potential(i, state)
        if ket != 0
            H_mat[id1, id1] += μ
        end
        #n↓
        ket = chemical_potential(i + Ns, state)
        if ket != 0
            H_mat[id1, id1] += μ
        end
    end
end

function calc_hermiltonian_dense(basis, link_list)
    L = length(basis)
    println("dim = ",L)
    H_mat = zeros(Float64, L, L)
    reverse_basis = make_reverse_basis(basis)

    t = -1.0
    U = 0.0
    μ = -U/4.0

    println("Calculate Hermitian")
    @time for state in basis
        for i in 1:Ns
            link = link_list[i]
            for j in link
                calc_Hkij(i, j, state, t, H_mat, reverse_basis)
            end
        end
        calc_Hv(state, U, μ, H_mat, reverse_basis)
    end

    return Symmetric(H_mat)
end

function calc_excited_state(φ, kx, ky, pos, basis)
    dim = length(basis)
    reverse_basis = make_reverse_basis(basis)
    ϕ = zeros(Complex, dim)

    for state in basis
        for i in 1:Ns
            for j in 1:Ns
                rx = pos[i][1] - pos[j][1]
                ry = pos[i][2] - pos[j][2]
                sig, ket = ciaj(i, j, state)
                if ket != 0
                    c = sig*exp(im*(kx*rx + ky*ry))
                    ϕ[reverse_basis[ket]] += c*φ[reverse_basis[state]]
                end
                sig, ket = ciaj(i + Ns, j + Ns, state)
                if ket != 0
                    c = sig*exp(im*(kx*rx + ky*ry))
                    ϕ[reverse_basis[ket]] += c*φ[reverse_basis[state]]
                end  
            end
        end
    end

    ϕ /= norm(ϕ'*ϕ)
    return ϕ
end

function calc_spectral_function(ω, E, H, φ)
    #println("calc_spectral_function")
    η = 0.1
    dim = length(φ)
    z = (ω + E + im*η)*Matrix(I, dim, dim)

    maxite = 20000
    ϵ = 1E-5

    return (z-H) \ φ
    #return Krylov.CG(z - H, φ, maxite = maxite, ϵ = ϵ)
end

function make_spectral_function(E, φ, H, pos, unit_vec, basis)
    NΩ = 100
    Ω = range(0.0, stop = 10.0, length = NΩ)
    V = dot(unit_vec[1], cross(unit_vec[2], unit_vec[3]))
    #逆格子ベクトル
    g1 = 2*pi*cross(unit_vec[2], unit_vec[3])/V
    g2 = 2*pi*cross(unit_vec[3], unit_vec[1])/V
    g3 = 2*pi*cross(unit_vec[1], unit_vec[2])/V
    Lx = sqrt(dot(unit_vec[1], unit_vec[1]))*Nx
    Ly = sqrt(dot(unit_vec[2], unit_vec[2]))*Ny

    G = zeros(Complex, 2*Nb*Nx, NΩ)
    
    @time for i in 1:NΩ
        for m in -Nb*Nx + 1:Nb*Nx
            n = 0.
            k = m/Lx*g1 + n/Ly*g2
            kx = k[1]
            ky = k[2]
            ω = Ω[i]
            ϕ = calc_excited_state(φ, kx, ky, pos, basis)
            x = calc_spectral_function(ω, E, H, ϕ)
            G[m + Nb*Nx, i] = ϕ'*x
        end
        println((i-1)*(2*Nb*Nx),"/",NΩ*(2*Nb*Nx))
    end

    X = zeros(2*Nb*Nx,NΩ)
    Y = zeros(2*Nb*Nx,NΩ)
    for i in 1:NΩ
        for m in -Nb*Nx + 1:Nb*Nx
            X[m + Nb*Nx,i] = m - 0.5
            Y[m + Nb*Nx,i] = Ω[i]
        end
    end
    plt.pcolormesh(X, Y, -1.0/pi*imag(G), cmap="Blues")
    plt.show()
end

function testHubbardModelOnSquareLattice()
    link_mat, link_list, pos = Model.honeycomb_lattice(Nx, Ny)
    unit_vec = Model.get_square_lattice_unit_vec()
    Model.show_links(link_mat)

    basis = make_n_basis()
    basis = basis[Integer(Ns)]
    Ns_half = Ns/2
    basis = make_s_basis(Ns_half, Ns_half, basis)
    L = length(basis)
    H_mat = calc_hermiltonian_dense(basis, link_list)

    #for i in 1:length(basis)
    #    for j in 1:length(basis)
    #        print(H_mat[i, j]," ")
    #    end
    #    println("|",basis[i],">")
    #end

    E, U = @time eigen(H_mat)
    Egs = E[1]

    φgs = U[:,1]
    println(φgs'*H_mat*φgs)
    #for i in 1:length(basis)
    #    println(E[i])
    #end
    
    sns.set_style("white")
    plt.figure(figsize=(10, 8))
    plt.subplot(311)
    @time plt.hist(E, bins = 300)
    plt.subplot(312)
    @time plt.plot(1:length(E), E)
    #plt.subplot(313)
    #@time plt.plot(E_vec, DOS)
    plt.savefig("dos.png")
    plt.show()
end


function main()
    testHubbardModelOnSquareLattice()
end

main()
#test_calc_fermion_factor()