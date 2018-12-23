

module Fermion
    using LinearAlgebra
    using SparseArrays
    include("./Krylov.jl").lanczos_vector

    # Ns:格子数
    function make_basis(Ns::Integer)
        nstate = 4^Ns
        basis = Int64[]
        for i in 1:nstate
            push!(basis, i)
        end
        return basis
    end

    # Ns:格子数
    # Ne:電子数
    function make_n_basis(Ns::Integer, Ne::Integer)
        nstate = 4^Ns
        basis = Int64[]
        for i in 1:nstate
            if count_ones(i) == Ne
                push!(basis, i)
            end
        end
        return basis
    end

    # Ns:格子数
    # ndown,nup: down/up spin の数
    function make_s_basis(Ns, ndown, nup, basis)
        sbasis = Int64[]
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

    # fermion signの計算に使用
    bitfrag = Array{Int64}(undef, 32)
    for i in 1:32
        bitfrag[i] = ~((1 << i) - 1)
        #println(bits(bitfrag[i]))
    end

    # 右から1始まりでi番目より左の立っているbitを数える
    function calc_fermion_sign(i::Int64, state::Int64)
        return (-1.0)^(count_ones(state & bitfrag[i]))
    end
    
    function test_calc_fermion_sign()
        test1 = 21343154
        println(bits(test1))
        println(calc_fermion_sign(4, test1)," ",count_ones(test1 & bitfrag[4]))
        test2 =122321
        println(bits(test2))
        println(calc_fermion_sign(7, test2)," ",count_ones(test2 & bitfrag[7]))
    end

    # creation operator
    function ci(i::Int64, state::Int64)
        if (state & (1 << (i - 1))) == 0
            fermion_sign = calc_fermion_sign(i, state)
            state = xor(state , (1 << (i - 1)))
            return fermion_sign, state
        else
            return 0, 0
        end
    end

    # anihiration operator
    function ai(i::Int64, state::Int64)
        if (state & (1 << (i - 1))) == 0
            return 0,0
        else
            fermion_sign = calc_fermion_sign(i, state)
            state = xor(state , (1 << (i - 1)))

            return fermion_sign, state
        end
    end

    # c_i*a_j
    function ciaj(i::Int64, j::Int64, state::Int64)
        if (state & (1 << (j - 1))) == 0
            return 0, 0
        else
            fermion_sign = calc_fermion_sign(j, state)
            state = xor(state , (1 << (j - 1)))
            if (state & (1 << (i - 1))) == 0
                fermion_sign *= calc_fermion_sign(i, state)
                state = xor(state , (1 << (i - 1)))
                return fermion_sign, state
            else
                return 0, 0
            end
        end
    end

    function test_ciaj()
        println(bitstring(21343154))
        sign, state = ciaj(3,5,21343154)
        println(bitstring(state))
        println(sign)
    end

    function calc_Hkij(i::Int64, j::Int64, Ns::Int64, state::Int64, t::Float64, row::Array{Int64}, col::Array{Int64}, val::Array{Float64}, reverse_basis::Dict)
        # up spin
        sig, ket = ciaj(i, j, state)
        if ket != 0
            push!(row, reverse_basis[state])
            push!(col, reverse_basis[ket])
            push!(val, -t*sig)
        end

        # down spin
        sig, ket = ciaj(i + Ns, j + Ns, state)
        if ket != 0
            push!(row, reverse_basis[state])
            push!(col, reverse_basis[ket])
            push!(val, -t*sig)
        end
    end

    # n↓n↑
    function coulmb_repulsion(i::Int64, Ns::Int64, state::Int64)
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

    # n
    function onsite_energy(i::Int64, state::Int64)
        if state & (1 << (i - 1)) == 0
            return 0
        else
            return 1
        end
    end

    function calc_Hv(Ns::Int64 ,state::Int64, U::Float64, μ::Float64, row::Array{Int64}, col::Array{Int64}, val::Array{Float64}, reverse_basis::Dict)
        diag = 0.
        for i in 1:Ns
            #n↑n↓
            ket = coulmb_repulsion(i, Ns, state)
            if ket != 0
                diag += U
            end
            #n↑
            ket = onsite_energy(i, state)
            if ket != 0
                diag -= μ
            end
            #n↓
            ket = onsite_energy(i + Ns, state)
            if ket != 0
                diag -= μ
            end
        end

        id = reverse_basis[state]
        push!(col, id)
        push!(row, id)
        push!(val, diag)
    end

    function calc_Hubbard_model(H_parameter, Ns, basis, link_list)
        println("dim = ", length(basis))
    
        t = H_parameter[1]
        U = H_parameter[2]
        μ = H_parameter[3]

        reverse_basis = make_reverse_basis(basis)
    
        row = Int64[]
        col = Int64[]
        val = Float64[]
    
        for state in basis
            for i in 1:Ns
                link = link_list[i]
                for j in link
                    if i == j
                        continue
                    end
                    calc_Hkij(i, j, Ns, state, t, row, col, val, reverse_basis)
                end
            end
            calc_Hv(Ns, state, U, μ, row, col, val, reverse_basis)
        end
    
        H_mat = sparse(row, col, val)
        return H_mat
    end

    # basisにはあらかじめ粒子数-1の基底を含めておくように
    # Ne粒子系にup spinを一つ加えた時のスペクトル関数を計算
    function calc_ck_state(φ::Array{Float64}, kx::Float64, ky::Float64, Ne::Int64, pos:: Array{Array{Float64}} , basis::Array{Int64}, basis_Np::Array{Int64})
        Ns = length(pos)
        reverse_basis = make_reverse_basis(basis)
        reverse_basis_Np = make_reverse_basis(basis_Np)
        ϕ = zeros(Complex, length(basis_Np))

        # ck = ∑e^(ikri)ci/√Ns
        for i in 1:Ns
            rx = pos[i][1]
            ry = pos[i][2]
            a = exp(-im*kx*rx)
            #println(a'*a)
            for state in basis
                sign, ket = ci(i, state)
                if ket != 0 
                    ϕ[reverse_basis_Np[ket]] += sign*a*φ[reverse_basis[state]]
                end
            end
        end
        
        return ϕ/Ns
    end

    function calc_ak_state(φ::Array{Float64}, kx::Float64, ky::Float64, Ne::Int64, pos::Array{Array{Float64}}, basis::Array{Int64}, basis_Nm::Array{Int64})
        Ns = length(pos)
        reverse_basis = make_reverse_basis(basis)
        reverse_basis_Nm = make_reverse_basis(basis_Nm)
        ϕ = zeros(Complex, length(basis_Nm))

        for i in 1:Ns
            rx = pos[i][1]
            ry = pos[i][2]
            a = exp(im*kx*rx)
            #println(a'*a)
            for state in basis
                sign, ket = ai(i, state)
                if ket != 0 
                    ϕ[reverse_basis_Nm[ket]] += sign*a*φ[reverse_basis[state]]
                end
            end
        end

        return ϕ/Ns
    end

    function calc_ai_state(φ, i, Ns, Ne, basis, basis_Nm)
        reverse_basis = make_reverse_basis(basis)
        reverse_basis_Nm = make_reverse_basis(basis_Nm)
        ϕ = zeros(Complex, length(basis_Nm))

        for state in basis
            sign, ket = ai(i, state)
            if ket != 0 
                ϕ[reverse_basis_Nm[ket]] += sign*φ[reverse_basis[state]]
            end
        end

        return ϕ
    end

    function calc_1d_spectral_func(Egs, φgs, H_parameter, Ns, Ne, pos, unit_vec, basis, link_list)
        η = 1e-1/2
        Nb = 2
        NΩ = 400
        Ω = range(-7.5, stop=7.5, length=NΩ)
        G = zeros(Complex, Nb*Ns, NΩ)
    
        V = dot(unit_vec[1], cross(unit_vec[2], unit_vec[3]))
        #逆格子ベクトル
        g1 = 2*pi*cross(unit_vec[2], unit_vec[3])/V
        g2 = 2*pi*cross(unit_vec[3], unit_vec[1])/V
        g3 = 2*pi*cross(unit_vec[1], unit_vec[2])/V
        
        t = H_parameter[1]
        U = H_parameter[2]
        μ = H_parameter[3]
    
        basis_Np = Fermion.make_n_basis(Ns, Ne+1)
        # 波数表示の生成演算子を作用させるためNe+1の部分空間でのハミルトニアンを作る
        H = Fermion.calc_Hubbard_model(H_parameter, Ns, basis_Np, link_list)
        for m in 1:Nb*Ns
            k = (m - 1)*g1/Ns
            kx = k[1]; ky = k[2]
    
            φex = Fermion.calc_ck_state(φgs, kx, ky, Ne, pos, basis, basis_Np)
            #println(norm(φex))
            norm2_φex = φex'*φex
            φex /= norm(φex)
        
            # 連分数展開の準備
            # 励起状態を試行ベクトルとしてランチョスベクトルを計算
            tridiag = Krylov.lanczos_vector(H, φex; minite = 200)
            tridim = size(tridiag)[1]
            α = []
            for i in 1:tridim
                append!(α, tridiag[i, i])
            end
            β = []
            # βのindexの調整
            append!(β, 0.0)
            for i in 1:tridim - 1
                append!(β, tridiag[i, i+1])
            end
    
            # 各周波数について連分数展開によりG(k,ω)を計算
            for i in 1:NΩ
                ω = Ω[i]
                z = ω + Egs + im*η
                A = calc_continued_fraction_expansion(z, α, β)
                G[m, i] += norm2_φex/A
            end
            println(m*NΩ,"/",NΩ*Nb*Ns)
        end
    
        basis_Np = Fermion.make_n_basis(Ns, Ne-1)
        basis_Nm = basis_Np
        H = Fermion.calc_Hubbard_model(H_parameter, Ns, basis_Nm, link_list)
        @time for m in 1:Nb*Ns
            k = (m - 1)*g1/Ns
            kx = k[1]; ky = k[2]
    
            φex = Fermion.calc_ak_state(φgs, kx, ky, Ne, pos, basis, basis_Nm)
            #println(norm(φex))
            norm2_φex = φex'*φex
            φex /= norm(φex)
        
            # 連分数展開の準備
            # 励起状態を試行ベクトルとしてランチョスベクトルを計算
            tridiag = Krylov.lanczos_vector(H, φex; minite = 200)
            tridim = size(tridiag)[1]
            α = []
            for i in 1:tridim
                append!(α, tridiag[i, i])
            end
            β = []
            # βのindexの調整
            append!(β, 0.0)
            for i in 1:tridim - 1
                append!(β, tridiag[i, i+1])
            end
    
            for i in 1:NΩ
                ω = Ω[i]
                z = ω - Egs + im*η
                A = calc_continued_fraction_expansion(z, -α, -β)
                G[m, i] += norm2_φex/A
            end
            println(m*NΩ,"/",NΩ*Nb*Ns)
        end
    
        # kasikadousukka
        X = zeros(Nb*Ns+1, NΩ)
        Y = zeros(Nb*Ns+1, NΩ)
        for i in 1:NΩ
            for m in 1:Nb*Ns+1
                X[m, i] = (m-1)/Ns*g1[1] - 1/Ns*g1[1]/2
                Y[m, i] = Ω[i]
            end
        end

        return X, Y, G
    end
    
    function calc_continued_fraction_expansion(z::Complex, α, β)
        tridim = length(α)
        A = z - α[tridim]
        for j in tridim:-1:2
            A = z - α[j-1] - β[j]^2.0/A
        end
        
        return A
    end

end