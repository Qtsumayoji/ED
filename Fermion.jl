module Fermion
    using LinearAlgebra
    using SparseArrays
    using Base.Threads
    include("./Krylov.jl")

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
        diag = 0.0
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

    function calc_Hv_for_indirect_RIXS(Ns::Int64 ,state::Int64, U::Float64, μ::Float64, jd::Int64, Vd::Float64, row::Array{Int64}, col::Array{Int64}, val::Array{Float64}, reverse_basis::Dict)
        diag = 0.0
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

        # 1s-3d interaction
        ket = onsite_energy(j, state)
        if ket != 0
            diag -= Vd
        end
        ket = onsite_energy(j + Ns, state)
        if ket != 0
            diag -= Vd
        end

        id = reverse_basis[state]
        push!(col, id)
        push!(row, id)
        push!(val, diag)
    end

    function calc_Hubbard_model(H_para, system_para, basis)
        println("dim = ", length(basis))
        
        Ns = system_para.Ns
        link_list = system_para.link_list
        
        t = H_para.t
        U = H_para.U
        μ = H_para.μ

        reverse_basis = make_reverse_basis(basis)
    
        row = Int64[]
        col = Int64[]
        val = Float64[]
    
        @inbounds for state in basis
            @inbounds for i in 1:Ns
                link = link_list[i]
                @inbounds for j in link
                    calc_Hkij(i, j, Ns, state, t, row, col, val, reverse_basis)
                end
            end
            calc_Hv(Ns, state, U, μ, row, col, val, reverse_basis)
        end
    
        H_mat = sparse(row, col, val)
        return H_mat
    end

    function calc_Hubbard_model_for_indirect_RIXS(H_para, jd, Vd, system_para, basis)
        println("dim = ", length(basis))

        Ns = system_para.Ns
        link_list = system_para.link_list
    
        t = H_para[1]
        U = H_para[2]
        μ = H_para[3]

        reverse_basis = make_reverse_basis(basis)
    
        row = Int64[]
        col = Int64[]
        val = Float64[]
    
        @inbounds for state in basis
            @inbounds for i in 1:Ns
                link = link_list[i]
                @inbounds for j in link
                    calc_Hkij(i, j, Ns, state, t, row, col, val, reverse_basis)
                end
            end
            calc_Hv_for_indirect_RIXS(Ns, state, U, μ, row, col, val, jd, Vd, reverse_basis)
        end
    
        H_mat = sparse(row, col, val)
        return H_mat
    end

    # basisにはあらかじめ粒子数-1の基底を含めておくように
    # Ne粒子系にup spinを一つ加えた時のスペクトル関数を計算
    function calc_ck_state(φ::Array{Float64}, k, Ne::Int64, pos:: Array{Array{Float64}} , basis::Array{Int64}, basis_Np::Array{Int64})
        Ns = length(pos)
        reverse_basis = make_reverse_basis(basis)
        reverse_basis_Np = make_reverse_basis(basis_Np)
        ϕ = zeros(Complex, length(basis_Np))
        kx = k[1]
        ky = k[2]

        # ck = ∑e^(ikri)ci/√Ns
        @inbounds for i in 1:Ns
            r = pos[i]
            a = exp(-im*k'*r)
            #println(a'*a)
            @inbounds for state in basis
                sign, ket = ci(i, state)
                if ket != 0 
                    ϕ[reverse_basis_Np[ket]] += sign*a*φ[reverse_basis[state]]
                end
            end
        end
        
        return ϕ/Ns
    end

    function calc_ak_state(φ, k, Ne, pos, basis, basis_Nm)
        Ns = length(pos)
        reverse_basis = make_reverse_basis(basis)
        reverse_basis_Nm = make_reverse_basis(basis_Nm)
        ϕ = zeros(Complex, length(basis_Nm))
        kx = k[1]
        ky = k[2]

        @inbounds for i in 1:Ns
            r = pos[i]
            a = exp(im*k'*r)
            #println(a'*a)
            @inbounds for state in basis
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

    function calc_continued_fraction_expansion(z::Complex, α, β)
        tridim = length(α)
        A = z - α[tridim]
        
        for j in tridim:-1:2
            A = z - α[j-1] - β[j]^2.0/A
        end
        
        return A
    end

    function calc_1d_spectral_func(Egs, φgs, H_para, system_para, spectral_func_para, basis)
        ns = system_para.ns
        Nx = system_para.Nx
        Ny = system_para.Ny
        Ns = system_para.Ns
        Ne = system_para.Ne
        pos = system_para.pos
        unit_vec = system_para.unit_vec
        link_list = system_para.link_list

        η = spectral_func_para.η
        n_lanczos_vec = spectral_func_para.n_lanczos_vec
        Nb = 2
        NΩ = spectral_func_para.NΩ
        Ω = spectral_func_para.Ω
        G = spectral_func_para.G
    
        #逆格子ベクトル
        g = system_para.reciprocal_lattice_vec
        g1 = g[1]
        g2 = g[2]
        g3 = g[3]
        
        t = H_para.t
        U = H_para.U
        μ = H_para.μ
    
        basis_Np = Fermion.make_n_basis(Ns, Ne+1)
        # 波数表示の生成演算子を作用させるためNe+1の部分空間でのハミルトニアンを作る
        H = Fermion.calc_Hubbard_model(H_para, system_para, basis_Np)
        for m in 1:Nb*Ns
            k = (m - 1)*g1/Ns
    
            φex = Fermion.calc_ck_state(φgs, k, Ne, pos, basis, basis_Np)
            #println(norm(φex))
            norm2_φex = φex'*φex
            φex /= norm(φex)
        
            # 連分数展開の準備
            # 励起状態を試行ベクトルとしてランチョスベクトルを計算
            α, β = Krylov.lanczos_vector(H, φex; minite = n_lanczos_vec)
    
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
        H = Fermion.calc_Hubbard_model(H_para, system_para, basis_Nm)
        for m in 1:Nb*Ns
            k = (m - 1)*g1/Ns
            kx = k[1]; ky = k[2]
    
            φex = Fermion.calc_ak_state(φgs, k, Ne, pos, basis, basis_Nm)
            #println(norm(φex))
            norm2_φex = φex'*φex
            φex /= norm(φex)
        
            # 連分数展開の準備
            # 励起状態を試行ベクトルとしてランチョスベクトルを計算
            α, β = Krylov.lanczos_vector(H, φex; minite = n_lanczos_vec)
    
            for i in 1:NΩ
                ω = Ω[i]
                z = ω - Egs + im*η
                A = calc_continued_fraction_expansion(z, -α, -β)
                G[m, i] += norm2_φex/A
            end
            println(m*NΩ,"/",NΩ*Nb*Ns)
        end
    end

    # m:Gの波数のindex
    function calc_spectral_func_k(m, k, Egs, φgs, H_para, system_para, spectral_func_para, basis)
        ns = system_para.ns
        Nx = system_para.Nx
        Ny = system_para.Ny
        Ns = system_para.Ns
        Ne = system_para.Ne
        pos = system_para.pos
        unit_vec = system_para.unit_vec
        link_list = system_para.link_list

        η = spectral_func_para.η
        n_lanczos_vec = spectral_func_para.n_lanczos_vec
        NΩ = spectral_func_para.NΩ
        Ω = spectral_func_para.Ω
        G = spectral_func_para.G
        
        t = H_para.t
        U = H_para.U
        μ = H_para.μ
    
        basis_Np = Fermion.make_n_basis(Ns, Ne+1)
        # 波数表示の生成演算子を作用させるためNe+1の部分空間でのハミルトニアンを作る
        H = Fermion.calc_Hubbard_model(H_para, system_para, basis_Np)
        φex = Fermion.calc_ck_state(φgs, k, Ne, pos, basis, basis_Np)

        norm2_φex = φex'*φex
        φex /= norm(φex)
    
        # 連分数展開の準備
        # 励起状態を試行ベクトルとしてランチョスベクトルを計算
        α, β = Krylov.lanczos_vector(H, φex; minite =n_lanczos_vec)
    
        # 各周波数について連分数展開によりG(k,ω)を計算
        for i in 1:NΩ
            ω = Ω[i]
            z = ω + Egs + im*η
            A = calc_continued_fraction_expansion(z, α, β)
            G[m, i] += norm2_φex/A
        end
    
        basis_Np = Fermion.make_n_basis(Ns, Ne-1)
        basis_Nm = basis_Np
        H = Fermion.calc_Hubbard_model(H_para, system_para, basis_Nm)

        φex = Fermion.calc_ak_state(φgs, k, Ne, pos, basis, basis_Nm)
        #println(norm(φex))
        norm2_φex = φex'*φex
        φex /= norm(φex)
    
        # 連分数展開の準備
        # 励起状態を試行ベクトルとしてランチョスベクトルを計算
        α, β = Krylov.lanczos_vector(H, φex; minite = n_lanczos_vec)

        for i in 1:NΩ
            ω = Ω[i]
            z = ω - Egs + im*η
            A = calc_continued_fraction_expansion(z, -α, -β)
            G[m, i] += norm2_φex/A
        end
    end

    function calc_RIXS_spectra(Egs, φgs, H_parameter, system_para, basis)
        η = 1e-2
        n_lanczos_vec = 200
        Nb = 1
        # 入射光の振動数
        NΩ = 1000
        Ω = range(-3., stop=3., length=NΩ)
        # x-ray の運動量変化
        NQ = 10
        Q = range(0.0, stop=pi, lingth=NQ)

        G = zeros(Complex, NQ, NΩ)
        
        t = H_parameter[1]
        U = H_parameter[2]
        μ = H_parameter[3]
        ns = system_para.ns
        Nx = system_para.Nx
        Ny = system_para.Ny
        Ns = system_para.Ns
        Ne = system_para.Ne
        pos = system_para.pos
        unit_vec = system_para.unit_vec
        link_list = system_para.link_list
        H = Fermion.calc_Hubbard_model(H_parameter, Ns, basis, link_list)
    end

end