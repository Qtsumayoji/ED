module Fermion
    using LinearAlgebra
    using SparseArrays
    include("./Krylov.jl")

    using PyCall
    @pyimport pylab as plt
    @pyimport seaborn as sns

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
    function make_s_basis(Ns, nup, ndown, basis)
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

    # n↓n↑
    function nund(i::Int64, Ns::Int64, state::Int64)
        if state & (1 << (i - 1)) == 0
            return 0.0
        else
            if state & (1 << (i - 1 + Ns)) == 0
                return 0.0
            else
                return 1.0
            end
        end
    end

    # n
    function niσ(i::Int64, state::Int64)
        if state & (1 << (i - 1)) == 0
            return 0.0
        else
            return 1.0
        end
    end

    function number_op(i::Int64, Ns::Int64, state::Int64)
        return niσ(i, state) + niσ(i + Ns, state)
    end

    function ninj(i::Int64, j::Int64, Ns::Int64, state::Int64)
        return number_op(i, Ns, state)*number_op(j, Ns, state)
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

        sig, ket = ciaj(j, i, state)
        if ket != 0
            push!(row, reverse_basis[state])
            push!(col, reverse_basis[ket])
            push!(val, -t*sig)
        end

        # down spin
        sig, ket = ciaj(j + Ns, i + Ns, state)
        if ket != 0
            push!(row, reverse_basis[state])
            push!(col, reverse_basis[ket])
            push!(val, -t*sig)
        end
    end

    function calc_Hv(Ns::Int64 ,state::Int64, U::Float64, μ::Float64, row::Array{Int64}, col::Array{Int64}, val::Array{Float64}, reverse_basis::Dict)
        diag = 0.0
        for i in 1:Ns
            #n↑n↓
            diag += U*nund(i, Ns, state)

            #n↑ + n↓
            diag -= μ*number_op(i, Ns, state)
        end

        id = reverse_basis[state]
        push!(col, id)
        push!(row, id)
        push!(val, diag)
    end

    function calc_H_diag(Ns::Int64 ,state::Int64, diag::Float64, link_list, row::Array{Int64}, col::Array{Int64}, val::Array{Float64}, reverse_basis::Dict)
        for i in 1:Ns
            link = link_list[i][1]
            para = link[3]
            U = para[1]
            μ = para[2]
            #println(link)

            # n↑n↓
            diag += U*nund(i, Ns, state)
            # n↑ + n↓
            diag -= μ*number_op(i, Ns, state)
        end

        id = reverse_basis[state]
        push!(col, id)
        push!(row, id)
        push!(val, diag)
    end

    function calc_1s3d_for_indirect_RIXS(Ns::Int64 ,state::Int64, jd::Int64, Vd::Float64, row::Array{Int64}, col::Array{Int64}, val::Array{Float64}, reverse_basis::Dict)
        # 1s-3d interaction
        diag = -Vd*number_op(jd, Ns, state)

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

    function calc_H1s3d_for_indirect_RIXS(jd, Vd, system_para, basis)
        println("dim = ", length(basis))

        Ns = system_para.Ns   
        reverse_basis = make_reverse_basis(basis)
    
        row = Int64[]
        col = Int64[]
        val = Float64[]
    
        @inbounds for state in basis
            calc_1s3d_for_indirect_RIXS(Ns, state, jd, Vd, row, col, val, reverse_basis)
        end
    
        H1s3d = sparse(row, col, val)
        return H1s3d
    end

    function calc_extHubbard_model(system_para, basis)
        println("dim = ", length(basis))
        
        Ns = system_para.Ns
        link_list = system_para.link_list

        reverse_basis = make_reverse_basis(basis)
    
        row = Int64[]
        col = Int64[]
        val = Float64[]
    
        @inbounds for state in basis
            diag = 0.0
            @inbounds for i in 1:Ns
                links = link_list[i]
                @inbounds for link in links[2:end]
                    #println(link)
                    j = link[1]
                    #　neighbor = link[2]
                    para = link[3]
                    t = para[1]
                    V = para[2]
                    #println(para)

                    calc_Hkij(i, j, Ns, state, t, row, col, val, reverse_basis)
                    diag += V*ninj(i, j, Ns, state)
                end
            end

            calc_H_diag(Ns, state, diag, link_list, row, col, val, reverse_basis)
        end
    
        H_mat = sparse(row, col, val)
        return H_mat
    end

    function calc_3band_dp_model(system_para, basis)
        println("dim = ", length(basis))
        
        Ns = system_para.Ns
        link_list = system_para.link_list

        reverse_basis = make_reverse_basis(basis)
    
        row = Int64[]
        col = Int64[]
        val = Float64[]
    
        @inbounds for state in basis
            diag = 0.0
            @inbounds for i in 1:Ns
                links = link_list[i]
                @inbounds for link in links[2:end]
                    #println(link)
                    j = link[1]
                    #　neighbor = link[2]
                    para = link[3]
                    t = para[1]
                    #println(para)

                    calc_Hkij(i, j, Ns, state, t, row, col, val, reverse_basis)
                end
            end

            calc_H_diag(Ns, state, diag, link_list, row, col, val, reverse_basis)
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

    function calc_charge_correlation_func(φ, system_para, basis)
        reverse_basis = make_reverse_basis(basis)

        Ns = system_para.Ns
        Nx = system_para.Nx
        Ny = system_para.Ny
        pos = system_para.pos
        unit_vec = system_para.unit_vec
        link_list = system_para.link_list

        CCF = zeros(Ns, Ns)
        
        for i in 1:Ns
            for j in 1:Ns
                x = zeros(length(φ))
                # <ninj>
                @inbounds for state in basis
                    ninj = number_op(i, Ns, state)*number_op(j, Ns, state)
                    ind = reverse_basis[state]
                    x[ind] += ninj*φ[ind]
                end
                CCF[i, j] = φ'*x
            end
        end

        # 逆格子ベクトル
        g = system_para.reciprocal_lattice_vec
        g1 = g[1]
        g2 = g[2]
        g3 = g[3]

        Cq = zeros(Nx, Ny)
        kx = zeros(0)
        ky = zeros(0)
        for m in 1:Nx
            for n in 1:Ny
                q = (m-1)/Nx*g1 + (n-1)/Ny*g2
                qx = q[1]
                qy = q[2]
                push!(kx, qx)
                push!(ky, qy)
    
                tmp_Cq = 0.0
                for i in 1:Ns
                    for j in 1:i - 1
                        Δr = pos[i] - pos[j]
                        tmp_Cq += 2.0*CCF[i, j]*cos.(q'*Δr)
                    end
                    tmp_Cq += CCF[i, i]
                end
                Cq[m, n] = tmp_Cq/Ns
            end
        end

        sns.heatmap(Cq, cmap="magma")
        plt.show()
        return CCF, Cq
    end

    function calc_ci_state(φ, i, Ns, Ne, basis, basis_Np)
        reverse_basis = make_reverse_basis(basis)
        reverse_basis_Nm = make_reverse_basis(basis_Np)
        ϕ = zeros(Complex, length(basis_Np))

        for state in basis
            sign, ket = ci(i, state)
            if ket != 0 
                ϕ[reverse_basis_Nm[ket]] += sign*φ[reverse_basis[state]]
            end
        end

        return ϕ
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
end #end Fermion