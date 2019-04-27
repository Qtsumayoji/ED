module Spectrum
    using LinearAlgebra
    using SparseArrays
    include("./Krylov.jl")
    include("./Fermion.jl")

    # z - H
    function calc_continued_fraction_expansion(z::Complex, α, β)
        tridim = length(α)
        A = z - α[tridim]
        
        @inbounds for j in tridim:-1:2
            A = z - α[j-1] - β[j]^2.0/A
        end
        
        return A
    end

    # m:Gの波数のindex
    function calc_spectral_func(m, k, Egs, φgs, system_para, spectral_func_para, basis)
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

        basis_Np = Fermion.make_n_basis(Ns, Ne + 1)
        # 波数表示の生成演算子を作用させるためNe+1の部分空間でのハミルトニアンを作る
        #H = Fermion.calc_Hubbard_model(H_para, system_para, basis_Np)
        H = Fermion.calc_extHubbard_model(system_para, basis_Np)
        φex = Fermion.calc_ck_state(φgs, k, pos, basis, basis_Np)

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
    
        basis_Np = Fermion.make_n_basis(Ns, Ne - 1)
        basis_Nm = basis_Np
        #H = Fermion.calc_Hubbard_model(H_para, system_para, basis_Nm)
        H = Fermion.calc_extHubbard_model(system_para, basis_Nm)

        φex = Fermion.calc_ak_state(φgs, k, pos, basis, basis_Nm)
        norm2_φex = φex'*φex
        φex /= norm(φex)
        
        # 連分数展開の準備
        # 励起状態を試行ベクトルとしてランチョスベクトルを計算
        α, β = Krylov.lanczos_vector(H, φex; minite = n_lanczos_vec)
        for i in 1:NΩ
            ω = Ω[i]
            z = Egs - ω -  im*η
            A = -calc_continued_fraction_expansion(z, α, β)
            G[m, i] += norm2_φex/A
        end
    end

    # mは波数だが1 site問題なのでm=1
    function calc_single_site_spectral_func(m, i, Egs, φgs, system_para, spectral_func_para, basis)
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

        basis_Np = Fermion.make_n_basis(Ns, Ne + 1)
        # 波数表示の生成演算子を作用させるためNe+1の部分空間でのハミルトニアンを作る
        #H = Fermion.calc_Hubbard_model(H_para, system_para, basis_Np)
        H = Fermion.calc_extHubbard_model(system_para, basis_Np)
        φex = Fermion.calc_ci_state(φgs, i, Ns, basis, basis_Np)

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
    
        basis_Np = Fermion.make_n_basis(Ns, Ne - 1)
        basis_Nm = basis_Np
        #H = Fermion.calc_Hubbard_model(H_para, system_para, basis_Nm)
        H = Fermion.calc_extHubbard_model(system_para, basis_Nm)

        φex = Fermion.calc_ai_state(φgs, i, Ns, basis, basis_Np)
        norm2_φex = φex'*φex
        φex /= norm(φex)
        
        # 連分数展開の準備
        # 励起状態を試行ベクトルとしてランチョスベクトルを計算
        α, β = Krylov.lanczos_vector(H, φex; minite = n_lanczos_vec)
        for i in 1:NΩ
            ω = Ω[i]
            z = Egs - ω -  im*η
            A = -calc_continued_fraction_expansion(z, α, β)
            G[m, i] += norm2_φex/A
        end
    end

    # m:Gの波数のindex
    function calc_1s4p_RIXS_spectrum(m, q, Egs, φgs, H, system_para, RIXS_para, basis)
        Ns = system_para.Ns 
        Ne = system_para.Ne 
        pos = system_para.pos
        unit_vec = system_para.unit_vec 
        link_list = system_para.link_list 

        Vd = RIXS_para.Vd
        η = RIXS_para.η
        Γ = RIXS_para.Γ
        ωin = RIXS_para.ωin
        n_lanczos_vec = RIXS_para.n_lanczos_vec
        NΩ = RIXS_para.NΩ
        Ω = RIXS_para.Ω 
        G = RIXS_para.G

        dim = length(φgs)
        φex = zeros(Complex, dim)

        z = Diagonal([Egs + ωin + im*Γ for i in 1:dim])
        for id in 1:Ns
            r = pos[id]
            H1s3d = Fermion.calc_intermediate_H(id, Vd, system_para, basis)
            H_RIXS = H + H1s3d
            @time φex += exp(im*q'*r)*Krylov.BiCGSTAB_mod(H_RIXS - real(z), - imag(z), φgs, maxite = 400, ϵ = 1e-3)
        end
        
        norm2_φex = φex'*φex
        φex /= norm(φex)
        α, β = Krylov.lanczos_vector(H, φex; minite = n_lanczos_vec)
        for i in 1:NΩ
            ω = Ω[i]
            z = Egs + ω + im*η
            A = -calc_continued_fraction_expansion(z, α, β)
            # 2 = spin
            G[m, i] = 2.0*norm2_φex/A
        end

        # Δω=0 のRIXS強度の最大値からローレンチアンを作りなおして引き算
        a = 1.0/π*imag(G[m, 1])
        elas = pi*a*η^2.0*[1.0/(Ω[j]^2.0 + η^2.0) for j in 1:NΩ]
        G[m,:] -= im*elas
    end

    # m:Gの波数のindex
    function calc_O2p_RIXS_spectrum(m, Egs, φgs, system_para, RIXS_para, basis)
        Ns = system_para.Ns 
        Ne = system_para.Ne 
        nup = system_para.nup
        ndown = system_para.ndown
        unit_vec = system_para.unit_vec 
        link_list = system_para.link_list 

        Vd = RIXS_para.Vd
        η = RIXS_para.η
        Γ = RIXS_para.Γ
        ωin = RIXS_para.ωin
        n_lanczos_vec = RIXS_para.n_lanczos_vec
        NΩ = RIXS_para.NΩ
        Ω = RIXS_para.Ω 
        G = RIXS_para.G

        Ne = Ne - 1
        nup -= 1
        basis_Nm = Fermion.make_restricted_Hilbert_space(Ns, Ne, nup, ndown)
        dim_Nm = length(basis_Nm)

        Ne = Ne + 1
        nup += 1
        #basis_Np = Fermion.make_restricted_Hilbert_space(Ns, Ne, nup, ndown)
        basis_Np = Fermion.make_n_basis(Ns, Ne)
        dim_Np = length(basis_Np)

        dim = length(basis)

        φai = zeros(dim_Nm)
        φex = zeros(Complex, dim_Nm)
        φci = zeros(dim_Np)

        H = Fermion.calc_3band_dp_model(system_para, basis_Nm)

        z = Diagonal([Egs + ωin + im*Γ for i in 1:dim_Nm])
        for id in 1:Ns
            # Make_3dband_dp_model.jl参照
            # O siteの1体ポテンシャルは でラベルされている
            link = link_list[id][1]
            #println(link)
            if link[2] == 1
                φai = Fermion.calc_ai_state(φgs, id, Ns, basis, basis_Nm)
                H_int = Fermion.calc_intermediate_H(id, -Vd, system_para, basis_Nm)
                H_RIXS = H + H_int
                @time φex = Krylov.BiCGSTAB_mod(real(z) - H_RIXS, imag(z), φai, maxite = 400, ϵ = 1e-3)
                φci += Fermion.calc_ci_state(φex, id, Ns, basis_Nm, basis_Np)
                φci += Fermion.calc_ci_state(φex, id + Ns, Ns, basis_Nm, basis_Np)
            end
        end

        #println("φgs=", φgs)
        #println("φai=", φai)
        #φai = φai/z[1]

        for id in 1:Ns
            link = link_list[id][1]
            if link[2] == 1
            end
        end

        #φai, φex = [], []
        φex = φci

        #println("φex=", φex)
        for i in 1:dim_Np println(i," ",string(basis_Np[i],base=2)," ",φex[i]) end
        for i in 1:dim println(i," ",string(basis[i],base=2)," ",φgs[i]) end

        H = Fermion.calc_3band_dp_model(system_para, basis_Np)

        norm2_φex = φex'*φex
        φex /= norm(φex)
        α, β = Krylov.lanczos_vector(H, φex; minite = n_lanczos_vec)
        for i in 1:NΩ
            ω = Ω[i]
            z = Egs - ω - im*η
            A = calc_continued_fraction_expansion(z, α, β)
            # 2 = spin
            G[m, i] = 2.0*norm2_φex/A
        end        
        
    # down spinを1sから2pへ励起させる過程
#        Ne = Ne + 1
#        ndown += 1
#        basis_Np = Fermion.make_restricted_Hilbert_space(Ns, Ne, nup, ndown)
#        dim_Np = length(basis_Np)
#
#        Ne = Ne - 1
#        ndown -= 1
#        basis_Nm = Fermion.make_restricted_Hilbert_space(Ns, Ne, nup, ndown)
#        dim_Nm = length(basis_Nm)
#
#        φci = zeros(dim_Np)
#        φex = zeros(Complex, dim_Np)
#        φai = zeros(Complex, dim_Nm)
#
#        z = Diagonal([Egs + ωin + im*Γ for i in 1:dim_Np])
#        for id in 1:Ns
#            link = link_list[id][1]
#            if link[2] == 1
#                φci = Fermion.calc_ci_state(φgs, id + Ns, Ns, basis, basis_Np)
#                H = Fermion.calc_3band_dp_model(system_para, basis_Np)
#                H_int = Fermion.calc_intermediate_H(id + Ns, Vd, system_para, basis_Np)
#                H_RIXS = H + H_int
#                @time φex += Krylov.BiCGSTAB_mod(H_RIXS - real(z), - imag(z), φci, maxite = 400, ϵ = 1e-3)
#                φai += Fermion.calc_ai_state(φex, id + Ns, Ns, basis_Np, basis_Nm)
#            end
#        end
#        φci, φex = [], []
#        φex = φai
#        
#        norm2_φex = φex'*φex
#        φex /= norm(φex)
#        H = Fermion.calc_3band_dp_model(system_para, basis_Nm)
#        α, β = Krylov.lanczos_vector(H, φex; minite = n_lanczos_vec)
#        for i in 1:NΩ
#            ω = Ω[i]
#            z = Egs + ω + im*η
#            A = -calc_continued_fraction_expansion(z, α, β)
#            # 2 = spin
#            G[m, i] += norm2_φex/A
#        end

        # Δω=0 のRIXS強度の最大値からローレンチアンを作りなおして引き算
        # G[m, 1]がピークだと仮定している
        a = 1.0/π*imag(G[m, 1])
        elas = pi*a*η^2.0*[1.0/(Ω[j]^2.0 + η^2.0) for j in 1:NΩ]
        G[m,:] -= im*elas
    end

    # m:Gの波数のindex
    function calc_RIXS_spectrum_dense(m, q, Egs, φgs, H, system_para, RIXS_para, basis)
        Ns = system_para.Ns 
        Ne = system_para.Ne 
        pos = system_para.pos
        unit_vec = system_para.unit_vec 
        link_list = system_para.link_list 

        Vd = RIXS_para.Vd
        η = RIXS_para.η
        Γ = RIXS_para.Γ
        ωin = RIXS_para.ωin
        n_lanczos_vec = RIXS_para.n_lanczos_vec
        NΩ = RIXS_para.NΩ
        Ω = RIXS_para.Ω 
        G = RIXS_para.G

        dim = length(φgs)
        φex = zeros(Complex, dim)
        H_RIXS = zeros(Complex, dim, dim)

        z = Diagonal([Egs + ωin + im*Γ for i in 1:dim])
        for id in 1:Ns
            r = pos[id]
            H1s3d = Fermion.calc_H1s3d_for_indirect_RIXS(id, Vd, system_para, basis)
            H_RIXS .= H + H1s3d
            φex += exp(im*q'*r)*inv(H_RIXS - z)*φgs
        end
        
        H_RIXS .= H
        for i in 1:NΩ
            ω = Ω[i]
            z = Egs + ω + im*η
            z = Diagonal([Egs + ω + im*η for i in 1:dim])
            # 2 = spin
            G[m, i] = 2.0*φex'*inv(H_RIXS - z)*φex
        end

        # Δω=0 のRIXS強度の最大値からローレンチアンを作りなおして引き算
        a = 1.0/π*imag(G[m, 1])
        elas = pi*a*η^2.0*[1.0/(Ω[j]^2.0 + η^2.0) for j in 1:NΩ]
        G[m,:] -= im*elas
    end

    # m:Gの波数のindex
    function calc_XAS_spectrum(Egs, φgs, H, system_para, RIXS_para, basis)
        Ns = system_para.Ns 
        Ne = system_para.Ne
        link_list = system_para.link_list 

        Vd = RIXS_para.Vd
        Γ = RIXS_para.Γ
        n_lanczos_vec = RIXS_para.n_lanczos_vec
        NΩ = RIXS_para.NΩ
        Ω = RIXS_para.Ω 
        G = RIXS_para.G

        dim = length(φgs)

        for id in 1:Ns
            H_int = Fermion.calc_intermediate_H(id, Vd, system_para, basis)
            H_XAS = H + H_int
            norm2_φgs = φgs'*φgs
            φgs /= norm(φgs)
            α, β = Krylov.lanczos_vector(H_XAS, φgs; minite = n_lanczos_vec)
            for i in 1:NΩ
                ωin = Ω[i]
                z = Egs + ωin + im*Γ
                A = -calc_continued_fraction_expansion(z, α, β)
                # 2 = spin
                # core-holeの相互作用にスピン依存性がないので単純に2倍する
                G[1, i] = 2.0*norm2_φgs/A
            end
        end
    end

    # m:Gの波数のindex
    function calc_O2p_XAS_spectrum(Egs, φgs, system_para, RIXS_para, basis)
        Ns = system_para.Ns 
        Ne = system_para.Ne
        link_list = system_para.link_list 

        Vd = RIXS_para.Vd
        Γ = RIXS_para.Γ
        n_lanczos_vec = RIXS_para.n_lanczos_vec
        NΩ = RIXS_para.NΩ
        Ω = RIXS_para.Ω 
        G = RIXS_para.G

        basis_Np = Fermion.make_n_basis(Ns, Ne + 1)
        dim = length(basis_Np)
        φci = zeros(dim)
        φex = zeros(Complex, dim)
        dim = length(φgs)

        for id in 1:Ns
            link = link_list[id][1]
            if link[2] == 3
                φci = Fermion.calc_ci_state(φgs, id, Ns, basis, basis_Np)
                H = Fermion.calc_3band_dp_model(system_para, basis_Np)
                H_int = Fermion.calc_intermediate_H(id, Vd, system_para, basis_Np)
                H_XAS = H + H_int
                norm2_φci = φci'*φci
                φci /= norm(φci)
                α, β = Krylov.lanczos_vector(H_XAS, φci; minite = n_lanczos_vec)
                for i in 1:NΩ
                    ωin = Ω[i]
                    z = Egs + ωin + im*Γ
                    A = -calc_continued_fraction_expansion(z, α, β)
                    # 2 = spin
                    # core-holeの相互作用にスピン依存性がないので単純に2倍する
                    G[1, i] = 2.0*norm2_φci/A
                end
            end
        end
    end

    function calc_dynamical_structure_factor(m, q, E, φ, H, system_para, dynamical_structure_factor_para, basis)
        #reverse_basis = Fermion.make_reverse_basis(basis)

        Ns = system_para.Ns
        pos = system_para.pos

        η = dynamical_structure_factor_para.η
        n_lanczos_vec = dynamical_structure_factor_para.n_lanczos_vec
        NΩ = dynamical_structure_factor_para.NΩ
        Ω = dynamical_structure_factor_para.Ω
        G = dynamical_structure_factor_para.G

        dim = length(basis)
        x = zeros(Complex, dim)
        for i in 1:Ns
            # ∑iσ exp(iqr)niσ|x>
            r = pos[i]
            @inbounds for j in 1:dim
                state = basis[j]
                a = exp(im*q'*r)*Fermion.number_op(i, Ns, state)
                #a = Fermion.number_op(i, Ns, state)
                x[j] += a*φ[j]
            end
        end
        #println("<nq>=",φ'*x)
        
        norm2_x = x'*x
        x /= norm(x)

        α, β = Krylov.lanczos_vector(H, x; minite = n_lanczos_vec)
        for i in 1:NΩ
            ω = Ω[i]
            z = E + ω + im*η
            A = calc_continued_fraction_expansion(z, α, β)
            G[m, i] = norm2_x/A
        end
        a = 1.0/π*imag(G[m, 1])
        elas = pi*a*η^2.0*[1.0/(Ω[j]^2.0 + η^2.0) for j in 1:NΩ]
        G[m,:] -= im*elas
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
        Nb = spectral_func_para.Nb
        NΩ = spectral_func_para.NΩ
        Ω = spectral_func_para.Ω
        G = spectral_func_para.G
    
        g = system_para.reciprocal_lattice_vec
        g1 = g[1]
        g2 = g[2]
        g3 = g[3]
        
        t = H_para.t
        U = H_para.U
        μ = H_para.μ
    
        basis_Np = Fermion.make_n_basis(Ns, Ne+1)
        H = Fermion.calc_Hubbard_model(H_para, system_para, basis_Np)
        for m in 1:Nb*Ns
            k = (m - 1)*g1/Ns
    
            φex = Fermion.calc_ck_state(φgs, k, pos, basis, basis_Np)
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
    
            φex = Fermion.calc_ak_state(φgs, k, pos, basis, basis_Nm)
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

end # end Spectrum