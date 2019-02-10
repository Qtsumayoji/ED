module Spectrum
    using LinearAlgebra
    using SparseArrays
    include("./Krylov.jl")
    include("./Fermion.jl")

    function calc_continued_fraction_expansion(z::Complex, α, β)
        tridim = length(α)
        A = z - α[tridim]
        
        for j in tridim:-1:2
            A = z - α[j-1] - β[j]^2.0/A
        end
        
        return A
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
            z = Egs - ω -  im*η
            A = -calc_continued_fraction_expansion(z, α, β)
            G[m, i] += norm2_φex/A
        end
    end

    # m:Gの波数のindex
    function calc_RIXS_spectrum(m, q, Egs, φgs, H, system_para, RIXS_para, basis)
        Nx = system_para.Nx
        Ny = system_para.Ny
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

        φex = similar(φgs)

        for id in 1:Ns
            r = pos[id]
            H1s3d = Fermion.calc_H1s3d_for_indirect_RIXS(id, Vd, system_para, basis)
            H_RIXS = H + H1s3d
            dim = length(φgs)
            z = Diagonal([Egs + ωin + im*Γ for i in 1:dim])
            φex += exp(im*q'*r)*Krylov.CG(H_RIXS - z, φgs, maxite = 200)
        end
        
        # spin
        φex *= 2.0
        norm2_φex = φex'*φex
        φex /= norm(φex)
        α, β = Krylov.lanczos_vector(H, φex; minite = n_lanczos_vec)
        for i in 1:NΩ
            ω = Ω[i]
            z = Egs + ω + im*η
            A = calc_continued_fraction_expansion(z, -α, -β)
            G[m, i] += norm2_φex/A
        end
    end

    function calc_dynamical_structure_factor(m, q, E, φ, H, system_para, dynamical_structure_factor_para, basis)
        reverse_basis = Fermion.make_reverse_basis(basis)

        Ns = system_para.Ns
        Nx = system_para.Nx
        Ny = system_para.Ny
        Ne = system_para.Ne
        pos = system_para.pos

        η = dynamical_structure_factor_para.η
        n_lanczos_vec = dynamical_structure_factor_para.n_lanczos_vec
        NΩ = dynamical_structure_factor_para.NΩ
        Ω = dynamical_structure_factor_para.Ω
        G = dynamical_structure_factor_para.G

        x = zeros(Complex, length(φ))
        for i in 1:Ns
            # ∑iσ exp(iqr)niσ|x>
            r = pos[i]
            j = 1
            @inbounds for state in basis
                a = exp(im*q'*r)*Fermion.number_op(i, Ns, state)
                x[j] += a*φ[j]
                j += 1
            end
        end

        norm2_x = x'*x
        x /= norm(x)

        α, β = Krylov.lanczos_vector(H, x; minite = n_lanczos_vec)

        for i in 1:NΩ
            ω = Ω[i]
            z = E + ω + im*η
            A = calc_continued_fraction_expansion(z, α, β)
            G[m, i] += norm2_x/A
        end

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

end # end Spectrum