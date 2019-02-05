module Krylov
    using LinearAlgebra
    using SparseArrays
    
    function lanczos(A; minite=100, maxite = 2000, ϵ = 1E-5, nev = 1)
        dim = size(A)[1]
        α_st = Float64[]
        β_st = Float64[]
        E0 = [typemax(Float64) for i in 1:nev]
        u0 = rand(dim)
        u0 = u0/norm(u0)
        v = A*u0
        α = u0'*v
        β = norm(v - α*u0)
        u1 = (v - α*u0)/β
        push!(α_st, α)
        push!(β_st, β)

        for ite in 1:maxite
            v = A*u1
            α = u1'*v
            b = β
            β = norm(v - β*u0 - α*u1)
            v = (v - b*u0 - α*u1)/β
            u0 .= u1
            u1 .= v
            
            push!(α_st, α)
            tridiag = Matrix(Tridiagonal(β_st, α_st, β_st))
            E, U = eigen(tridiag)
            push!(β_st, β)

            is_convergence = true
            if nev <= ite
                for i in 1:nev
                    ΔE = norm(E0[i] - E[i])
                    if ΔE > ϵ
                        is_convergence = false
                    end
                end
            else
                is_convergence = false
            end

            if ite < minite
                is_convergence = false
            end

            if is_convergence
                println("lanczos numite = ",ite)
                println("lanczos ΔEgs = ", norm(E0[1] - E[1]))
                return E, tridiag
            end
            E0 = similar(E)
            E0 .= E
        end

        println("*****FAILED*****")
        return 0,0
    end

    function lanczos_vector(A, ϕ; minite = 200)
        dim = size(A)[1]
        tridiag = []
        α_st = Complex[]
        β_st = Complex[]
        u0 = similar(ϕ)
        u0 .= ϕ

        u0 = u0/norm(u0)
        v = A*u0
        α = u0'*v
        β = norm(v - α*u0)
        u1 = (v - α*u0)/β
        push!(α_st, α)
        push!(β_st, β)

        for ite in 1:minite
            v = A*u1
            α = u1'*v
            b = β
            β = norm(v - β*u0 - α*u1)
            v = (v - b*u0 - α*u1)/β
            u0 .= u1
            u1 .= v
            
            push!(α_st, α)
            tridiag = Matrix(Tridiagonal(β_st, α_st, β_st))
            push!(β_st, β)
        end
        return tridiag
    end

    #列ベクトルとしてsiz個のベクトルが入っているuを直交化する
    function gram_schmidt(u::Array)
        siz = size(u)
        n = siz[2]
        for i in 1:n
            for j in 1:i - 1
                a = u[:, j]'*u[:, i]
                u[:, i] += -a*u[:, j]
            end
            u[:, i] /= norm(u[:, i]) 
        end
    end

    #function get_random_vector(dim::Integer, x)
    #    x .= rand(dim)
    #    x /= sqrt(x'*x)
    #end

    function get_random_vector(dim::Integer, x::Array{Float64})
        x .= rand(dim)
    end

    function get_random_vector(dim::Integer, x::Array{Complex{Float64}})
        x .= rand(dim) + im*rand(dim)
    end

    function CG(A, b; maxite = 2000, ϵ = 1E-5)
        println("maxite=",maxite," ϵ=",ϵ)

        x = similar(b)
        r = similar(b)
        p = similar(b)

        x .= b
        r .= b - A*x
        p .= r

        for i in 1:maxite
            y = A*p
            c = r'*r
            α = c/(p'*y)
            x = x + α*p
            r = r - α*y
            β = r'*r/c
            p = r + β*p
            if norm(r) < ϵ
                println("CG numite = ", i)
                println("CG residual = ",norm(r))
                return x
            end
        end

        println("CG residual = ",norm(r))
        println("*****FAILED*****")
        return x
    end

    function test_conv_CG(A, b; maxite = 10000, ϵ = 1E-5)
        println("maxite=",maxite," ϵ=",ϵ)

        res = []

        x = similar(b)
        r = similar(b)
        p = similar(b)
        get_random_vector(length(b), x)
        #x .= b
        r .= b - A*x
        p .= r

        for i in 1:maxite
            y = A*p
            α = r'*r/(p'*y)
            x = x + α*p
            c = r'*r
            r = r - α*y
            β = r'*r/c
            p = r + β*p
            push!(res, norm(r))
            if i%100 == 0
                println("ite=",i," norm(r)=",norm(r))
            end
        end

        return res
    end

    function BiCG(A, b; maxite = 1000, ϵ = 1E-5)
        println("maxite=",maxite," ϵ=",ϵ)
        dim = length(b)
        x = similar(b)
        r = similar(b)
        p = similar(b)
        q = similar(b)
        s = similar(b)
        x *= 0.0
        r .= b
        p .= r
        get_random_vector(dim, q)
        s .= q

        for i in 1:maxite
            y = A*p
            c = s'*r

            α = c/(q'*y)
            x = x + α*p
            r = r - α*y
            s = s - α*A'*q

            β = s'*r/c
            p = r + β*p
            q = s + β*q
            if norm(r) < ϵ
                println("BiCG numite = ", i)
                return x
            end
        end

        println("residual = ",norm(r))
        println("*****FAILED******")
        return x
    end

    #A:matrix
    #σ:approximate eigenvalue
    function inverse_iteration(A, σ; maxite = 1000, ϵ_inv = 1E-9, ϵ_cg = 1e-5)
        dim = size(A)[1]
        for i in 1:dim
            A[i, i] -= σ
        end
        G = A
        φ = zeros(dim)
        get_random_vector(dim, φ)
    
        for ite in 1:maxite
            x = CG(G, φ, maxite = 200, ϵ = ϵ_cg)
            x /= norm(x)

            # ΔE = <x'|A|x'> - σ = σ' - σ
            ΔE = norm(x'*G*x)
            println("InvIte:ΔE = ",ΔE)
            if ΔE < ϵ_inv
                println("inverse_iteration Converged!")
                println("numite = ",ite)
                return x
            end
            φ .= x
        end
        
        return φ
    end
end