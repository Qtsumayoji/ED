module bit_operations
    using StatsBase

    function bit_flip(bit::Int64, i::Int64)
        # 1 << i-1 は左にiずらすのではなく i番目にbitを立てるために1引いている
        return xor(bit , (1 << (i - 1)))
    end

    function bit_on(bit::Int64, i::Int64)
        return bit | (1 << (i - 1)) 
    end

    function bit_off(bit::Int64, i::Int64)
        return bit & ~(1 << (i - 1)) 
    end

    function bit_shift(bit::Int64, i::Int64, Δ::Int64)
        return bit_on(bit_off(bit, i), i + Δ)
    end

    function two_particle_Sz_0_basis(Ns)
        dim = binomial(Ns, 1)*binomial(Ns, 1)
        basis = zeros(Int64, dim)

        cnt = 1
        for i in 1:Ns
            bit = 0
            bit = bit_on(bit, 1)
            bit = bit_on(bit, i + Ns)
            basis[cnt] = bit
            cnt += 1

            for j in 1:Ns-1
                bit = bit_shift(bit, j, 1)
                basis[cnt] = bit
                cnt += 1
            end
        end

        return basis
    end
end