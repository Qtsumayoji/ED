module Spin
    using LinearAlgebra
    using SparseArrays

    # sz=constとなる基底を作成
    function make_sz_basis(Ns::Int64, totSz::Int64)
        nstate = 2^Ns
        basis = []
        for i in 1:nstate
            pSz = count_ones(i)
            mSz = Ns - count_ones(i)
            Sz = pSz - mSz
            if Sz == totSz
                push!(basis, i)
                #println(bitstring(i)," ",Sz)
            end
        end
        
        return basis
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

    # s+s-
    function spin_pm(i::Int64, j::Int64, state::Int64)
        if (state & (1 << (j - 1))) == 0
            return 0
        else
            state = xor(state , (1 << (j - 1)))
            if (state & (1 << (i - 1))) == 0
                state = xor(state , (1 << (i - 1)))
                return state
            else
                return 0
            end
        end
    end

    # s-_i s+_j
    function spin_mp(i::Int64, j::Int64, state::Int64)
        # j番目のbitが1ならばreturn 0
        if (state & (1 << (j - 1))) == 0
            # j番目のbitを反転させる
            state = xor(state , (1 << (j - 1)))
            if (state & (1 << (i - 1))) == 0
                return 0
            else
                state = xor(state , (1 << (i - 1)))
                # ここまで通れば行列要素は非ゼロ
                return state
            end
        else
            return 0
        end
    end

    function spin_zz(i::Int64, j::Int64, state::Int64)
        if (state & (1 << (j - 1))) == 0
            if (state & (1 << (i - 1))) == 0
                return 0.25
            else
                return -0.25
            end
        else
            if (state & (1 << (i - 1))) == 0
                return -0.25
            else
                return 0.25
            end
        end
    end

    function calc_sxsy(i::Int64, j::Int64, state::Int64, J::Float64, row, col, val, reverse_basis::Dict)
        ket = spin_pm(i, j, state)
        if ket != 0
            push!(row, reverse_basis[state])
            push!(col, reverse_basis[ket])
            push!(val, -J*0.5)
        end

        ket = spin_mp(i, j, state)
        if ket != 0
            push!(row, reverse_basis[state])
            push!(col, reverse_basis[ket])
            push!(val, -J*0.5)
        end
    end
    
    # coo形式の疎行列を作成
    # link_list[i]:i番目のサイトとつながっているサイトの番号が格納されている
    # 例:1d chainのとき3番目のサイトとつながっているものは link_list[3] = [2,5]
    # から取得できる
    # 
    # J         :Heisenberg相互作用の係数,J>0で強磁性を想定
    # Ns        :総サイト数
    # basis     :szの保存した基底
    # link_list :上記
    function calc_Heisenberg_model(J, Ns, basis, link_list)
        reverse_basis = make_reverse_basis(basis)
    
        # 行、列のindexとvalueを格納
        row = Int64[]
        col = Int64[]
        val = Float64[]
    
        for state in basis
            diag = 0.0

            # 対角要素の計算と非対角要素をrow,col,valに格納していく
            for i in 1:Ns
                link = link_list[i]
                for j in link
                    calc_sxsy(i, j, state, J, row, col, val, reverse_basis)
                    diag += -J*spin_zz(i, j, state)
                end
            end

            id = reverse_basis[state]
            push!(col, id)
            push!(row, id)
            push!(val, diag)
        end

        H_mat = sparse(row, col, val)
        return H_mat
    end
end