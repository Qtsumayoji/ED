module Model 
    include("./Parameter.jl")
    using LinearAlgebra

    function show_links(link_mat)
        #println(size(link_mat))
        #println(size(link_mat)[1])
        Ns = size(link_mat)[1]
        for i in 1:Ns
            for j in 1:Ns
                if link_mat[i][j] != []
                    print("* ")
                    #print(j," ")
                else
                    print("  ")
                end
            end
            println("(",i," site)")
        end
    end

    function make_link_list(link_mat)
        Ns = size(link_mat)[1]
        link_list = [[] for i in 1:Ns]
        for i in 1:Ns
            tmp = [i ,link_mat[i][i][1], link_mat[i][i][2:end]]
            push!(link_list[i], tmp)
            for j in 1:i
                #println(link_mat[i][j])
                if link_mat[i][j] != [] && i != j
                    tmp = [j ,link_mat[i][j][1], link_mat[i][j][2:end]]
                    push!(link_list[i], tmp)
                end
            end

            println(i," ",link_list[i])
        end

        return link_list
    end

    function read_model(filename)
        fp = open(filename, "r" )
        # 1-3行目はサイト数の読み込み
        ns = readline(fp)
        ns = split(ns, "=")[2]
        ns = parse(Int64, ns)

        Nx = readline(fp)
        Nx = split(Nx, "=")[2]
        Nx = parse(Int64, Nx)

        Ny = readline(fp)
        Ny = split(Ny, "=")[2]
        Ny = parse(Int64, Ny)
    
        Nz = readline(fp)
        Nz = split(Nz, "=")[2]
        Nz = parse(Int64, Nz)

        Ns = readline(fp)
        Ns = split(Ns, "=")[2]
        Ns = parse(Int64, Ns)

        system_size = [ns; Nx; Ny; Nz; Ns]
        #println(system_size)

        pos = Array{Float64}[zeros(3) for i in 1:Ns]
        for i in 1:Ns
            line = readline(fp)
            line = rstrip(line, '\n')
            str = split(line, " ")
            # str[1]はサイトの番号
            str = str[2:end]
            for j in 1:3
                pos[i][j] = parse(Float64, str[j])
            end
            #println(pos[i])
        end

        link_mat = [[[] for i in 1:Ns] for i in 1:Ns]
        for line in eachline(fp)
            if line == ""
                continue
            end
            str = split(line, " ")
            # はじめの2つはsite_i,site_j
            i_site = parse(Int64, str[1])
            j_site = parse(Int64, str[2])
            # 3つめはn次近接を識別するパラメータ
            push!(link_mat[i_site][j_site], parse(Int64, str[3]))
            for a in str[4:end]
                push!(link_mat[i_site][j_site], parse(Float64, a))
            end
            #println(i_site," ",j_site," ",link_mat[i_site][j_site])
        end
        show_links(link_mat)

        link_list = make_link_list(link_mat)

        return link_mat, link_list, pos, system_size
    end

    #そのうち
    function read_RIXS_para(filename)
        core_hole_para = []
        fp = open(filename, "r" )
        # 1-3行目はサイト数の読み込み
        ns = readline(fp)
        ns = split(ns, "=")[2]
        ns = parse(Int64, ns)

        return core_hole_para
    end

    # fileから読み込めるようにする
    function get_square_lattice_unit_vec()
        unit_vec = Array{Float64}[zeros(3) for i in 1:3]
        a = 1.0
        unit_vec[1] = [ a ; 0.0; 0.0]
        unit_vec[2] = [0.0;   a; 0.0]
        unit_vec[3] = [0.0; 0.0;   a]
        return unit_vec
    end

    function push_para(link, i, para)
        if link == []
            push!(link, i)
            push!(link, para)
        end
    end

    function make_sq_lattice(filename, ns, Nx, Ny, ax, ay, H0 = [], H1 = [], H2 = [], H3 = [])
        Ns = Nx*Ny
        link_mat = [[[] for i in 1:Ns] for i in 1:Ns]
        link_list = [[] for i in 1:Ns]
        pos = Array{Float64}[zeros(3) for i in 1:Ns]

        a1 = [ ax; 0.0; 0.0]
        a2 = [0.0;  ay; 0.0]
        a3 = [0.0; 0.0; 1.0]

        fp = open(filename, "w" )
        write(fp, "ns="*string(ns)*"\n")
        write(fp, "Nx="*string(Nx)*"\n")
        write(fp, "Ny="*string(Ny)*"\n")
        write(fp, "Nz="*string(1)*"\n")
        write(fp, "Ns="*string(Ns)*"\n")

        # 1-Nsまでのサイトに振られた番号を計算
        # x正方向に1ずつ増えてy正方向に積み重なっていく
        compute_n(ix, iy, Nx) = ix + Nx*(iy - 1) 
        for iy in 1:Ny
            NU = iy%Ny + 1
            ND = (iy+Ny-2)%Ny + 1
            NUU = (iy+1)%Ny + 1
            NDD = (iy+Ny-3)%Ny + 1

            for ix in 1:Nx
                NR = ix%Nx + 1
                NL = (ix+Nx-2)%Nx + 1
                NRR = (ix+1)%Nx + 1
                NLL = (ix+Nx-3)%Nx + 1
                #println(NLL," ",NL," ",ix," ",NR," ",NRR)

                p   = compute_n( ix,  iy, Nx)
                pL  = compute_n( NL,  iy, Nx)
                pLL = compute_n(NLL,  iy, Nx)
                pR  = compute_n( NR,  iy, Nx)
                pRR = compute_n(NRR,  iy, Nx)
                pU  = compute_n( ix,  NU, Nx)
                pUU = compute_n( ix, NUU, Nx)
                pD  = compute_n( ix,  ND, Nx)
                pDD = compute_n( ix, NDD, Nx)
                pLU = compute_n( NL,  NU, Nx)
                pLD = compute_n( NL,  ND, Nx)
                pRU = compute_n( NR,  NU, Nx)
                pRD = compute_n( NR,  ND, Nx)
                #println(pLL," ",pL," ",p," ",pR," ",pRR)

                if H0 != []
                    push_para(link_mat[p][p], 0, H0)
                end
                # 近接
                if H1 != []
                    push_para(link_mat[p][pL], 1, H1)
                    push_para(link_mat[p][pR], 1, H1)
                    # iy == 1 だと上下方向の周期境界条件でlink_mat[p][p]に値が入る
                    if Ny != 1
                        push_para(link_mat[p][pU], 1, H1)
                        push_para(link_mat[p][pD], 1, H1)
                    end
                end
                # 斜め方向
                if H2 != []
                    push_para(link_mat[p][pLU], 2, H2)
                    push_para(link_mat[p][pLD], 2, H2)
                    push_para(link_mat[p][pRU], 2, H2)
                    push_para(link_mat[p][pRD], 2, H2)
                end
                # 次近接?
                if H3 != []
                    push_para(link_mat[p][pLL], 3, H3)
                    push_para(link_mat[p][pRR], 3, H3)
                    if Ny != 1
                        push_para(link_mat[p][pUU], 3, H3)
                        push_para(link_mat[p][pDD], 3, H3)
                    end
                end

                pos[p] = ix*a1 + iy*a2
                #println(link_mat[p][p])
            end
        end

        for i in 1:Ns
            tmp = string(i)*" "*string(pos[i][1])*" "*string(pos[i][2])*" "*string(pos[i][3])
            write(fp, tmp*"\n")
        end

        for i in 1:Ns
            for j in 1:Ns
                #println(link_mat[i][j])
                if link_mat[i][j] != []
                    write(fp, string(i)*" "*string(j))
                    for linkmat in link_mat[i][j]
                        for para in linkmat
                            write(fp, " "*string(para))
                        end
                    end
                    if i == Ns & j == Ns
                    else
                        write(fp, "\n")
                    end
                end
            end
        end
        close(fp)
    end
end # end Model