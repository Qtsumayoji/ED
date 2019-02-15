module Model 
    using LinearAlgebra

    function show_links(link_mat)
        #println(size(link_mat))
        #println(size(link_mat)[1])
        Ns = size(link_mat)[1]
        for i in 1:Ns
            for j in 1:Ns
                if link_mat[i][j] != []
                    print(1," ")
                else
                    print(0," ")
                end
            end
            println("(",i," site)")
        end
    end

    function make_link_list(link_mat)
        Ns = size(link_mat)[1]
        link_list = [[] for i in 1:Ns]
        for i in 1:Ns
            for j in 1:Ns
                #println(link_mat[i][j])
                if link_mat[i][j] != [] && i != j
                    tmp = [j ,link_mat[i][j]]
                    push!(link_list[i], tmp)
                end
            end
        end

        return link_list
    end

    function read_model(filename)
        fp = open(filename, "r" )
        # 1行目はサイト数の読み込み
        Ns = readline(fp)
        Ns = split(Ns, "=")[2]
        Ns = parse(Int64, Ns)
        println("Num site=",Ns)

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
            str = split(line, " ")
            # はじめの2つはsite_i,site_j
            i_site = parse(Int64, str[1])
            j_site = parse(Int64, str[2])

            for a in str[3:end]
                push!(link_mat[i_site][j_site], parse(Float64, a))
            end
            #println(i_site," ",j_site," ",link_mat[i_site][j_site])
        end

        link_list = make_link_list(link_mat)

        return link_mat, link_list, pos
    end

    function push_para(link, para)
        if link == []
            push!(link, para)
        end
    end

    function make_sq_lattice(filename, Nx, Ny, ax, ay, H0 = [], H1 = [], H2 = [], H3 = [])
        Ns = Nx*Ny
        link_mat = [[[] for i in 1:Ns] for i in 1:Ns]
        link_list = [[] for i in 1:Ns]
        pos = Array{Float64}[zeros(3) for i in 1:Ns]

        a1 = [ ax; 0.0; 0.0]
        a2 = [0.0;  ay; 0.0]
        a3 = [0.0; 0.0; 1.0]

        fp = open(filename, "w" )
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
                println(pLL," ",pL," ",p," ",pR," ",pRR)

                if H0 != []
                    push_para(link_mat[p][p], H0)
                end
                # 近接
                if H1 != []
                    push_para(link_mat[p][pL], H1)
                    push_para(link_mat[p][pR], H1)
                    # iy == 1 だと上下方向の周期境界条件でlink_mat[p][p]に値が入る
                    if iy != 1
                        push_para(link_mat[p][pU], H1)
                        push_para(link_mat[p][pD], H1)
                    end
                end
                # 斜め方向
                if H2 != []
                    push_para(link_mat[p][pLU], H2)
                    push_para(link_mat[p][pLD], H2)
                    push_para(link_mat[p][pRU], H2)
                    push_para(link_mat[p][pRD], H2)
                end
                # 次近接?
                if H3 != []
                    push_para(link_mat[p][pLL], H3)
                    push_para(link_mat[p][pRR], H3)
                    if iy != 1
                        push_para(link_mat[p][pUU], H3)
                        push_para(link_mat[p][pDD], H3)
                    end
                end

                pos[p] = ix*a1 + iy*a2
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
                        write(fp,"\n")
                    end
                end
            end
        end
        close(fp)
    end

    function test_make_sq_lattice()
        Nx = 6
        Ny = 1
        ax = 1.0
        ay = 1.0

        t1 = 1.0
        U = 10.0*t1
        V = 2*t1
        t2 = -0.35*t1
        H0 = [0; U]
        H1 = [1; t1; V]
        H3 = [3; t2]

        make_sq_lattice("test.txt", Nx, Ny, ax, ay, H0, H1,[], H3)
    end
    #test_make_sq_lattice()

    function test()
        test_make_sq_lattice()
        link_mat, link_list, pos = read_model("test.txt")
        println(link_list)
        show_links(link_mat)
    end
    test()
end # end Model