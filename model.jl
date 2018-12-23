module Model
    using LinearAlgebra
    const no_link = 0

    function get_square_lattice_unit_vec()
        unit_vec = Array{Float64}[zeros(3) for i in 1:3]
        a = 1.0
        unit_vec[1] = [ a ; 0.0; 0.0]
        unit_vec[2] = [0.0;   a; 0.0]
        unit_vec[3] = [0.0; 0.0;   a]
        return unit_vec
    end

    function square_lattice(width::Int, height::Int)
        Ns = width*height
        link_mat = zeros(Float64, Ns, Ns)
        link_list = [[] for i in 1:Ns]
        pos = Array{Float64}[zeros(3) for i in 1:Ns]

        a = 1.0
        a1 = [  a; 0.0; 0.0]
        a2 = [0.0;   a; 0.0]
        a3 = [0.0; 0.0;   a]

        for i in 1:height
            for j in 1:width
                site = (i - 1)*width + j
                r_site = site + 1
                l_site = site - 1
                u_site = site + width
                d_site = site - width

                #周期境界条件
                if j == width
                    r_site = site - width + 1
                end
                if j == 1
                    l_site = site + width - 1
                end
                if i == height
                    u_site = j
                end
                if i == 1
                    d_site = (height - 1)*width + j
                end
                #println(site," ",r_site," ",l_site," ",u_site," ",d_site)

                link_mat[site, r_site] = 1
                link_mat[site, l_site] = 1
                link_mat[site, u_site] = 1
                link_mat[site, d_site] = 1

                pos[site] = i*a2 + j*a1 
            end
        end
        
        if height == 1
            for i in 1:Ns
                link_mat[i, i] = no_link
            end
        end

        make_link_list(link_mat, link_list, Ns)

        return link_mat, link_list, pos
    end

    # 未完成
    function triangular_lattice(links, height, width)
        Ns = height*width*2
        a1 = [1., 0., 0.]
        a2 = [cos(π/3), sin(π/3), 0.]
        a3 = [0., 0., 1.]
    
        for i in 1:height
            for j in 1:width
                site = (i-1)*width + j
    
                tunitcell_l = (i-1)*width + j - 1
                #左端
                if j == 1
                    tunitcell_l = (i-1)*width + width
                end
    
                tunitcell_r = (i-1)*width + j + 1
                #右端
                if j == width
                    tunitcell_r = (i-1)*width + 1
                end
    
                #上端
                tunitcell_d1 = (i-2)*width + j
                if i == 1
                    tunitcell_d1 = (height-1)*width + j
                end
    
                #上端
                tunitcell_d2 = (i-2)*width + j + 1
                if i == 1
                    tunitcell_d2 = (height-1)*width + j + 1
                end
    
                tunitcell_u1 = i*width + j
                #下端
                if i == height
                    tunitcell_u1 = j
                end
    
                tunitcell_u2 = i*width + j + 1
                #下端
                if i == height
                    tunitcell_u2 = j + 1
                end
    
                links[site, tunitcell_l] = 1
                links[site, tunitcell_r] = 1
                links[site, tunitcell_u1] = 1
                links[site, tunitcell_d1] = 1
                links[site, tunitcell_u2] = 1
                links[site, tunitcell_d2] = 1
    
            end
        end
    end

    function honeycomb_lattice(width, height)
        Ns = 2*width*height
        link_mat = zeros(Float64, Ns, Ns)
        link_list = [[] for i in 1:Ns]
        pos = Array{Float64}[zeros(3) for i in 1:Ns]

        a1 = [1.0, 0.0, 0.0]
        a2 = [1.0/2.0, sqrt(3.0)/2.0, 0.0]
        a3 = [0.0, 0.0, 1.0]
    
        a21 = a2 + a1
        unitcell1 = a1 + a21/3.0
        unitcell2 = a1 + 2.0/3.0*a21
        for i in 1:height
            for j in 1:width
                tunitcell = (i-1)*width + j
                tsite1 = 2*tunitcell - 1
                tsite2 = tsite1 + 1
    
                tunitcell_l = tunitcell - 1
                #左端
                if j == 1
                    tunitcell_l = (i-1)*width + width
                end
    
                tunitcell_r = tunitcell + 1
                #右端
                if j == width
                    tunitcell_r = (i-1)*width + 1
                end
    
                #上端
                tunitcell_d = (i - 2)*width + j
                if i == 1
                    tunitcell_d = (height - 1)*width + j
                end
    
                tunitcell_u = i*width + j
                #下端
                if i == height
                    tunitcell_u = j
                end
    
                link_mat[tsite1, tsite2] = 1
                link_mat[tsite1, 2*tunitcell_l] = 1
                link_mat[tsite1, 2*tunitcell_d] = 1
                link_mat[tsite2, tsite1] = 1
                link_mat[tsite2, 2*tunitcell_r - 1] = 1
                link_mat[tsite2, 2*tunitcell_u - 1] = 1
            end
        end

        make_link_list(link_mat, link_list, Ns)

        return link_mat, link_list, pos
    end

    function show_links(link_mat)
        Ns = length(link_mat[1,:])
        for i in 1:Ns
            for j in 1:Ns
                if link_mat[i, j] != no_link
                    print(1," ")
                else
                    print(0," ")
                end
            end
            println("(",i," site)")
        end
    end

    function make_link_list(link_mat, link_list, Ns)
        for i in 1:Ns
            for j in 1:Ns
                if link_mat[i, j] != no_link & i != j
                    push!(link_list[i], j)
                end
            end
        end
    end
end