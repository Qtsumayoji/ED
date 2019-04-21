module Make_singleband_Hubbard
    include("./Model.jl")

    function main(t1::Float64 = 1.0, t2::Float64 = 0., U::Float64 = 0., V::Float64 = 0., μ::Float64 = 0.)
        filename = "1d_extHubbard.txt"
        # 単位胞あたりのサイト数
        ns = 1
        Nx = 8
        Ny = 1
        Ns = ns*Nx*Ny
        ax = 1.0
        ay = 1.0

        if ARGS != []
            t1 = parse(Float64, ARGS[1])
            t2 = parse(Float64, ARGS[2])
            U  = parse(Float64, ARGS[3])
            V  = parse(Float64, ARGS[4])
            μ  = parse(Float64, ARGS[5]) 
        end

        # 0次近接
        H0 = [U; μ]
        # 1次元鎖の1次近接
        H1 = [t1; V]
        # 1次元鎖の2次近接,正方格子の3次近接
        H3 = [t2; 0.0]

        Model.make_sq_lattice(filename, ns, Nx, Ny, ax, ay, H0, H1,[], H3)

        link_mat, link_list, pos = Model.read_model(filename)
        #println(link_list)
        #show_links(link_mat)
        #println(typeof(link_mat))
        #println(typeof(link_list))
    end
end #end Module

#main()