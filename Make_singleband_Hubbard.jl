include("./Model.jl")

function main()
    filename = "1d_extHubbard.txt"
    # 単位胞あたりのサイト数
    ns = 1
    Nx = 6
    Ny = 1
    Ns = ns*Nx*Ny
    ax = 1.0
    ay = 1.0
    
    t1 = 1.0
    U  = 10.0*t1
    V  = 1.5*t1
    t2 = -0.34*t1
    μ  = 0.0#U/2
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

main()