include("./Model.jl")

function main()
    filename = "corner_sharing_dp.txt"
    # 単位胞あたりのサイト数
    ns = 4
    Nx = 2
    Ny = 1
    Ns = ns*Nx + 1

    ax = 1.0
    ay = 1.0

    Δ = 2.0
    pdσ = -1.5
    pdπ = -pdσ/2.0
    ppσ = 0.5
    ppπ = -0.3*ppσ
    Udd = 8.0
    Upp = Udd/2.0
    Upc = Upp/0.8
    J = 0.0#1.2

    tpp = -0.5*ppπ + 0.5*ppσ
    tdp = sqrt(3.0)/2.0*pdσ

    Hdp = Any[1;tdp]
    Hpp = Any[2;tpp]
    Hp  = Any[3;Upp;0.0]
    Hd  = Any[4;Udd;Δ]
    Hdd = Any[5;J]

    # corner_sharing_unitt_cell
    a = 1.0
    Opx_1 = [0.0; a]
    Opy_2 = [a; 0.0]
    Cux2y2_3 = [a; a]
    Opy_4 = [a; 2a]
    vx = [2a; 0.0]

    fp = open(filename, "w" )

    write(fp, "ns="*string(ns)*"\n")
    write(fp, "Nx="*string(Nx)*"\n")
    write(fp, "Ny="*string(Ny)*"\n")
    write(fp, "Nz="*string(1)*"\n")
    write(fp, "Ns="*string(Ns)*"\n")

    for i in 1:Nx
        write(fp, string((i - 1)*4 + 1)*" "*string(Opx_1[1])*" "*string(Opx_1[2])*" 0\n")
        write(fp, string((i - 1)*4 + 2)*" "*string(Opy_2[1])*" "*string(Opy_2[2])*" 0\n")
        write(fp, string((i - 1)*4 + 3)*" "*string(Cux2y2_3[1])*" "*string(Cux2y2_3[2])*" 0\n")
        write(fp, string((i - 1)*4 + 4)*" "*string(Opy_4[1])*" "*string(Opy_4[2])*" 0\n")
        Opx_1    += vx
        Opy_2    += vx
        Cux2y2_3 += vx
        Opy_4    += vx
        if i == Nx
            write(fp, string(i*4 + 1)*" "*string(Opx_1[1])*" "*string(Opx_1[2])*" 0\n")            
        end
    end

    for i in 1:Nx
        site1 = (i - 1)*4 + 1
        site2 = site1 + 1
        site3 = site2 + 1
        site4 = site3 + 1
        site5 = site4 + 1
        site7 = site5 + 2
        if i == Nx
            site7 = (Nx - 1)*3
        end

        write(fp,string(site1)*" "*string(site1)*" "*string(Hp[1])*" "*string(Hp[2])*" "*string(Hp[3])*"\n")
        write(fp,string(site1)*" "*string(site2)*" "*string(Hpp[1])*" "*string(-Hpp[2])*"\n")
        write(fp,string(site1)*" "*string(site3)*" "*string(Hdp[1])*" "*string(+Hdp[2])*"\n")
        write(fp,string(site1)*" "*string(site4)*" "*string(Hpp[1])*" "*string(+Hpp[2])*"\n")

        write(fp,string(site2)*" "*string(site2)*" "*string(Hp[1])*" "*string(Hp[2])*" "*string(Hp[3])*"\n")
        write(fp,string(site2)*" "*string(site1)*" "*string(Hpp[1])*" "*string(-Hpp[2])*"\n")
        write(fp,string(site2)*" "*string(site3)*" "*string(Hdp[1])*" "*string(-Hdp[2])*"\n")
        write(fp,string(site2)*" "*string(site5)*" "*string(Hpp[1])*" "*string(+Hpp[2])*"\n")

        write(fp,string(site3)*" "*string(site3)*" "*string(Hd[1])*" "*string(Hd[2])*" "*string(Hd[3])*"\n")
        write(fp,string(site3)*" "*string(site1)*" "*string(Hdp[1])*" "*string(+Hdp[2])*"\n")
        write(fp,string(site3)*" "*string(site2)*" "*string(Hdp[1])*" "*string(-Hdp[2])*"\n")
        write(fp,string(site3)*" "*string(site4)*" "*string(Hdp[1])*" "*string(+Hdp[2])*"\n")
        write(fp,string(site3)*" "*string(site5)*" "*string(Hdp[1])*" "*string(-Hdp[2])*"\n")
        write(fp,string(site3)*" "*string(site7)*" "*string(Hdd[1])*" "*string(Hdd[2])*"\n")

        write(fp,string(site4)*" "*string(site4)*" "*string(Hp[1])*" "*string(Hp[2])*" "*string(Hp[3])*"\n")
        write(fp,string(site4)*" "*string(site1)*" "*string(Hpp[1])*" "*string(+Hpp[2])*"\n")
        write(fp,string(site4)*" "*string(site3)*" "*string(Hdp[1])*" "*string(+Hdp[2])*"\n")
        write(fp,string(site4)*" "*string(site5)*" "*string(Hpp[1])*" "*string(-Hpp[2])*"\n")

        write(fp,string(site5)*" "*string(site5)*" "*string(Hp[1])*" "*string(Hp[2])*" "*string(Hp[3])*"\n")
        write(fp,string(site5)*" "*string(site2)*" "*string(Hpp[1])*" "*string(+Hpp[2])*"\n")
        write(fp,string(site5)*" "*string(site3)*" "*string(Hdp[1])*" "*string(-Hdp[2])*"\n")
        write(fp,string(site5)*" "*string(site4)*" "*string(Hpp[1])*" "*string(-Hpp[2])*"\n")
    end

    close(fp)
end

main()