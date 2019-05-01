include("./Model.jl")

# 3d軌道の縮退度
deg_3d = 1

function write_inv(fp, site_i, site_j, H)
    write(fp,string(site_i)*" "*string(site_j)*" "*H*"\n")
    write(fp,string(site_j)*" "*string(site_i)*" "*H*"\n")
end

#     p    p    
#   /        \
# p            p
function pp_element_p(fp, site_i, site_j, Hpp)
    Hx_x = Hpp[1]
    Hx_y = Hpp[2]
    Hz_z = Hpp[3]

    # px-px
    write_inv(fp, site_i[1], site_j[1], string(Hx_x[1])*" "*string(Hx_x[2]))
    # px-py
    write_inv(fp, site_i[1], site_j[2], string(Hx_y[1])*" "*string(Hx_y[2]))

    # py-px
    write_inv(fp, site_i[2], site_j[1], string(Hx_y[1])*" "*string(Hx_y[2]))
    # py-py
    write_inv(fp, site_i[2], site_j[2], string(Hx_x[1])*" "*string(Hx_x[2]))

    # pz-pz
    write_inv(fp, site_i[3], site_j[3], string(Hz_z[1])*" "*string(Hz_z[2]))
end

function pp_element_m(fp, site_i, site_j, Hpp)
    Hx_x = Hpp[1]
    Hx_y = Hpp[2]
    Hz_z = Hpp[3]

    # px-px
    write_inv(fp, site_i[1], site_j[1], string(Hx_x[1])*" "*string(Hx_x[2]))
    # px-py
    write_inv(fp, site_i[1], site_j[2], string(Hx_y[1])*" "*string(-Hx_y[2]))

    # py-px
    write_inv(fp, site_i[2], site_j[1], string(Hx_y[1])*" "*string(-Hx_y[2]))
    # py-py
    write_inv(fp, site_i[2], site_j[2], string(Hx_x[1])*" "*string(Hx_x[2]))

    # pz-pz
    write_inv(fp, site_i[3], site_j[3], string(Hz_z[1])*" "*string(Hz_z[2]))
end

function p_element(fp, site_i, Hp)
    Hpx = Hp[1]
    Hpy = Hp[2]
    Hpz = Hp[3]

    # px
    write(fp,string(site_i[1])*" "*string(site_i[1])*" "*string(Hpx[1])*" "*string(Hpx[2])*" "*string(Hpx[3])*"\n")

    # py
    write(fp,string(site_i[2])*" "*string(site_i[2])*" "*string(Hpy[1])*" "*string(Hpy[2])*" "*string(Hpy[3])*"\n")

    # pz
    write(fp,string(site_i[3])*" "*string(site_i[3])*" "*string(Hpz[1])*" "*string(Hpz[2])*" "*string(Hpz[3])*"\n")
end

function dp_element_l(fp, site_Cu, site_O, Hdp)
    Hx_x2y2 = Hdp[1]
    Hx_3z2r2 = Hdp[2]
    Hx_xy = Hdp[3]

    # dx2y2-px
    write_inv(fp, site_Cu[1], site_O[1], string(Hx_x2y2[1])*" "*string(Hx_x2y2[2]))

    # d3z2r2-px
    if deg_3d == 3
        write_inv(fp, site_Cu[2], site_O[1], string(Hx_3z2r2[1])*" "*string(Hx_3z2r2[2]))
        write_inv(fp, site_Cu[3], site_O[2], string(Hx_xy[1])*" "*string(Hx_xy[2]))
    end
end

function dp_element_r(fp, site_Cu, site_O, Hdp)
    Hx_x2y2 = Hdp[1]
    Hx_3z2r2 = Hdp[2]
    Hx_xy = Hdp[3]

    # dx2y2-px
    write_inv(fp, site_Cu[1], site_O[1], string(Hx_x2y2[1])*" "*string(-Hx_x2y2[2]))

    # d3z2r2-px
    if deg_3d == 3
        write_inv(fp, site_Cu[2], site_O[1], string(Hx_3z2r2[1])*" "*string(-Hx_3z2r2[2]))
        write_inv(fp, site_Cu[2], site_O[2], string(Hx_xy[1])*" "*string(-Hx_xy[2]))
    end
end

function dp_element_u(fp, site_Cu, site_O, Hdp)
    Hx_x2y2 = Hdp[1]
    Hx_3z2r2 = Hdp[2]
    Hx_xy = Hdp[3]

    # dx2y2-px
    write_inv(fp, site_Cu[1], site_O[2], string(Hx_x2y2[1])*" "*string(-Hx_x2y2[2]))

    # d3z2r2-px
    if deg_3d == 3
        write_inv(fp, site_Cu[2], site_O[2], string(Hx_3z2r2[1])*" "*string(-Hx_3z2r2[2]))
        write_inv(fp, site_Cu[2], site_O[1], string(Hx_xy[1])*" "*string(-Hx_xy[2]))
    end
end

function dp_element_d(fp, site_Cu, site_O, Hdp)
    Hx_x2y2 = Hdp[1]
    Hx_3z2r2 = Hdp[2]
    Hx_xy = Hdp[3]

    # dx2y2-px
    write_inv(fp, site_Cu[1], site_O[2], string(Hx_x2y2[1])*" "*string(+Hx_x2y2[2]))

    # d3z2r2-px
    if deg_3d == 3
        write_inv(fp, site_Cu[2], site_O[2], string(Hx_3z2r2[1])*" "*string(+Hx_3z2r2[2]))
        write_inv(fp, site_Cu[2], site_O[1], string(Hx_xy[1])*" "*string(+Hx_xy[2]))
    end
end


function d_element(fp, site_Cu, Hd)
    Hdx2y2 = Hd[1]
    Hd3z2r2 = Hd[2]
    Hdxy = Hd[3]

    # dx2y2
    write(fp,string(site_Cu[1])*" "*string(site_Cu[1])*" "*string(Hdx2y2[1])*" "*string(Hdx2y2[2])*" "*string(Hdx2y2[3])*"\n")

    if deg_3d == 3
        # d3z2r2
        write(fp,string(site_Cu[2])*" "*string(site_Cu[2])*" "*string(Hd3z2r2[1])*" "*string(Hd3z2r2[2])*" "*string(Hd3z2r2[3])*"\n")
        write(fp,string(site_Cu[3])*" "*string(site_Cu[3])*" "*string(Hdxy[1])*" "*string(Hdxy[2])*" "*string(Hdxy[3])*"\n")
    end
end

function main()
    filename = "tanioka.txt"


    # 単位胞あたりのサイト数
    ns = 9 + deg_3d
    Nx = 2
    Ny = 1
    Ns = ns*Nx + 3

    ax = 1.0
    ay = 1.0

    #Δ = 2.0
    #pdσ = -1.5
    #pdπ = -pdσ/2.0
    #ppσ = 0.5
    #ppπ = -0.3*ppσ
    #Udd = 8.0
    #Upp = Udd/2.0

    Δ = 2.0
    pdσ = -1.5
    pdπ = 0.7
    ppσ = 0.4
    ppπ = -0.12
    Udd = 8.0
    Upp = Udd/2.0
    
    tx_x =  0.5*ppπ + 0.5*ppσ
    tx_y = -0.5*ppπ + 0.5*ppσ
    tz_z = ppπ

    tx_x2y2  = sqrt(3.0)/2.0*pdσ
    tx_3z2r2 = -0.5*pdσ
    tx_xy    = pdπ

    Hx_x     = Any[11;tx_x]
    Hx_y     = Any[12;tx_y]
    Hz_z     = Any[33;tz_z]

    Hx_x2y2  = Any[14;tx_x2y2]
    Hx_3z2r2 = Any[15;tx_3z2r2]
    Hy_xy    = Any[16;tx_xy]
    
    Hpx      = Any[1;Upp;0.0]
    Hpy      = Any[2;Upp;0.0]
    Hpz      = Any[3;Upp;0.0]

    Hdx2y2   = Any[4;Udd;Δ]
    Hd3z2r2  = Any[5;Udd;0.0]
    Hdxy     = Any[5;Udd;0.0]

    Hpp = [Hx_x, Hx_y, Hz_z]
    Hdp = [Hx_x2y2, Hx_3z2r2, Hy_xy]
    Hp = [Hpx, Hpy, Hpz]
    Hd = [Hdx2y2, Hd3z2r2, Hdxy]
    #Hdd = Any["dd";J]

    # corner_sharing_unitt_cell
    a = 1.0
    O1 = [0.0; a]
    O2 = [a; 0.0]
    Cu3 = [a; a]
    O4 = [a; 2a]
    vx = [2a; 0.0]

    fp = open(filename, "w" )

    write(fp, "ns="*string(ns)*"\n")
    write(fp, "Nx="*string(Nx)*"\n")
    write(fp, "Ny="*string(Ny)*"\n")
    write(fp, "Nz="*string(1)*"\n")
    write(fp, "Ns="*string(Ns)*"\n")

    for i in 1:Nx
        cnt = 0
        write(fp, string((i - 1)*ns + 1)*" "*string(O1[1])*" "*string(O1[2])*" 0\n")
        write(fp, string((i - 1)*ns + 2)*" "*string(O1[1])*" "*string(O1[2])*" 0\n")
        write(fp, string((i - 1)*ns + 3)*" "*string(O1[1])*" "*string(O1[2])*" 0\n")
        write(fp, string((i - 1)*ns + 4)*" "*string(O2[1])*" "*string(O2[2])*" 0\n")
        write(fp, string((i - 1)*ns + 5)*" "*string(O2[1])*" "*string(O2[2])*" 0\n")
        write(fp, string((i - 1)*ns + 6)*" "*string(O2[1])*" "*string(O2[2])*" 0\n")
        cnt = 7
        for j in 1:deg_3d
            write(fp, string((i - 1)*ns + cnt)*" "*string(Cu3[1])*" "*string(Cu3[2])*" 0\n")
            cnt += 1
        end
        write(fp, string((i - 1)*ns + cnt)*" "*string(O4[1])*" "*string(O4[2])*" 0\n")
        write(fp, string((i - 1)*ns + cnt + 1)*" "*string(O4[1])*" "*string(O4[2])*" 0\n")
        write(fp, string((i - 1)*ns + cnt + 2)*" "*string(O4[1])*" "*string(O4[2])*" 0\n")

        O1  += vx
        O2  += vx
        Cu3 += vx
        O4  += vx
        if i == Nx
            write(fp, string(i*ns + 1)*" "*string(O1[1])*" "*string(O1[2])*" 0\n")            
            write(fp, string(i*ns + 2)*" "*string(O1[1])*" "*string(O1[2])*" 0\n")            
            write(fp, string(i*ns + 3)*" "*string(O1[1])*" "*string(O1[2])*" 0\n")            
        end
    end

    for i in 1:Nx
        site1 = (i - 1)*ns + 1
        site_id = [i for i in site1:site1 + ns + 2]
        tmp = 7 + deg_3d - 1
        site1 = site_id[1:3]
        site2 = site_id[4:6]
        site3 = site_id[7:tmp]
        site4 = site_id[tmp + 1:tmp + 3]
        site5 = site_id[tmp + 4:tmp + 6]

        p_element(fp, site1, Hp)
        p_element(fp, site2, Hp)
        d_element(fp, site3, Hd)
        p_element(fp, site4, Hp)
        p_element(fp, site5, Hp)

        pp_element_m(fp, site1, site2, Hpp)
        pp_element_p(fp, site1, site4, Hpp)
        pp_element_p(fp, site2, site5, Hpp)
        pp_element_m(fp, site4, site5, Hpp)

        dp_element_l(fp, site3, site1, Hdp)
        dp_element_d(fp, site3, site2, Hdp)
        dp_element_u(fp, site3, site4, Hdp)
        dp_element_r(fp, site3, site5, Hdp)
    end

    close(fp)
end

main()