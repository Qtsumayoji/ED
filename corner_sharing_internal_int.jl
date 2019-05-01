include("./Model.jl")

function write_inv(fp, site_i, site_j, H)
    write(fp,string(site_i)*" "*string(site_j)*" ")
    write_H(fp, H)

    write(fp,string(site_j)*" "*string(site_i)*" ")
    write_H(fp, H)
end

function write_H(fp, H)
    for h in H
        write(fp,string(h)*" ")
    end
    write(fp, "\n")
end    

#     p    p    
#   /        \
# p            p
function pp_element(fp, α, site_i, site_j, Hpp)
    Hx_x = Hpp[1]
    Hx_y = Hpp[2]
    Hx_y_α = Any[Hx_y[1], α*Hx_y[2], Hx_y[3]]
    Hz_z = Hpp[3]

    # px-px
    write_inv(fp, site_i[1], site_j[1], Hx_x)
    # px-py
    write_inv(fp, site_i[1], site_j[2], Hx_y_α)

    # py-px
    write_inv(fp, site_i[2], site_j[1], Hx_y_α)
    # py-py
    write_inv(fp, site_i[2], site_j[2], Hx_x)

    # pz-pz
    write_inv(fp, site_i[3], site_j[3], Hz_z)
end

function p_element(fp, site_i, Hp)
    Hpx = Hp[1]
    Hpy = Hp[2]
    Hpz = Hp[3]

    # px
    write(fp,string(site_i[1])*" "*string(site_i[1])*" ")
    write_H(fp, Hpx)

    # py
    write(fp,string(site_i[2])*" "*string(site_i[2])*" ")
    write_H(fp, Hpy)

    # pz
    write(fp,string(site_i[3])*" "*string(site_i[3])*" ")
    write_H(fp, Hpz)
end

function dp_element_h(fp, α, site_Cu, site_O, Hdp)
    Hx_x2y2 = Hdp[1]
    Hx_x2y2_α = Any[Hx_x2y2[1], α*Hx_x2y2[2], Hx_x2y2[3]]
    Hx_3z2r2 = Hdp[2]
    Hx_3z2r2_α = Any[Hx_3z2r2[1], α*Hx_3z2r2[2], Hx_3z2r2[3]]
    Hz_3z2r2 = Hdp[3]

    # dx2y2-px
    write_inv(fp, site_Cu[1], site_O[1], Hx_x2y2_α)

    # d3z2r2-px
    write_inv(fp, site_Cu[2], site_O[1], Hx_3z2r2_α)

    # d3z2r2-pz
    #write_inv(fp, site_Cu[2], site_O[3], string(Hz_3z2r2[1])*" "*string(Hz_3z2r2[2]))
end

function dp_element_v(fp, α, site_Cu, site_O, Hdp)
    Hx_x2y2 = Hdp[1]
    Hx_x2y2_α = Any[Hx_x2y2[1], α*Hx_x2y2[2], Hx_x2y2[3]]
    Hx_3z2r2 = Hdp[2]
    Hx_3z2r2_α = Any[Hx_3z2r2[1], α*Hx_3z2r2[2], Hx_3z2r2[3]]
    Hz_3z2r2 = Hdp[3]

    # dx2y2-py
    write_inv(fp, site_Cu[1], site_O[2], Hx_x2y2_α)

    # d3z2r2-py
    write_inv(fp, site_Cu[2], site_O[2], Hx_3z2r2_α)

    # d3z2r2-pz
    #write_inv(fp, site_Cu[2], site_O[3], string(Hz_3z2r2[1])*" "*string(Hz_3z2r2[2]))
end

function d_element(fp, site_Cu, Hd)
    Hdx2y2 = Hd[1]
    Hd3z2r2 = Hd[2]

    # dx2y2
    write(fp,string(site_Cu[1])*" "*string(site_Cu[1])*" ")
    write_H(fp, Hdx2y2)

    # d3z2r2
    write(fp,string(site_Cu[2])*" "*string(site_Cu[2])*" ")
    write_H(fp, Hd3z2r2)
end

function d_element_internal(fp, site_Cu, H_internal)
    H_internal_dx2y2d3z2r2 = H_internal[1]

    write_inv(fp, site_Cu[1], site_Cu[2], H_internal_dx2y2d3z2r2)
end

function p_element_internal(fp, site_O, H_internal)
    H_internal_pxpy = H_internal[1]
    H_internal_pypz = H_internal[2]
    H_internal_pzpx = H_internal[3]

    write_inv(fp, site_O[1], site_O[2], H_internal_pxpy)
    write_inv(fp, site_O[2], site_O[3], H_internal_pypz)
    write_inv(fp, site_O[3], site_O[1], H_internal_pzpx)
end

function main()
    filename = "corner_sharing_dp_deg.txt"
    # 単位胞あたりのサイト数
    ns = 11
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
    Upp = 0.0#4.0
    
    #J = 0.0#1.2

    tx_x = 0.5*ppπ + 0.5*ppσ
    tx_y = -0.5*ppπ + 0.5*ppσ
    tz_z = ppπ

    tx_x2y2 = sqrt(3.0)/2.0*pdσ
    tx_3z2r2 = 0.0#-0.5*pdσ
    tz_3z2r2 = 0.0#-0.5*pdσ
    ty_xy = pdπ

    Hx_x = Any[0;tx_x;0.0]
    Hx_y = Any[0;tx_y;0.0]
    Hz_z = Any[0;tz_z;0.0]

    Hx_x2y2 =  Any[0;tx_x2y2 ;0.0]
    Hx_3z2r2 = Any[0;tx_3z2r2;0.0]
    Hy_xy = Any[0;ty_xy;0.0]
    
    # 軌道の判定に必要
    Hpx = Any[1;Upp;0]
    Hpy = Any[2;Upp;0]
    Hpz = Any[3;Upp;0]

    #　同一軌道のクーロン反発と1体エネルギー
    Hdx2y2  = Any[4;Udd;Δ]
    Hd3z2r2 = Any[5;0.0;Δ]

    # 同一サイト内のクーロン反発
    H_internal_pxpy = Any[12;0.0;Upp]
    H_internal_pypz = Any[23;0.0;Upp]
    H_internal_pzpx = Any[31;0.0;Upp]
    H_internal_dx2y2d3z2r2 = Any[31;0.0;Udd]

    Hpp = [Hx_x, Hx_y, Hz_z]
    Hdp = [Hx_x2y2, Hx_3z2r2, Hy_xy]
    Hp = [Hpx, Hpy, Hpz]
    Hd = [Hdx2y2, Hd3z2r2]
    #Hdd = Any["dd";J]
    H_internal_pp = [H_internal_pxpy, H_internal_pypz, H_internal_pzpx]
    H_internal_dd = [H_internal_dx2y2d3z2r2]

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
        write(fp, string((i - 1)*ns + 1)*" "*string(O1[1])*" "*string(O1[2])*" 0\n")
        write(fp, string((i - 1)*ns + 2)*" "*string(O1[1])*" "*string(O1[2])*" 0\n")
        write(fp, string((i - 1)*ns + 3)*" "*string(O1[1])*" "*string(O1[2])*" 0\n")
        write(fp, string((i - 1)*ns + 4)*" "*string(O2[1])*" "*string(O2[2])*" 0\n")
        write(fp, string((i - 1)*ns + 5)*" "*string(O2[1])*" "*string(O2[2])*" 0\n")
        write(fp, string((i - 1)*ns + 6)*" "*string(O2[1])*" "*string(O2[2])*" 0\n")
        write(fp, string((i - 1)*ns + 7)*" "*string(Cu3[1])*" "*string(Cu3[2])*" 0\n")
        write(fp, string((i - 1)*ns + 8)*" "*string(Cu3[1])*" "*string(Cu3[2])*" 0\n")
        write(fp, string((i - 1)*ns + 9)*" "*string(O4[1])*" "*string(O4[2])*" 0\n")
        write(fp, string((i - 1)*ns + 10)*" "*string(O4[1])*" "*string(O4[2])*" 0\n")
        write(fp, string((i - 1)*ns + 11)*" "*string(O4[1])*" "*string(O4[2])*" 0\n")
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

        
        site1 = site_id[1:3]
        site2 = site_id[4:6]
        site3 = site_id[7:8]
        site4 = site_id[9:11]
        site5 = site_id[12:14]

        p_element(fp, site1, Hp)
        p_element(fp, site2, Hp)
        d_element(fp, site3, Hd)
        p_element(fp, site4, Hp)
        p_element(fp, site5, Hp)

        p_element_internal(fp, site1, H_internal_pp)
        p_element_internal(fp, site2, H_internal_pp)
        d_element_internal(fp, site3, H_internal_dd)
        p_element_internal(fp, site4, H_internal_pp)
        p_element_internal(fp, site5, H_internal_pp)

        pp_element(fp, -1.0, site1, site2, Hpp)
        pp_element(fp, +1.0, site1, site4, Hpp)
        pp_element(fp, +1.0, site2, site5, Hpp)
        pp_element(fp, -1.0, site4, site5, Hpp)

        dp_element_h(fp, +1.0, site3, site1, Hdp)
        dp_element_v(fp, -1.0, site3, site2, Hdp)
        dp_element_v(fp, +1.0, site3, site4, Hdp)
        dp_element_h(fp, -1.0, site3, site5, Hdp)
    end

    close(fp)
end

main()