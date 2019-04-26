module Make_Anderson_model
    include("./Model.jl")

    function main(μ, U, Ef, Vl)
        filename = "Anderson_model.txt"
        # 単位胞あたりのサイト数
        ns = 6
        Nx = 1
        Ny = 1
        Ns = ns

        ax = 1.0
        ay = 1.0

        H11 = Any[1;μ;U]
        H1i = Any[2;Vl]
        Hii = Any[3;Ef]

        fp = open(filename, "w")

        write(fp, "ns="*string(ns)*"\n")
        write(fp, "Nx="*string(Nx)*"\n")
        write(fp, "Ny="*string(Ny)*"\n")
        write(fp, "Nz="*string(1)*"\n")
        write(fp, "Ns="*string(Ns)*"\n")

        write(fp, string(1)*" "*string(0)*" "*string(0)*" 0\n")

        for i in 2:Ns
            write(fp, string(i)*" "*string(cos(2(i-2)/(Ns-1)*pi))*" "*string(sin(2(i-2)/(Ns-1)*pi))*" 0\n")
        end

        write(fp, string(1)*" "*string(1)*" "*string(H11[1])*" "*string(H11[2])*" "*string(H11[3])*"\n")

        for i in 2:Ns
            write(fp, string(i)*" "*string(i)*" "*string(Hii[1])*" "*string(Hii[2])*" 0.0\n")
            write(fp, string(1)*" "*string(i)*" "*string(H1i[1])*" "*string(H1i[2])*" 0.0\n")
            write(fp, string(i)*" "*string(1)*" "*string(H1i[1])*" "*string(H1i[2])*" 0.0\n")
        end
        
        close(fp)
    end
end
    #main(1.0,2.0,3.0,4.0)