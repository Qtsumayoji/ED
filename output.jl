module output

    function write_data(X, Y, filename)
        dim_X = length(X)
        dim_Y = length(Y)

        if dim_X != dim_Y
            println("dimension error. dim_X=",dim_X," dim_Y=",dim_Y)
            return
        end

        dim = dim_X

        fp = open(filename, "w" )

        for i in 1:dim - 1
            write(fp, string(X[i])*" "*string(Y[i])*"\n")
        end
        write(fp, string(X[end])*" "*string(Y[end]))

        close(fp)
    end

    function write_data_serial(X, Y, Z, filename)
        dim_X = length(X)
        dim_Y = length(Y)
        dim_Z = size(Z)
        row = dim_Z[1]
        col = dim_Z[2]
        

        if dim_Z != (dim_X, dim_Y)
            println("dimension error. dim_X=",dim_X," dim_Y=",dim_Y," dim_Z=",dim_Z)
            return
        end

        fp = open(filename, "w" )

        write(fp, "\n"*"k,omega,reG,imG")
        for i in 1:row
            for j in 1:col
                x = string(X[i])
                y = string(Y[j])
                rez = string(real(Z[i, j]))
                imz = string(imag(Z[i, j]))
                write(fp,"\n"*x*","*y*","*rez*","*imz)
            end
        end

        close(fp)
    end

    function write_data_mat(X, Y, Z, filename)
        dim_X = length(X)
        dim_Y = length(Y)
        dim_Z = size(Z)
        row = dim_Z[1]
        col = dim_Z[2]
        

        if dim_Z != (dim_X, dim_Y)
            println("dimension error. dim_X=",dim_X," dim_Y=",dim_Y," dim_Z=",dim_Z)
            return
        end

        fp = open(filename, "w" )

        for i in 1:col
            y = string(Y[i])
            write(fp, " ,"*y)
        end

        for i in 1:row
            x = string(X[i])
            write(fp, "\n"*x)
            for j in 1:col
                rez = string(real(Z[i, j]))
                imz = string(imag(Z[i, j]))
                write(fp, rez*","*imz)
            end
        end

        close(fp)
    end
end