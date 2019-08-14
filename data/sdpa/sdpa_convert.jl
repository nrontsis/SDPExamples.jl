# Helper Script that reads in the problems in the SDPA sparse file format and recreates and
# stores the problem data as matrices in Julia compatible JLD format

# The following is stored in the JLD file:
# m: dimension of x
# nblocks: number of blocks in the diagonal structure of the matrices
# blockvec: vector of numbers that give the sizes on the individual blocks
# c: objective vector
# F: array containing the constraint matrices, F[1]=F0, F[2]=F1, F[m+1]=Fm
# optval: Optimal objective value (extracted from README file)

# based on a earlier version from Michael Garstka
# University of Oxford, Control Group


using JLD2
using SparseArrays
# Specify path to .dat-s files
dirPath = "./"

# create array of all relevant filenames that contain problems
filenames = []
for f in filter(x -> endswith(x, "dat-s"), readdir(dirPath))
    f = split(f,".")[1]
    push!(filenames,f)
end

# filenames = ["example"]
for file in filenames
    # open file and read everything into memory
    f = open(dirPath*file*".dat-s");
    lines = readlines(f)
    close(f)

    # initialize variables
    counter = 1
    m = 0
    nblocks = Inf
    blocks = nothing
    c = nothing
    F = nothing # This is gonna hold the block matrices
    currentM = -1
    for ln in lines
        # dont do anything if line contains a comment
        if ln[1] == '"' || ln[1] == "*"
            counter == 1
        else
            # m: number of constraint matrices (in SDPA it starts from 0 -> +1)
            if counter == 1
                m = parse(Int64,strip(ln))
            # nblocks: number of blocks in the diagonal structure of the matrices
            elseif counter == 2
                nblocks = parse(Int64,strip(ln))

            # vector of numbers that give the sizes on the individual blocks
            # negative number indicates diagonal submatrix
            elseif counter == 3
                blocks_string = split(replace(strip(ln,['{','}','(',')']),"+" => ""))
                blocks = [parse(Int, s) for s in blocks_string]
                
                F = Array{SparseMatrixCSC{Float64, Int}}(undef, m + 1, nblocks)
                for i = 1:m + 1, j = 1:nblocks
                    F[i, j] = spzeros(abs(blocks[j]), abs(blocks[j]))
                end 

            # objective function vector c
            elseif counter == 4
                # try first to split by whitespace as delimiter (also remove certain trailing, ending characters)
                cvec = split(replace(strip(ln,['{','}','(',')']),"+" => ""))
                # otherwise try comma as delimiter
                if length(cvec) == 1
                    cvec = split(replace(strip(ln,['{','}','(',')']),"+" => ""),",")
                end
                c = [parse(Float64, ss) for ss in cvec]

            # all other lines contain constraint matrices with one entry per line
            # save them directly as sparse matrix
            else
                # FIXME: Accuracy
                line = [parse(Float64, ss) for ss in split(ln)]
                matrix_number = Int(line[1]) + 1
                block_number = Int(line[2])
                i = Int(line[3])
                j = Int(line[4])
                entry = line[5]

                F[matrix_number, block_number][i, j] = entry
                F[matrix_number, block_number][j, i] = entry
            end
            counter += 1
        end
    end

    # extract solution for objective value from README file
    f = open(dirPath*"README")
    lines = readlines(f)
    close(f)
    lines = lines[20:113]
    ln = filter(s -> !isa(match(Regex(file), s), Nothing), lines)
    if length(ln) == 0
        optVal = "Not provided"
        minimum_value = NaN
    else
        str = split(ln[1])[4]
        if !isa(match(Regex(str), "primal"), Nothing)
            minimum_value = Inf
        elseif !isa(match(Regex(str), "dual"), Nothing)
            minimum_value = -Inf
        else
            minimum_value = parse(Float64, str)
        end
    end

    # save to JLD2 file
    filepath = dirPath*file*".jld2"
    @save filepath m nblocks blocks c F minimum_value
    println("Saved problem: $(file) to "*dirPath*file*".jld2, Optimal Value: $(minimum_value)")
end