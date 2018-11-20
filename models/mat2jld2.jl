# Convert mat to jld2 with v0.6 julia
using MAT
using FileIO, JLD2
function main()
    for file in readdir(".")
        if file[(end-3):end] == ".mat"
            @show file
            save(string(join(split(file, ".")[1:(end-1)], '.'), ".jld2"), matread(file))
        end
    end
end
