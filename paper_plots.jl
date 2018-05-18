include("common.jl")

using PyPlot
using StatsBase: countmap

# Distribution of set sizes
function set_size_dists()
    close()
    for row in dataset_info()
        dataset = row[1]
        packed_seqs = read_packed_data(dataset)
        all_sizes = Int64[]
        for seq in packed_seqs
            curr_sizes = seq.sizes
            append!(all_sizes, curr_sizes[curr_sizes .> 0])
        end
        nums, counts = zip(sort(countmap(all_sizes))...)
        tot = sum(counts)
        fracs = [count / tot for count in counts]
        @show nums, fracs
        semilogy(collect(nums), fracs, lw=2, marker=row[2],
                 ms=8, label=dataset)
    end
    fsz = 20
    ylabel("Fraction", fontsize=fsz)
    xlabel("Size of set", fontsize=fsz)
    legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.,
           frameon=false, fontsize=fsz-2)
    ax = gca()
    ax[:set_xticks](1:5)
    ax[:tick_params](axis="both", labelsize=fsz-2, length=7, width=1.5)
    ax[:tick_params](axis="y",    which="minor", length=4, width=1)    
    
    savefig("set_size_dists.pdf", bbox_inches="tight")
    show()
end

# Repeat behavior
function repeat_behavior()
    close()
    fsz=12

    function adjust_axis(ax, ylabel::String)
        ax[:set_xlabel]("Set size", fontsize=fsz)
        ax[:set_ylabel](ylabel, fontsize=fsz)
        ax[:set_xticks](1:5)
        ax[:tick_params](axis="both", labelsize=fsz, length=7, width=1.5)
        ax[:set_ylim](0, 1.0)
    end

    info = dataset_info()
    all_seqs = [read_data(row[1]) for row in info]

    # Exact repeats
    subplot(221)
    for (ind, row) in enumerate(info)
        seqs = all_seqs[ind]
        total, repeated = zeros(Int64, 5), zeros(Int64, 5)
        for seq in seqs
            for i in 1:length(seq)
                curr_set = seq[i]
                curr_size = length(curr_set)
                for j in (i - 1):-1:1
                    if seq[j] == curr_set
                        repeated[curr_size] += 1
                        break
                    end
                end
                total[curr_size] += 1
            end
        end
        plot(collect(1:5), repeated ./ total, marker=row[2], label=row[1])
    end
    adjust_axis(gca(), "Fraction exact repeats")

    # Entirely novel
    subplot(222)
    for (ind, row) in enumerate(info)
        seqs = all_seqs[ind]
        total, disjoint = zeros(Int64, 5), zeros(Int64, 5)
        for seq in seqs
            for i in 1:length(seq)
                curr_set = seq[i]
                curr_size = length(curr_set)
                is_disjoint = true
                for j in (i - 1):-1:1
                    if length(curr_set ∩ seq[j]) > 0
                        is_disjoint = false
                        break
                    end
                end
                if is_disjoint; disjoint[curr_size] += 1; end
                total[curr_size] += 1
            end
        end
        plot(collect(1:5), disjoint ./ total, marker=row[2], label=row[1])
    end
    adjust_axis(gca(), "Fraction entirely novel")

    # Subsets
    subplot(223)
    for (ind, row) in enumerate(info)
        seqs = all_seqs[ind]
        total, subsets = zeros(Int64, 5), zeros(Int64, 5)
        for seq in seqs
            for i in 1:length(seq)
                curr_set = seq[i]
                curr_size = length(curr_set)
                for j in (i - 1):-1:1
                    if issubset(curr_set, seq[j])
                        subsets[curr_size] += 1
                        break
                    end
                end
                total[curr_size] += 1
            end

        end
        plot(collect(1:5), subsets ./ total, marker=row[2], label=row[1])
    end
    adjust_axis(gca(), "Fraction subset")    

    # Supersets    
    subplot(224)
    for (ind, row) in enumerate(info)
        seqs = all_seqs[ind]
        total, supersets = zeros(Int64, 5), zeros(Int64, 5)
        for seq in seqs
            for i in 1:length(seq)
                curr_set = seq[i]
                curr_size = length(curr_set)
                for j in (i - 1):-1:1
                    if issubset(seq[j], curr_set)
                        supersets[curr_size] += 1
                        break
                    end
                end
                total[curr_size] += 1
            end
        end
        plot(collect(1:5), supersets ./ total, marker=row[2], label=row[1])
    end
    adjust_axis(gca(), "Fraction superset")

    tight_layout()
    savefig("repeat_behavior.pdf", bbox_inches="tight")
    show()
end
