include("common.jl")
include("SmallFixedSizeSets.jl")

using Combinatorics
using DataStructures: counter
using Printf

# Summary statistics of a dataset.
function summary_stats(dataset::String)
    all_sets = Set{SmallFixedSizeSet}()
    seqs = read_data(dataset)
    for seq in seqs
        for subset in seq
            v = sort(collect(subset))
            push!(all_sets, SFSSFromOrderedVec(v))
        end
    end
    packed_seqs = pack_data(seqs)
    all_items = Int64[]
    num_sets = 0
    for seq in packed_seqs
        append!(all_items, seq.sequence)
        num_sets += length(seq.sizes)
    end
    num_seqs = length(packed_seqs)
    univ_size = length(unique(all_items))
    num_uniq_sets = length(all_sets)
    println("dataset: $dataset")
    println("  number of sequences: $num_seqs")
    println("  universe size: $univ_size")
    println("  number of sets: $num_sets")
    println("  number of unique sets: $num_uniq_sets")
end

# Subset correlations in a dataset.
function correlation_behavior(dataset::String, num_rand_samples::Int64=100)
    function get_counts(sizes::Vector{Int64}, items::Vector{Int64},
                        tupsize::Int64, shuffle::Bool=false)
        counts = counter(NTuple{tupsize, Int64})
        new_items = copy(items)
        if shuffle; shuffle!(new_items); end
        curr_sz_ind = 1
        for sz in sizes
            vec = new_items[curr_sz_ind:(curr_sz_ind + sz - 1)]
            curr_sz_ind += sz
            for tup in combinations(vec, tupsize)
                push!(counts, NTuple{tupsize, Int64}(sort(tup)))
            end
        end
        return [v for (k, v) in counts]
    end

    println("$dataset...")
    packed_seqs = read_packed_data(dataset)
    
    w2, w3 = Int64[], Int64[]
    w2_null = Vector{Vector{Int64}}(num_rand_samples)
    w3_null = Vector{Vector{Int64}}(num_rand_samples)
    for j in 1:num_rand_samples
        w2_null[j], w3_null[j] = Vector{Int64}(), Vector{Int64}()
    end        
    
    for (ind, seq) in enumerate(packed_seqs)
        print(@sprintf("%d of %d \r", ind, length(packed_seqs)))
        flush(STDOUT)
        sz2, sz3 = Int64[], Int64[]
        itm2, itm3 = Int64[], Int64[]
        curr_ind = 1
        sequence = seq.sequence
        for sz in seq.sizes
            subset_vec = sequence[curr_ind:(curr_ind + sz - 1)]
            curr_ind += sz
            if sz >= 2
                push!(sz2, sz)
                append!(itm2, subset_vec)
            end
            if sz >= 3
                push!(sz3, sz)
                append!(itm3, subset_vec)
            end
        end
        # For every pair, count how many times it appears
        for v in get_counts(sz2, itm2, 2, false); push!(w2, v); end
        for v in get_counts(sz3, itm3, 3, false); push!(w3, v); end
        for j in 1:num_rand_samples
            for v in get_counts(sz2, itm2, 2, true); push!(w2_null[j], v); end
            for v in get_counts(sz3, itm3, 3, true); push!(w3_null[j], v); end
        end
    end

    println(@sprintf("mean (2): %f", mean(w2)))
    println(@sprintf("mean (3): %f", mean(w3)))
    m2 = [mean(w2_null[j]) for j in 1:num_rand_samples]
    m3 = [mean(w3_null[j]) for j in 1:num_rand_samples]
    println(@sprintf("mean, null (2): %f +/- %f", mean(m2), std(m2)))
    println(@sprintf("mean, null (3): %f +/- %f", mean(m3), std(m3)))
    println("------")
end
;
