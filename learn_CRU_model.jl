include("common.jl")
include("model_common.jl")

using Base.Threads
using Combinatorics
using MAT
using SmallFixedSizeSets

# maximum look back for repeats
const MAX_BACK = 200

# All unordered partitions of a list of integers.
function unordered_partitions(A::Vector{Int64})
    s = length(A)
    partitions = Set{Set{Set{Int32}}}()
    # base case
    if s == 0;
        push!(partitions, Set{Set{Int32}}())
        return partitions
    end
    
    for int_part in integer_partitions(s)
        for part_size in int_part
            for part_item in combinations(collect(1:s), part_size)
                curr_set = Set{Int32}(A[part_item])
                remainder = ones(Bool, s)
                remainder[part_item] = false
                # Recursion
                for remaining_part in unordered_partitions(A[remainder])
                    push!(remaining_part, curr_set)
                    push!(partitions, remaining_part)
                end
            end
        end
    end
    
    return partitions
end

# All ordered partitions of a given size.
function ordered_partitions(size::Int64)
    
    partitions = Vector{Vector{Vector{Int32}}}()
    for set_collection in unordered_partitions(collect(1:size))
        curr_part = Vector{Vector{Int32}}()
        for set in set_collection
            push!(curr_part, sort(collect(set)))
        end
        # account for orderings
        append!(partitions, collect(permutations(curr_part)))
    end
    return partitions
end

"""
Compute probability of repeating a sample

A: accumulation set
p: correlation probability
w: recency weights
pos: position in sequence to evaluate likelihood
"""
function repeat_prob(seqdata::SequenceData, A::SmallFixedSizeSet, pos::Int64,
                     p::Float64)
    sequence = seqdata.sequence
    sizes = seqdata.sizes
    curr_sz_ind = 1
    positional_probs = zeros(Float64, pos - 1)
    start = max(pos - MAX_BACK, 1)
    for k in 1:(start - 1); curr_sz_ind += sizes[k]; end
    for ind = start:(pos - 1)
        sz = sizes[ind]
        subset_vec = sequence[curr_sz_ind:(curr_sz_ind + sz - 1)]
        curr_sz_ind += sz
        subset = SFSSFromOrderedVec(subset_vec)

        s = length(subset_vec)
        positional_prob = (1 - p)^s  # initialize with size-0 case

        # Performance optimization: skip if A does not intersect with the subset.
        if length(SFSS_intersect(A, subset)) > 0
            end_ind = min(length(A), s)
            for r = 1:end_ind
                count = 0
                for R in combinations(subset_vec, r)
                    count += SFSS_issubset(SFSSFromVec(R), A)
                end
                positional_prob += count * p^r * (1 - p)^(s - r)
            end
        end
        positional_probs[ind] = positional_prob
    end
    return positional_probs
end

"""
Compute the probability of success.

A: accumulation set
B: next set to add
pos: position in sequence
w: recency weights
p: correlation probability
is_last: whether or not B is the last one in the partition
"""
function success_prob(seqdata::SequenceData, A::SmallFixedSizeSet,
                      B::SmallFixedSizeSet, pos::Int64, p::Float64,
                      is_last::Bool)
    AB = SFSS_union(A, B)
    bsize = length(B)
    sequence = seqdata.sequence
    sizes = seqdata.sizes
    curr_sz_ind = 1
    positional_probs = zeros(Float64, pos - 1)
    start = max(pos - MAX_BACK, 1)
    for k in 1:(start - 1); curr_sz_ind += sizes[k]; end
    for ind = start:(pos - 1)
        sz = sizes[ind]
        subset_vec = sequence[curr_sz_ind:(curr_sz_ind + sz - 1)]
        curr_sz_ind += sz
        subset = SFSSFromOrderedVec(subset_vec)
        # Performance optimization: skip if B is not contained in the subset
        positional_prob = 0.0        
        if SFSS_issubset(B, subset)
            s = length(subset_vec)
            # B is a subset, so we have to consider it start at
            # size of B (otherwise, B won't be subset)
            start = max(1, length(B))
            for r = start:s
                count = 0
                for R in combinations(subset_vec, r)
                    Rs = SFSSFromVec(R)
                    if SFSS_issubset(B, Rs)
                        if is_last
                            C = SFSS_setdiff(Rs, AB)
                            csize = length(C)
                            bcsize = bsize + csize
                            count += factorial(bsize) * factorial(csize) / factorial(bcsize)
                        else
                            count += SFSS_issubset(Rs, AB)
                        end
                    end
                end
                positional_prob += count * p^r * (1 - p)^(s - r)
            end
        end
        positional_probs[ind] = positional_prob
    end
    return positional_probs
end

# Compute the log likelihood for the model parameters.
function log_likelihood(seqs::Vector{SequenceData},
                        repeat_seqs::Vector{SequenceData},
                        p::Float64,
                        w::Vector{Float64})
    # Get all sorted partitions of a set
    partitions = [ordered_partitions(i) for i in 1:5]
    function sorted_set_partitions(S::SmallFixedSizeSet)
        set_parts = Vector{Vector{SmallFixedSizeSet}}()
        for partition in partitions[length(S)]
            push!(set_parts, [S[inds] for inds in partition])
        end
        return set_parts
    end

    gradient = zeros(Float64, length(w), Threads.nthreads())
    ll = zeros(Float64, Threads.nthreads())
    total_choices = zeros(Int64, Threads.nthreads())    
    num_seqs = length(seqs)
    Threads.@threads for ind = 1:num_seqs
        seq = seqs[ind]
        repeat_seq = repeat_seqs[ind]
        if Threads.threadid() == 1
            print(@sprintf("%d of %d \r", ind, num_seqs))
        end
        flush(STDOUT)

        local_choice_count = 0
        seq_ll = 0.0
        curr_sz_ind = 1
        seq_gradient = zeros(Float64, length(w))
        for pos = 1:length(seq.sizes)
            sz = repeat_seq.sizes[pos]
            ran = curr_sz_ind:(curr_sz_ind + sz - 1)
            curr_sz_ind += sz
            if sz == 0; continue; end
            repeat_subset = SFSSFromOrderedVec(repeat_seq.sequence[ran])
            scaling_sizes = seq.sizes[1:(pos - 1)]
            w_curr = w[1:(pos - 1)] .* reverse(scaling_sizes)
            weight_normalization = sum(w_curr)
            local_prob = 0.0
            local_gradient = zeros(Float64, pos - 1)

            for partition in sorted_set_partitions(repeat_subset)
                A = Set0()  # accumulation set
                prob = 1.0
                num_parts = length(partition)
                part_gradient = zeros(Float64, pos - 1, num_parts)
                weight_probs = zeros(Float64, num_parts)
                for j = 1:num_parts
                    B = partition[j]
                    
                    p_s = success_prob(seq, A, B, pos, p, j == num_parts)
                    p_r = repeat_prob(seq, A, pos, p)
                    q_s = sum(w_curr .* reverse(p_s))
                    q_r = sum(w_curr .* reverse(p_r))
                    
                    num = q_s
                    denom = weight_normalization - q_r
                    iterate_prob = num / denom
                    
                    prob *= iterate_prob
                    A = SFSS_union(A, B)

                    weight_probs[j] = iterate_prob
                    part_gradient[:, j] =
                        (p_s * denom - num * (scaling_sizes - p_r)) / denom^2
                end
                if prob > 0.0
                    local_prob += prob
                    rel_weights = prob ./ weight_probs
                    local_gradient += part_gradient * rel_weights
                else
                    zero_inds = find(weight_probs .== 0.0)
                    if length(zero_inds) == 1
                        zero_ind = zero_inds[1]
                        other_prob = 1.0
                        for (prob_ind, wp) in enumerate(weight_probs)
                            if prob_ind != zero_ind; other_prob *= wp; end
                        end
                        local_gradient += part_gradient[:, zero_ind] * other_prob
                    end
                end
                    
            end
            if local_prob > 0.0
                # Local prob could be 0 if more than MAX_BACK back
                local_choice_count += 1
                seq_ll += log(local_prob)
                seq_gradient[1:(pos - 1)] += reverse(local_gradient) / local_prob
            end
        end
        total_choices[Threads.threadid()] += local_choice_count
        ll[Threads.threadid()] += seq_ll
        gradient[:, Threads.threadid()] += seq_gradient
    end
    
    return sum(ll), vec(sum(gradient, 2)), sum(total_choices)
end

"""
Learn the CRU model for a dataset.

dataset: name of dataset
p: correlation probability

Optional parameters
max_iter: maximum number of gradient steps
init_step_size: initial step size of gradient descent
reduction: reduction in step size with unsuccessful steps
min_step_size: minimum step size
tol: stopping tolerance (relative likelihood change)
"""
function learn(dataset::String, p::Float64,
               max_iter::Int64=100,
               init_step_size::Float64=1.0,
               reduction::Float64=0.5,
               min_step_size::Float64=1e-5,
               tol::Float64=1e-3)
    # Preprocess data
    seqs = read_data(dataset)
    packed_seqs = pack_data(seqs)
    packed_repeat_seqs = packed_repeat_set_seqs(seqs)

    # get total number of choices
    max_seq_len = maximum([length(seq) for seq in seqs])
    recency_weights = ones(Float64, max_seq_len - 1)
    recency_weights /= sum(recency_weights)
    curr_ll, grad, total_choices =
        log_likelihood(packed_seqs, packed_repeat_seqs, p, recency_weights)
    println("total choices: $(total_choices)")
    
    grad /= (norm(grad, 1) + 1)
    step_size = init_step_size
    println("log likelihood: $(curr_ll)")
    iter = 1
    minimal_change = false
    for iter = 1:max_iter
        println("iteration $iter")
        while true
            test_recency_weights =
                simplex_projection(recency_weights + step_size * grad)
            next_ll, next_grad =
                log_likelihood(packed_seqs, packed_repeat_seqs, p,
                               test_recency_weights)[1:2]
            next_grad /= (norm(next_grad, 1) + 1)
            if next_ll > curr_ll
                change = exp((next_ll - curr_ll) / total_choices)
                minimal_change = change < 1 + tol
                println("ll: $(next_ll), change: $(change)")
                curr_ll = next_ll
                grad = next_grad
                recency_weights = test_recency_weights
                break
            else
                step_size *= reduction
            end
            if step_size < min_step_size; break; end
        end
        if minimal_change; step_size *= reduction; end
        if step_size < min_step_size; break; end
    end

    matwrite("output/$(dataset)-$(p).mat",
             Dict("w"              => recency_weights,
                  "p"              => p,
                  "log_likelihood" => curr_ll,
                  "iterations"     => iter,
                  "total_choices"  => total_choices))
end
;
