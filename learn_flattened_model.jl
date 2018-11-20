include("common.jl")
include("SmallFixedSizeSets.jl")
include("model_common.jl")

const MAX_BACK = 200

# Log-likelihood of the flattened model.
function log_likelihood(seqs::Vector{SequenceData},
                        repeat_seqs::Vector{SequenceData},
                        w::Vector{Float64})
    gradient = zeros(Float64, length(w), Threads.nthreads())
    ll = zeros(Float64, Threads.nthreads())
    total_choices = zeros(Int64, Threads.nthreads())    
    num_seqs = length(seqs)

    Threads.@threads for ind = 1:num_seqs
        seq = seqs[ind]
        repeat_seq = repeat_seqs[ind]
        if Threads.threadid() == 1
            @printf("%d of %d \r", ind, num_seqs)
            flush(stdout)
        end

        local_choice_count = 0
        seq_ll = 0.0
        seq_gradient = zeros(Float64, length(w))

        curr_sz_ind = 1
        for pos = 1:length(seq.sizes)
            
            sz = repeat_seq.sizes[pos]
            ran = curr_sz_ind:(curr_sz_ind + sz - 1)
            curr_sz_ind += sz
            if sz == 0; continue; end
            repeat_subset = SFSSFromOrderedVec(repeat_seq.sequence[ran])

            end_ind = 0
            for k in 1:(pos - 1); end_ind += seq.sizes[k]; end
            elmt_seq = seq.sequence[1:end_ind]

            # Only consider occurrences in MAX_BACK
            start = max(pos - MAX_BACK, 1)
            other_ind = 0
            for k in start:(pos - 1); other_ind += seq.sizes[k]; end
            elmt_seq = elmt_seq[(end_ind - other_ind + 1):end]

            local_weights = zeros(Float64, length(repeat_subset))
            occurrences = zeros(Float64, length(elmt_seq), length(repeat_subset))
            for (i, curr_elmt) in enumerate(reverse(elmt_seq))
                for (ind, elmt) in enumerate(repeat_subset)
                    if elmt == curr_elmt
                        local_weights[ind] += w[i]
                        occurrences[i, ind] = 1
                        break
                    end
                end
            end

            local_prob = 0.0
            grad_update = zeros(Float64, length(elmt_seq))
            for perm in permutations(collect(1:length(repeat_subset)))
                start_occurrences = ones(Float64, length(elmt_seq))
                part_gradient = zeros(Float64, length(elmt_seq), length(repeat_subset))
                w_norm = sum(w[1:end_ind])
                prob = 1.0
                probs = Float64[]
                for (part_ind, i) in enumerate(perm)
                    f = local_weights[i]
                    g = w_norm
                    new_prob = f / g
                    prob *= new_prob
                    push!(probs, new_prob)
                    w_norm -= local_weights[i]
                    part_gradient[:, part_ind] =
                        (g * occurrences[:, i] - f * start_occurrences) / g^2
                    start_occurrences -= occurrences[:, i]
                end
                rel_weights = prob ./ probs
                grad_update += part_gradient * rel_weights
                local_prob += prob
            end

            if local_prob > 0.0
                # Local prob could be 0 if more than MAX_BACK back
                local_choice_count += 1
                seq_ll += log(local_prob)
                seq_gradient[1:length(elmt_seq)] += grad_update / local_prob
            end
        end
        total_choices[Threads.threadid()] += local_choice_count
        ll[Threads.threadid()] += seq_ll
        gradient[:, Threads.threadid()] += seq_gradient
    end
    
    return sum(ll), vec(sum(gradient, dims=2)), sum(total_choices)
end

"""
Learn the flattened model for a dataset.

dataset: name of dataset

Optional parameters
max_iter: maximum number of gradient steps
init_step_size: initial step size of gradient descent
reduction: reduction in step size with unsuccessful steps
min_step_size: minimum step size
tol: stopping tolerance (relative likelihood change)
"""
function learn(dataset::AbstractString,
               max_iter::Int64=100,
               init_step_size::Float64=1.0,
               reduction::Float64=0.5,
               min_step_size::Float64=1e-15,
               tol::Float64=1e-3)
    # Preprocess data
    seqs = read_data(dataset)
    packed_seqs = pack_data(seqs)
    packed_repeat_seqs = packed_repeat_set_seqs(seqs)

    # get total number of choices
    max_seq_len = maximum([length(seq.sequence) for seq in packed_seqs])    
    recency_weights = ones(Float64, max_seq_len)
    recency_weights /= sum(recency_weights)
    curr_ll, grad, total_choices =
        log_likelihood(packed_seqs, packed_repeat_seqs, recency_weights)
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
            next_ll, next_grad = log_likelihood(packed_seqs, packed_repeat_seqs,
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
        if step_size < min_step_size; break; end
        if minimal_change; break; end
    end

    save("models/$(dataset)-flattened.jld2",
         Dict("w"              => recency_weights,
              "log_likelihood" => curr_ll,
              "iterations"     => iter,
              "total_choices"  => total_choices))
end
;
