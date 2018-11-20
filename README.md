# Sequences of Sets

This code and data repository accompanies the paper

- [Sequences of Sets](http://cs.cornell.edu/~arb/sequences-of-sets-KDD-2018.pdf). Austin R. Benson, Ravi Kumar, and Andrew Tomkins. Proceedings of the 24th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining (KDD '18'), 2018.

All of the code is written in Julia.

For questions, please email Austin at arb@cs.cornell.edu.

### Data

The datasets are in the `data/` directory. The file `data/dataset-seqs.txt` is a list of sequences. Each line of the file has the following form:

size1,size2,…,sizeN;elmt1,elmt2,…,elmtM

- size1,size2,…,sizeN are the number of elements in the N sets in the sequence.
- elmt1,elmt2,…,elmtM are the M elements (given as integer identifiers) in the N sets in order. The first size1 elements are in the first set, the next size2 elements are in the second set, and so on.
  For each sequence, size1 + size2 + … + sizeN = M.

The files `email-Enron-core-element-labels.txt`, `tags-math-sx-element-labels.txt`, and `tags-mathoverflow-element-labels.txt` contain labels for the element.

```
bash-3.2$ head -5 email-Enron-core-element-labels.txt 
1 phillip.allen@enron.com
2 john.arnold@enron.com
3 harry.arora@enron.com
4 robert.badeer@enron.com
5 susan.bailey@enron.com
```



### Learning models

##### Correlated Repeated Unions (CRU) model

We first show how to learn the Correlated Repeated Unions (CRU) model. The learning has some built-in parallelism, which you can use by setting the JULIA_NUM_THREADS environment variable.

```julia
include("learn_CRU_model.jl")
dataset = "email-Enron-core"  # data file at data/$dataset.txt
p = 0.9                       # correlation probability
# Learning takes several minutes
learn(dataset, p)  # --> model saved to models/$dataset-CRU-$p.jld2
```

All of the learned CRU models used in the paper are pre-computed and saved in the `output/` directory.

##### Flattened model

We use a "flattened" model as the baseline.

```julia
include("learn_flattened_model.jl")
dataset = "email-Enron-core"  # data file at data/$dataset.txt
# Learning takes several minutes
learn(dataset)  # --> model saved to models/$dataset-flattened.jld2
```



### Reproduce the figures and tables in the paper

##### Figure 1: Distribution of set sizes.

```julia
include("paper_figures.jl")
set_size_dist_fig()  # --> set_size_dist.pdf
```

##### Figure 2: Repeat behavior in the datasets

```julia
include("paper_figures.jl")
# The following takes a minute or so
repeat_behavior_fig()  # --> repeat_behavior.pdf
```

##### Figure 3: Distribution of the number of repeatss in sets containing at least one repeat.

```julia
include("paper_figures.jl")
num_repeats_dist_fig()  # --> num_repeats_dist.pdf 
```

##### Figure 4: Evidence of recency bias in set selection.

```julia
include("paper_figures.jl")
# The following takes a minute or so
recency_bias_fig()  # --> recency_bias.pdf
```

##### Figure 5: Likelihoods.

```julia
include("paper_figures.jl")
for row in dataset_info()
    dataset = row[1]
    likelihoods_fig(dataset)  # --> $dataset-rel-likelihoods.pdf
end
```

##### Figure 6: Recency weights.

```julia
include("paper_figures.jl")
for row in dataset_info()
    dataset = row[1]
    recency_weights_fig(dataset)  # --> weights-$dataset.pdf
end
```

##### Table 1: Summary statistics of datasets. 

```julia
include("paper_tables.jl")
for row in dataset_info()
    summary_stats(row[1])
end
```

##### Table 2: Subset correlations.

```julia
include("paper_tables.jl")
# This takes several minutes
for row in dataset_info()
    correlation_behavior(row[1])
end 
```