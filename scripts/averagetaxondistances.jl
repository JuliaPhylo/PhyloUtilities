"""
Utilities for calculating pairwise taxon distances from multiple gene trees,
normalizing them to account for rate variation across genes
(slow vs fast genes)
using the median distance between ingroup and outgroup taxa,
then averaging these normalized taxon distances across genes.
Some taxa could be missing for some genes, but each gene should have at
least 1 ingroup and 1 outgroup taxon.

by: Cécile Ané, 2025, using PhyloNetworks v1.1.0.

Functions defined here:

    getpairwisedistances
    normalizedistances_outgroup2ingroup!
    averagepairwisedistances

These functions were originally developed for
[Karimi et al. (2020)](https://doi.org/10.1093/sysbio/syz073).
Their [code](https://github.com/nkarimi/Adansonia_HybSeq/blob/master/trait-evolution/calibration.jl)
was slightly modified here.

See their documentation below for details and examples.
"""

"""
    getpairwisedistances(genetrees, taxa)

Return a tuple of 3 objects:
1. A vector `D` of matrices, one per gene tree, containing the pairwise distance
   between all pairs of `taxa`. If taxon i is missing from tree `k`, then the
   distance matrix `D[k]` for that tree will have zeros on its ith row and ith column.
   In each matrix, row & column `i` correspond to `taxa[i]`, that is, taxa are
   listed in the same order as in the input `taxa`.
2. A matrix `ngenes` containing the number of genes with data on the pair
   of taxa (i,j) in row i, column j
3. A vector of integers, giving the index of gene trees with some missing taxa.

This function uses `pairwisetaxondistancematrix(tree)` from `PhyloNetworks`,
which outputs a matrix in which the rows correspond to taxa in the order in
which they come in `tipLabels(tree)`.
It then takes care of the fact that taxa may not be listed in the same order by
`tipLabels` across all gene trees.

# example

```@repl
julia> using DataFrames, PhyloNetworks

julia> include("scripts/averagetaxondistances.jl")

julia> genetrees = readnewick.([ # read 3 gene trees
"(P_micranthum:0.0180,P_pulcherimum:0.0077,(P_carneum:0.0022,P_elegans:0.0009):0.0084,((P_eddyense:0.0022,P_chartaceum:0.0022):0.0034,((P_reptans:0.0009,P_occidentale:0.0074):0.0038,P_pectinatum:0.0077):0.0021,(P_pauciflorum:0.0009,P_filicinum:0.0008):0.0048):0.0047,(P_delicatum:0.0034,(P_viscosum:0.0035,P_elusum:0.0009):0.0025):0.0043);",
"(P_micranthum:0.0175,P_pauciflorum:0.0046,P_molle:0.0030,P_viscosum:0.0009,P_filicinum:0.0032,P_eddyense:0.0008,P_viscosum:0.0008,P_brandegeii:0.0008,P_chartaceum:0.0009,P_occidentale:0.0024,P_pectinatum:0.0051,P_pulcherimum:0.0034,P_carneum:0.0033,P_elegans:0.0009,P_delicatum:0.0009,P_elusum:0.0022,P_reptans:0.0022);",
"(P_micranthum:0.0095,((P_viscosum:0.0009,P_delicatum:0.0067,(P_pulcherimum:0.0031,P_elegans:0.0008):0.0042,(P_eddyense:0.0046,P_elusum:0.0019):0.0017,(P_brandegeii:0.0031,P_molle:0.0008,P_viscosum:0.0008):0.0053,P_reptans:0.0020):0.0020,(P_carneum:0.0031,P_pauciflorum:0.0020,P_chartaceum:0.0008,P_filicinum:0.0015):0.0029):0.0021,(P_pectinatum:0.0020,P_occidentale:0.0022):0.0020);",
]);

julia> taxa = sort!(union(tiplabels.(genetrees)...));

julia> D, ngenes, geneind = getpairwisedistances(genetrees, taxa);

julia> geneind # tree 1 is missing some taxa
1-element Vector{Int64}:
 1

```
"""
function getpairwisedistances(genetrees, taxa)
  ntips = length(taxa)
  D = Array{Float64,2}[]; # empty vector. will contain all distance matrices
  ngenes = zeros(Int, ntips, ntips) # number of genes that have data for each pair
  geneind = Int[];        # indices of genes with missing taxa
  istaxonmissing = Vector{Bool}(undef, ntips) # to be modified in place for each gene
  for (g_index,g) in enumerate(genetrees)
    M = zeros(ntips,ntips) # initialized at 0.0: for missing pairs
    taxnames = tipLabels(g)
    tipind = Int[]
    for k in 1:ntips
      j = findfirst(isequal(taxa[k]), taxnames)
      notfound = isnothing(j)
      istaxonmissing[k] = notfound # modified in place
      notfound || push!(tipind, j) # add j to tipind if taxa[k] was found
    end
    M[.!istaxonmissing, .!istaxonmissing] = pairwisetaxondistancematrix(g)[tipind,tipind]
    ngenes[.!istaxonmissing, .!istaxonmissing] .+= 1
    any(istaxonmissing) && push!(geneind,g_index)
    push!(D, M)
  end
  return D, ngenes, geneind
end

"""
    normalizedistances_outgroup2ingroup!(D; taxa, ingroup, outgroup)

Rescale each input distance matrix `D[k]`, such that all have the same
median patristic distance between outgroup taxa and ingroup taxa.
Input: `D` should be a vector of pairwise distances matrices, one per gene
(modified by the function).
Output: vector of original median ingroup-outgroup distance, one per gene.

Why the *median*? So that one taxon or one small clade with an unusually large
(or low) substitution rate does have an undue influence on the scaling factor.

Assumptions:
- all trees have at least 1 outgroup and 1 ingroup
- row & column `i` in D[k] (for gene k) correspond to `taxa[i]`
- `D[k][i,j]` = 0 if gene `k` doesn't have both taxa `i` and `j`
- `ingroup` and `outgroup` are sets. The function does *not* check whether they
  are subsets of `taxa`, or don't overlap, or cover the full set of `taxa`.
"""
function normalizedistances_outgroup2ingroup!(D; taxa, ingroup, outgroup)
  ntax = length(taxa)
  inding = findall(in(ingroup),  taxa) # indices of ingroup  taxa
  indout = findall(in(outgroup), taxa) # indices of outgroup taxa
  medianingroup2outgroup = Float64[]   # will contain 1 median per gene
  for dm in D # dm = distance matrix
    size(dm) = (ntax,ntax) || error("there's a distance matrix with wrong dimensions: $(size(dm))")
    absent = findall([all(dm[:,i] .== 0.0) for i in 1:ntax])
    push!(medianingroup2outgroup,
          median(dm[setdiff(inding, absent), setdiff(indout, absent)]) )
  end
  mi2o = mean(medianingroup2outgroup)
  for k in eachindex(D)
    D[k] .*= mi2o/medianingroup2outgroup[k]
  end
  return medianingroup2outgroup
end

"""
    averagepairwisedistances(D, ngenes)

Matrix `M` containing the average pairwise distances, weighted by number of
genes with data on the pair: M[i,j] = (sum_k D[k][i,j] ) / ngenes[i,j].
This is because for each pair of taxa `i,j`, it is assumed that a number
`ngenes[i,j]` of genes (indexed by k) contributed data for the pair, and
the other genes without both taxa i and j had D[k][i,j]=0.

# example

```@repl
julia> using CSV, DataFrames, PhyloNetworks, StatsBase

julia> include("scripts/averagetaxondistances.jl")

julia> genetrees = readnewick.([ # read 3 gene trees
"(P_micranthum:0.0180,P_pulcherimum:0.0077,(P_carneum:0.0022,P_elegans:0.0009):0.0084,((P_eddyense:0.0022,P_chartaceum:0.0022):0.0034,((P_reptans:0.0009,P_occidentale:0.0074):0.0038,P_pectinatum:0.0077):0.0021,(P_pauciflorum:0.0009,P_filicinum:0.0008):0.0048):0.0047,(P_delicatum:0.0034,(P_viscosum:0.0035,P_elusum:0.0009):0.0025):0.0043);",
"(P_micranthum:0.0175,P_pauciflorum:0.0046,P_molle:0.0030,P_viscosum:0.0009,P_filicinum:0.0032,P_eddyense:0.0008,P_viscosum:0.0008,P_brandegeii:0.0008,P_chartaceum:0.0009,P_occidentale:0.0024,P_pectinatum:0.0051,P_pulcherimum:0.0034,P_carneum:0.0033,P_elegans:0.0009,P_delicatum:0.0009,P_elusum:0.0022,P_reptans:0.0022);",
"(P_micranthum:0.0095,((P_viscosum:0.0009,P_delicatum:0.0067,(P_pulcherimum:0.0031,P_elegans:0.0008):0.0042,(P_eddyense:0.0046,P_elusum:0.0019):0.0017,(P_brandegeii:0.0031,P_molle:0.0008,P_viscosum:0.0008):0.0053,P_reptans:0.0020):0.0020,(P_carneum:0.0031,P_pauciflorum:0.0020,P_chartaceum:0.0008,P_filicinum:0.0015):0.0029):0.0021,(P_pectinatum:0.0020,P_occidentale:0.0022):0.0020);",
]);

julia> taxa = sort!(union(tiplabels.(genetrees)...));

julia> D, ngenes, geneind = getpairwisedistances(genetrees, taxa);

julia> outgroup = ["P_micranthum"] # there could be more than 1

julia> ingroup = setdiff(taxa, outgroup)

julia> med_in2out = normalizedistances_outgroup2ingroup!(D,
            taxa=taxa, ingroup=ingroup, outgroup=outgroup);

julia> avD = averagepairwisedistances(D, ngenes);

julia> avD[1:2,1:4] # distance betwen taxa 1:2 and taxa 1:4
2×4 Matrix{Float64}:
 0.0        0.0124633  0.00970258  0.0103238
 0.0124633  0.0        0.00855209  0.0123552

julia> avD_df = DataFrame( # convert to a data frame with column names
        [avD[:,j] for j in 1:size(avD,2)], # column data
        [Symbol(t) for t in taxa]);        # column names

julia> first(avD_df, 2) # first 2 rows
2×16 DataFrame
 Row │ P_brandegeii  P_carneum  P_chartaceum  P_delicatum   ⋯
     │ Float64       Float64    Float64       Float64       ⋯
─────┼───────────────────────────────────────────────────────
   1 │    0.0        0.0124633    0.00970258    0.0103238   ⋯
   2 │    0.0124633  0.0          0.00855209    0.0123552
                                           12 columns omitted

julia> CSV.write("averagedist.csv", avD_df); # save to a csv file
```
"""
function averagepairwisedistances(D, ngenes)
  return sum(D) ./ ngenes
end
