# PhyloNetworks v0.16 to v1 refactoring: re-naming guide

This document lists the main changes to names in `PhyloNetworks` when
the package was refactored, in Nov. 2024.
It is intended to help users update scripts developed under `PhyloNetworks`
v0.16 or older, to work with v1 and later.

While [PhyloNetwork's documentation](https://juliaphylo.github.io/PhyloNetworks.jl/stable/)
should be help: **each package version** has its **own documentation**.
We can **switch** between documentations by clicking on the desired package version
in the "version" menu on the lower left corner.
This documentation should help to update code that uses *exported* functions,
especially using the search box (near the top left corner) when we don't know
the exact function name.
But changes in *internal* field names, in particular, are not typically documented,
so here they are.

## types and functions moved away

Many functions were moved to [`SNaQ.jl`](https://github.com/juliaphylo/SNaQ.jl).
Obviously, `snaq!` is one of them.

Internal utilities specific to quartets were moved also, including
types `DataCF` and `Quartet`, and methods that used these objects,
possibly renamed. Here is a non-exhaustive list.

| old name        | new name
| ----------------|---------
|`readTopologyLevel1`     | `readnewicklevel1`
|`readMultiTopologyLevel1`| `readmultinewicklevel1`
| `unionTaxa`     | `tiplabels`
| `unionTaxaTree` | `tiplabelsTree`
| `readTrees2CF`  | `readtrees2CF`

and these internals, not renamed:

`readInputData`,
`allQuartets`, `randQuartets`, `readListQuartets`,
`taxadiff`, `sameTaxa`,
`calculateObsCFAll!`, `calculateObsCFAll_noDataCF!`,
`taxaTreesQuartets`, `taxonTrees`, `taxonQuartets`, `descData`, `summarizeDataCF`
`updateBL!`,
`edgesParts`, `getDescendants!`, `makeTable`,
and more.

## renamed functions: CamelCase to lowercase

To follow [julia conventions](https://docs.julialang.org/en/v1/manual/style-guide/#Use-naming-conventions-consistent-with-Julia-base/),
most exported methods and many internal methods
had their names changed to be "lowercase" if they were not already.

### exported functions

An old name (typically "CamelCase") was deprecated in v0.1 if
it was exported in v0.16, to make it usable without breaking old code.
In many cases, the function name was also changed to be more descriptive.

Most notably:

old name             | new name
---------------------|---------
`tipLabels`          | `tiplabels`
`readTopology`       | `readnewick`
`writeTopology`      | `writenewick`
`writeMultiTopology` | `writemultinewick`
`readMultiTopology`  | `readmultinewick`

For the full list, see
[v1.0.0 deprecated.jl](https://github.com/JuliaPhylo/PhyloNetworks.jl/blob/v1.0.0/src/deprecated.jl).
Note that using a deprecated name does not trigger a deprecation message by default, so users may not
know that their scripts should be updated.
To see deprecation warnings, start julia with option `depwarn=yes` like this:

```sh
julia --depwarn=yes
```

### internal functions

Some methods that were internal in v0.16 were renamed to be "lowercase",
but not deprecated in v1 (because not exported, and not expected to be used in scripts).
Some are now exported in v1.

`blobInfo`, `displayedNetworks!`,  
`resetEdgeNumbers!`, `resetNodeNumbers!`,  
`sampleBootstrapTrees` → `samplebootstrap_multiloci`,  
`inheritanceWeight`,
`pairwiseTaxonDistanceMatrix!`,  
`pairwiseTaxonDistanceGrad` → `pairwisetaxondistance_gradient`  
`readFastaToArray` (exported in v1), `readCSVtoArray`  
`writeTableCF` (internal in v0.16) → `tablequartetCF` (exported in v1)  
`writeSubTree!` (exported, but should not: so not deprecated),  
`readSubtree!`           → `readnewick_subtree!`,  
`parseRemainingSubtree!` → `parsenewick_remainingsubtree!`,  
`readFloat`              → `readnewick_float`,  
`getDataValue!`          → `parsenewick_getfloat!`,  
`parseEdgeData!`         → `parsenewick_edgedata!`,  
`parseHybridNode!`       → `parsenewick_hybridnode!`,   
`readnodename`           → `readnewick_nodename`,  
`parseTreeNode!`         → `parsenewick_treenode!`,
`synchronizePartnersData!`,
`updatePostOrder!`       → `postorder_nodeupdate!`,  
`updatePreOrder!`        → `preorder_nodeupdate!`,  
`remove_edgeLengthsGammas!`,  
`solvePolytomyRecursive!` and `solvePolytomy!` → re-written and replaced by `resolvetreepolytomy!`

also many others in `parsimony.jl`.

Some keyword arguments were renamed to be lowercase also.
Notably: `checkPreorder` → `checkpreorder` or simply `preorder`.

## renamed internal fields

These field names are *very* internal! They should not be used. If they are,
document which package version is used to make the code reproducible,
and expect un-documented changes.
Use exported methods instead, such as `getroot` etc.

Edge's:

| old field name  | new field name
|-----------------|---------------
| `isChild1`      | `ischild1`
| `isMajor`       | `ismajor`
| `containRoot`   | `containroot`
| `inCycle`       | `inte1`
|`istIdentifiable`| `boole1`
|`fromBadDiamondI`| `boole2`

Node's:

| old field name  | new field name
|-----------------|---------------
| `gammaz`        | `fvalue`
| `hasHybEdge`    | `booln1`
| `isBadDiamondI` | `booln2`
| `isBadDiamondII`| `booln3`
| `isExtBadTriangle` | `booln4`
| `isVeryBadTriangle`| `booln5`
| `isBadTriangle`    | `booln6`
| `inCycle`       | `intn1`
| `k`             | `intn2`
| `typeHyb`       | `int8n3`

HybridNetwork's:

| old field name  | new field name
|-----------------|---------------
|`numTaxa`| `numtaxa`
| `numNodes`      | `numnodes`
| `numEdges`      | `numedges`
| `root`          | `rooti`
| `numHybrids`    | `numhybrids`
| `cladewiseorder_nodeIndex`| `vec_int1`
| `visited`       | `vec_bool`
| `edges_changed` | `vec_edge`
| `nodes_changed` | `vec_node`
| `ht`            | `vec_float`
| `numht`         | `vec_int2`
| `numBad`        | `intg1`
| `hasVeryBadTriangle`| `boolg1`
| `index`         | `vec_int3`
| `loglik`        | `fscore`
| `blacklist`     | `vec_int4`
| `cleaned`       | `boolg2`
| `isRooted`      | `isrooted` (don't trust it!)
