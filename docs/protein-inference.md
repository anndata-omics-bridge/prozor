# Protein inference

## Peptide--protein topology

`PeptideProteinMatrix` represents peptides as rows and proteins as columns in a
SciPy CSR matrix. Construct it directly from edges:

```python
from prozor.sparse_matrix import PeptideProteinMatrix

matrix = PeptideProteinMatrix.from_edges(
    [
        ("PEP1", "P1"),
        ("PEP1", "P2"),
        ("PEP2", "P1"),
    ]
)
```

Repeated occurrence sites collapse to one peptide--protein edge. Labels are
sorted deterministically, and matrix shape is checked against the labels.

Two weightings are available:

| Weighting | Meaning |
| --- | --- |
| `binary` | Every observed peptide--protein edge has weight 1. |
| `inverse` | Each peptide distributes total row weight 1 across its proteins. |

Topology methods such as `proteins_per_peptide`, `peptides_per_protein`, and
greedy inference count populated edges rather than summing weights. Inverse
weighting therefore does not change group topology.

## Greedy parsimony

```python
from prozor.greedy import greedy_parsimony

result = greedy_parsimony(matrix, subsume=True)
```

At each iteration the algorithm:

1. selects the active protein evidence covering the most active peptides;
2. groups proteins with identical remaining peptide evidence;
3. optionally includes proteins whose remaining evidence is a subset of the
   selected evidence;
4. assigns the covered peptides to that group; and
5. removes the assigned evidence before the next iteration.

Protein and peptide identifiers inside each group are sorted. Ties between
candidate groups are resolved deterministically from their identifiers.

## Subsumed proteins

With the default `subsume=True`, a protein supported only by a subset of the
selected group's peptide evidence is retained in that group. Set
`subsume=False` to report only proteins with the winning evidence signature.

```python
groups = greedy_parsimony(matrix, subsume=False)
peptide_to_group = groups.to_dict()
```

## Interpretation boundary

Prozor returns inferred group membership and peptide assignments. It does not
decide whether a consumer should replace vendor-reported protein groups, mix
target and decoy evidence, quantify groups, or persist inference results. Those
are workflow policies and must be recorded by the consuming application.
