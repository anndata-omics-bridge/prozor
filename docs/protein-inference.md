# Protein inference

## Peptide--protein edges

Inference consumes ordinary `(peptide, protein)` string pairs. A set is the
natural boundary between occurrence matching and topology: repeated sites and
duplicate inputs collapse to one edge without constructing a numeric matrix.

```python
from prozor.inference.greedy import greedy_parsimony

edges = {
    ("PEP1", "P1"),
    ("PEP1", "P2"),
    ("PEP2", "P1"),
}
result = greedy_parsimony(edges)
```

Prozor needs only identity topology for parsimony. It therefore has no NumPy or
SciPy dependency and does not assign numeric weights to edges.

## Greedy parsimony

At each iteration the algorithm:

1. finds the active protein evidence covering the most unexplained peptides;
2. groups proteins with identical remaining evidence;
3. separates independent candidates into disjoint overlap components;
4. sends only overlapping, non-identical candidates to the tie resolver;
5. optionally includes proteins whose remaining evidence is a subset of the
   selected evidence; and
6. assigns the covered peptides before continuing.

Independent components are processed together. Final groups use the same
stable ordering as the sequential algorithm: descending peptide count,
descending protein-group size, then protein accession.

## Three different equal-score situations

- Proteins with identical evidence are one indistinguishable protein group;
  they are not a tie to resolve.
- Candidates with disjoint evidence are independent. Stable ordering is enough
  because choosing either first cannot change the other.
- Overlapping candidates with different evidence form a consequential tie.
  Only this case reaches the injected `TieResolver`.

```python
from collections.abc import Sequence

from prozor.inference.greedy import greedy_parsimony
from prozor.inference.ties import TieCandidate


def prefer_b(candidates: Sequence[TieCandidate]) -> TieCandidate:
    return next(candidate for candidate in candidates if candidate.proteins == ("B",))


result = greedy_parsimony(edges, resolve_tie=prefer_b)
```

The default preserves Prozor's deterministic group-size and accession rule.
Sequence coverage, database scores, and quantitative correlation are promising
experimental evidence, but are not silently applied by the production API.

## Subsumed proteins

With the default `subsume=True`, a protein supported only by a subset of the
selected group's peptide evidence is retained in that group. Set
`subsume=False` to report only proteins with the winning evidence signature.

```python
groups = greedy_parsimony(edges, subsume=False)
peptide_to_group = groups.to_dict()
```

## Interpretation boundary

Prozor returns inferred group membership and peptide assignments. It does not
decide whether a consumer should replace vendor-reported protein groups, mix
target and decoy evidence, quantify groups, or persist inference results. Those
are workflow policies and must be recorded by the consuming application.
