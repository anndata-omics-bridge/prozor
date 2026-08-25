# Matching peptides

## Aho--Corasick semantics

Prozor builds one automaton from the distinct peptide patterns and searches each
protein sequence in a single pass. It returns every occurrence, including:

- a peptide repeated at multiple positions in one protein;
- a short peptide nested inside a longer peptide; and
- overlapping matches.

Coordinates use Python's half-open convention: `start` is inclusive and `end`
is exclusive. The matched length is therefore `end - start`.

Duplicate input peptide patterns are removed while preserving their first-seen
order. Occurrences are not deduplicated; repeated sites remain separate
`PeptideAnnotation` records.

## Backend selection

Three backend requests are accepted:

| Request | Behavior |
| --- | --- |
| `auto` | Use `ahocorasick_rs` when installed, otherwise `ahocorapy`. |
| `ahocorapy` | Require the portable pure-Python backend. |
| `ahocorasick_rs` | Require the Rust backend; raise `ImportError` if unavailable. |

All public matching operations default to `auto`. Both implementations are
installed by Prozor, so a normal installation selects Rust. The pure-Python
implementation remains an explicit choice and the fallback when Rust cannot be
imported.

```python
from prozor.matching.annotation import annotate_peptides_streaming

result = annotate_peptides_streaming(
    ["PEPTIDE"],
    [("P1", "MPEPTIDEX")],
    backend="auto",
)

print(result.requested_backend)  # "auto"
print(result.resolved_backend)   # "ahocorapy" or "ahocorasick_rs"
```

Both implementations are required to return the same nested and overlapping
matches. The Rust adapter explicitly enables overlapping search.

## Case handling

Annotation functions are case-sensitive and do not normalize their inputs.
Consumers should supply peptide and protein strings in a consistent alphabet.
For lower-level use, `create_automaton(..., case_sensitive=False)` provides
case-insensitive matching while returning the original keyword spelling.

## Streaming contract

`annotate_peptides_streaming` accepts an iterable of `(protein_id, sequence)`
tuples. The iterable may be consumed only once, which allows direct integration
with FASTA readers and database cursors. The resulting occurrence annotations
are retained in memory.

Prozor deliberately does not own:

- FASTA parsing or header/accession interpretation;
- ProForma parsing or modified-to-unmodified sequence conversion;
- decoy and contaminant classification;
- enzyme settings or identification filtering; or
- persistence into analysis containers.

## Optional N-terminal tryptic context

`annotate_peptides(..., filter_tryptic=True)` keeps matches whose start is the
protein N-terminus, follows an initiator methionine, or follows one of the
configured prefix residues (by default `R` or `K`). This is an N-terminal
context filter, not a complete in-silico digestion model.
