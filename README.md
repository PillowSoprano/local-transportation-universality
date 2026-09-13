# Local universality: code and exact certificates

Code and exact certificates for **Local universality in three-way transportation polytopes**, by **Xiyao Yu**, Nanyang Technological University, Singapore.

[ORCID](https://orcid.org/0009-0006-2673-5269) · [GitHub](https://github.com/PillowSoprano)

Python 3.10+ standard library suffices. Run commands from this directory.

## Replay the supplied evidence

```bash
python3 code/verify_klein.py
python3 code/replay_unit_compiler.py
```

The first command verifies the 11-vertex, 22-triangle Klein-bottle example.
The second reads the arrays directly, rebuilds all pair rows, checks exact unit
margins and Q rank, then eliminates actual two-entry rows over the integers and
computes the reduced Smith form. It also verifies all full-array completions.
It does not import the construction code.

## Regenerate the construction and polygon examples

```bash
python3 code/construct_unit_torsion.py --orders 2 3 4 5 7 8 11 13 31 97 257 65537 --smith-limit 400
python3 code/verify_unit_compiler.py
python3 code/complete_unit_vertex.py data/unit_torsion_97.json
python3 code/complete_unit_vertex.py data/unit_compiler_pentagon.json --output data/full_unit_pentagon.json
python3 code/complete_unit_vertex.py data/taranenko/fractional_core_3.json --output data/full_unit_core_3.json
python3 code/replay_unit_compiler.py
```

The binary-family certificates cover twelve denominators through 65537.
The full completion in that family was computed for m=97, with order 644.
For m=65537, the partial array was computed and verified; the predicted
completion order 1472 was calculated, but that full completion was not run.
Additional full certificates are the denominator-three core in order 8 and
the pentagon face in order 494.

## File formats

Partial records list cells, a common denominator, and integer numerators.
Unlisted cells are zero. Every pair occurring in these supports has margin one.

Full records retain original cells and numerators and provide an
integer_completion_by_AB table: an entry c >= 0 denotes the added coordinate
(a,b,c) with weight 1; -1 denotes a pair served by the original support.
Every added column has three pair rows disjoint from all other columns.
This gives the full rank from the exact rank of the original support.

Polygon records store the integer-affine map using coordinate_variables and
coordinate_signs. A sign +1 gives t_i; -1 gives 1-t_i. Flag-encoded records
also specify the integer margins and their implied structural zeros.

## Sources

The constructions and verification code accompany the manuscript
Local universality in three-way transportation polytopes.
The completion uses the Latin-framework technique of Colbourn,
https://doi.org/10.1016/0166-218X(84)90075-1, as explained in Bartlett's
dissertation, Section 3.3:
https://thesis.caltech.edu/7819/1/caltech_dissertation_padraic_draft.pdf

The small denominator-three core derives from the public order-four catalogue
by Vladimirov, Taranenko, and Krotov, https://doi.org/10.5281/zenodo.14854809.
Its source is recorded in data/taranenko/PROVENANCE.json.
The catalogue is not itself reproduced in this archive.

The Klein-bottle triangulation is isomorphic to a configuration in Venturello's
balanced-complex library: https://doi.org/10.37236/8394.
