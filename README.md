# Turán densities of complete 4-graphs

Exact certificates and independent checkers accompanying Bruno Mérola Corrêa's
*Turán densities of complete 4-graphs via flag algebras on seven vertices* (preprint, 2026).

**[Read the paper](paper/turan4.pdf)** · **[All results](RESULTS.md)** ·
**[Verification guide](VERIFICATION.md)** · **[Zenodo archive](https://doi.org/10.5281/zenodo.23245169)**

## Main results

Using flag algebras on seven vertices, the paper establishes:

- **Turán density:** $\pi(K_5^{(4)}) < 0.691600$, with an exact rational certificate.
- **Codegree-squared density:** $\sigma(K_5^{(4)}) = 31/64$, resolving a conjecture of
  Balogh, Clemen and Lidický, together with a stability theorem.

It also gives further hypergraph bounds, limits of the plain flag algebra method on six and seven vertices,
and applications to lottery numbers. The [results index](RESULTS.md) links each result to its certificates
and checks. Exact fractions are the certified bounds; displayed decimals are rounded in the safe direction.

## Run the checks

Download the complete package. You need Python 3.11 or newer, NumPy 2.0 or newer, SciPy, Clarabel,
GCC with OpenMP, Bash, and `sha256sum` or `shasum`.
See the [setup notes](VERIFICATION.md#how-to-verify) for interpreter and compiler options, including macOS.

From the package root, start with:

```bash
verifier/verify_all.sh --fast --threads 2
```

**This is not the full verification:** it runs the 161 fast checks but omits the long raw scans.
To complete the independent checks supplied in the package, also run:

```bash
verifier/verify_all.sh --raw all --threads 2
```

The raw scans can take many hours. Both commands must finish with exit status 0.
The driver prints PASS/FAIL for each step and saves detailed output in `work/logs/`.
The [verification guide](VERIFICATION.md) explains individual checks, measured timings and memory use;
[EXPECTED.txt](verifier/EXPECTED.txt) records expected outputs. The search code is not needed for verification.

## What the verification establishes

The checkers test the certificates in exact arithmetic and, where required, enumerate every admissible graph
using the raw scans. They were implemented separately from the search code and share no code with it.
All fast checks and raw scans were run on clean cloud machines before the initial public releases;
the author also ran the fast checks.

These are computational checks, not a formal verification of the paper in a proof assistant.
The argument also depends on the mathematical reductions in the paper and the correctness of the programs,
compiler and hardware. The [full trust model](VERIFICATION.md#trust-model) describes the scope of each check,
external mathematical dependencies and limitations.

## Package contents

| Location | Contents |
|---|---|
| [paper/](paper/) | Preprint PDF; arXiv identifier pending |
| [certificates/](certificates/) | Exact certificates and dual points, with documentation for each result |
| [verifier/](verifier/) | Independent checkers, driver and expected outputs |
| [search/](search/) | Code used to produce the certificates |
| [SHA256SUMS](SHA256SUMS) | File integrity manifest |

Program adaptations and package changes are recorded in [verifier/CHANGES.txt](verifier/CHANGES.txt)
and [search/CHANGES.txt](search/CHANGES.txt).

## AI assistance

Code was written with AI coding agents using several models, including Claude and GPT.
The verifiers were written in separate agent sessions that did not share code with the search code.
The paper describes AI assistance in the broader work.

## Citation and licences

Please cite the paper and the version of the data package used. Citation metadata is in
[CITATION.cff](CITATION.cff); archived releases and their DOIs are available on
[Zenodo](https://doi.org/10.5281/zenodo.23245169).

Code is licensed under [MIT](LICENSE-CODE); certificates and documentation under [CC BY 4.0](LICENSE-DATA).
