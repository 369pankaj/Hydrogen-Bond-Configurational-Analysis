# H-Bond Strength Pipeline

Hydrogen-bond **free energies** (strengths) computed directly from GROMACS
trajectories with the reversible-work / PMF method of **Sapir & Harries**
(*J. Chem. Theory Comput.* **2017**, *13*, 2851) and the modified convention of
**Olgenblum et al.** (*JCTC* **2020**, *16*, 1249).

Developed and **validated on 4 M urea/water** (reproduces the committed Fortran
reference grids to 4 decimal places), and generalized through a chemistry
registry to arbitrary **binary osmolyte–water mixtures** (TMAO, DMSO, betaine, …),
**ternary mixtures**, and **protein-containing systems**.

---

## Table of contents

- [Contents](#contents) — files in this repository
- [Theory](#theory) — the full method, from (r, θ) to `dG`
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Configuration — `config.sh`](#configuration--configsh)
- [Generalizing to other mixtures & proteins](#generalizing-to-other-mixtures--proteins)
- [Running on a cluster (SLURM)](#running-on-a-cluster-slurm)
- [Using the tools directly](#using-the-tools-directly)
- [Output grid format](#output-grid-format)
- [Validation](#validation)
- [Notes & caveats](#notes--caveats)
- [References](#references)

---

## Contents

| File | Purpose |
|------|---------|
| `hb_pmf_traj.py` | Build the 6-column `P(r,θ) / P_rand` grid for a donor–acceptor pair from an MD trajectory (MDAnalysis). |
| `hb_strength.py` | Compute the H-bond free energy `dG` (Eq. 2) from a grid file; optional PMF landscape plot. |
| `validate.py` | Regression test: reproduces the reference 4 M urea/water values from committed grids. |
| `molecules.yaml` | **Chemistry registry** — donor/acceptor atoms per molecule (water, urea, TMAO, DMSO, betaine, protein). |
| `systems.yaml` | **Job list** — systems, their components, and run defaults for the automated front-end. |
| `hb_run.py` | **Automated driver** — auto-generates every valid donor→acceptor species from the registry and runs the pipeline. Use for arbitrary binary/ternary/protein systems. |
| `config.sh` | Simpler bash config (systems + explicit species table) for the hand-curated urea/water case. |
| `run_hb.sh` | Bash driver over `config.sh`. Equivalent to `hb_run.py` for the fixed urea/water species. |
| `submit.slurm` | SLURM submission *template* for HPC clusters. |
| `results/` | Auto-created; per-system grids, logs and CSV summaries. |

There are **two front-ends**: `hb_run.py` (registry-driven, generalizes to any
mixture/protein — recommended) and `run_hb.sh` (bash, the original hand-listed
urea/water species). Both call the same validated `hb_pmf_traj.py` +
`hb_strength.py` core.

---

## Theory

> **TL;DR.** A hydrogen bond is treated as a *statistical correlation* between two
> molecules in a specific geometry, not as a bond you count. We measure how much
> more likely a donor–acceptor pair is to be found in H-bonding geometry than it
> would be by chance, and convert that likelihood ratio into a free energy via the
> reversible work of forming the observed distribution from a random one
> (Sapir & Harries, *JCTC* 2017). The result `dG` (kJ/mol) is the **strength** of
> the bond; it is negative for a favourable interaction.

### 1. The pair reaction coordinate (r, θ)

Every ordered donor(heavy)–acceptor(heavy) pair is described by two numbers:

```
r  = |D − A|                                   distance between the heavy atoms (Å)
θ  = min over the donor's hydrogens H of  ∠(D→A , D→H)      (degrees)
```

`θ` is the **deviation from linearity**: when a donor hydrogen points straight at
the acceptor, the vectors `D→A` and `D→H` are parallel and `θ → 0°`. Taking the
*minimum* over the donor's hydrogens selects the best-aligned H, i.e. the one
actually capable of bonding. Distances use the minimum-image convention so
periodic boundaries and molecule wrapping are handled correctly.

For **like–like pairs** (water–water) the roles are symmetric: molecule A can
donate to B *or* accept from B. `--symmetric` mode therefore also evaluates the
acceptor's hydrogens along the reversed axis `A→D` and keeps the smaller `θ`, so a
pair counts as H-bonded if *either* partner donates. Each unordered pair is
processed once (no double counting). For **directional pairs** (water→urea-O,
urea-N→water, water→TMAO-O …) only the true donor's hydrogens are used.

### 2. Observed density  P(r, θ)

Accumulating `(r, θ)` for all pairs over all frames and normalising the 2-D
histogram gives the joint probability **density**

```
P(r,θ)  ≥ 0 ,        ∬ P(r,θ) dr dθ = 1
```

`P` concentrates in a sharp ridge at the H-bonding distance (~2.8 Å for O···O) and
small `θ`. It is the empirical fingerprint of the interaction.

### 3. Random reference  P_rand(r, θ)

To know whether that ridge is *significant*, compare it to the density expected if
the two molecules were uncorrelated (ideal gas). Position and orientation
decouple in the reference, so it **factorises**:

```
P_rand(r,θ) = P_rand(r) · P_rand(θ)
```

- **Radial part** — a uniform density of acceptors around the donor puts
  probability in a shell ∝ its volume:

  ```
  P_rand(r) ∝ 4π r²           (code: shell-volume ratio (r_hi³ − r_lo³)/r_max³, per Å)
  ```

- **Angular part** — for the angle between two vectors, a random orientation is
  distributed by the solid angle, `P_rand(θ) ∝ sin θ`. This analytic form is the
  default; a tabulated `--prand-theta` file can be supplied to match a specific
  reference convention. **The *shape* of `P_rand(θ)` does not cancel later** (§6),
  so this choice matters.

### 4. Potential of mean force  (Eq. 1)

The ratio of observed to random density is a free-energy surface — the potential
of mean force along `(r, θ)`:

```
PMF(r,θ) = −RT · ln[ P(r,θ) · F_avg / P_rand(r,θ) ]
```

`RT` is the thermal energy (`R = 8.314 × 10⁻³ kJ/mol/K`; at 300 K, RT ≈ 2.494
kJ/mol). Wherever the pair is *more* likely than random the PMF is negative (a
well); the H-bond minimum is the deepest point, and `dG_PMF = min PMF` is the
Eq. 1 estimate of the bond strength. `F_avg` is an optional bulk correction from
the radial distribution `g(r)` (supplied with `--rdf`) that references the PMF to
the bulk liquid; it shifts the PMF column but, as shown next, **cancels** in the
Eq. 2 strength.

### 5. H-bond strength as reversible work  (Eq. 2 — the reported value)

The quantity this pipeline reports as `dG` is the **reversible work of forming the
observed pair distribution from the random one, within the first solvation
shell** (Sapir & Harries 2017, Eq. 2). Restrict both densities to the shell
`r ∈ [r_min, r_max]` and renormalise them there:

```
ζ      = 1 / ∬_shell P(r,θ)      dr dθ           N₂(r,θ) = ζ · P(r,θ)
ζ_rand = 1 / ∬_shell P_rand(r,θ) dr dθ           N₁(r,θ) = ζ_rand · P_rand(r,θ)
```

so `N₂` and `N₁` are proper probability densities over the shell. The strength is

```
dG = −RT ∬_shell  N₂ · ln( N₂ / N₁ )  dr dθ  =  −RT · D_KL( N₂ ‖ N₁ )
```

This is `−RT` times the Kullback–Leibler divergence of the observed shell
distribution from the random one. Because `D_KL ≥ 0`, **`dG ≤ 0`**: the more the
real geometry departs from random, the more negative (stronger) the bond.

**Why the bulk `g(r)`/`F_avg` cancels.** `F_avg` multiplies `P` by a constant.
`N₂ = ζP` with `ζ = 1/∬P` is invariant under any constant rescaling of `P`, so the
strength is identical whether or not the bulk correction was applied. That is why
`g(r)` is needed only for the *PMF column* (§4), never for the strength — an RDF
is optional for `dG`.

### 6. Two conventions

```
dG_2017 = −RT ∬ N₂ ln(N₂/N₁)      forward KL  (Sapir & Harries 2017; the reported dG)
dG_2020 = −RT ∬ N₁ ln(N₁/N₂)      reverse KL  (Olgenblum et al. 2020)
```

They weight the same ratio differently — 2017 weights by the *observed* density
(emphasising populated H-bond geometries), 2020 by the *random* density. Both are
printed with `--all`; `dG_2017` is the pipeline's headline number.

### 7. From histogram to number (the estimator)

In practice the integrals are trapezoidal sums over the grid bins. Empty observed
bins contribute nothing (`lim_{x→0} x ln x = 0`) and empty reference bins are
dropped, so a finite grid gives a finite `dG` (reported as a lower bound if any
bin was dropped). Fully-populated grids — like the committed Fortran references —
have every bin positive and reproduce bit-for-bit.

### 8. Assumptions & limitations

- **First-shell definition.** `dG` depends on the shell window `[r_min, r_max]`;
  it should bracket the first `g(r)` peak (contact → first minimum). The reference
  used `[2.49, 3.25]` Å for O···O; `shell: auto` estimates it per pair.
- **Angular reference.** As in §3, the `P_rand(θ)` shape enters the log directly
  and affects `dG`; keep it consistent across systems you compare.
- **Geometric H-bond, not energetic.** The method uses geometry (r, θ) only; it
  never inspects the force field's electrostatics. The "strength" is the
  free-energy signature of the geometric correlation, in the spirit of the cited
  papers, not an enthalpic bond energy.
- **Two-body.** Cooperativity beyond the pair correlation is not modelled.

---

## Requirements

- Python ≥ 3.8
- [`MDAnalysis`](https://www.mdanalysis.org/) and `numpy`
- `matplotlib` (only for `hb_strength.py --plot`)

```bash
pip install MDAnalysis numpy matplotlib
```

---

## Quick start

```bash
# 1. sanity check against the committed reference grids
python validate.py
#   -> ALL VALUES REPRODUCED (PASS)
```

**Recommended front-end — `hb_run.py`** (registry-driven; works for any
mixture/protein defined in `molecules.yaml` + `systems.yaml`):

```bash
./hb_run.py --dry-run 4M     # print the species it will generate + commands
./hb_run.py 4M               # run one system
./hb_run.py                  # run every system listed in systems.yaml
./hb_run.py --stride 100 4M  # override a default for this invocation
```

**Bash front-end — `run_hb.sh`** (equivalent for the fixed urea/water species,
no PyYAML needed):

```bash
./run_hb.sh 4M               # one system
./run_hb.sh                  # all systems (or: ./run_hb.sh all)
STRIDE=100 ./run_hb.sh 4M    # env override
```

Either way, results land in `results/<system>/`:

```
results/4M/
├── <species>_grid.dat            # full 6-column grid, per species
├── <species>_grid_1st_shell.dat  # first-shell slice used for dG
├── strengths_4M.csv              # dG_2017, dG_2020, dG_PMF per species
└── run.log
```

Species names depend on the front-end: `hb_run.py` generates them from the
registry (`water_sym`, `waterD_ureaA`, `ureaD_waterA`, `urea_self`), while
`run_hb.sh` uses the hand-listed `config.sh` names (`WW`, `OuW`, `NuW`).

---

## Configuration — `config.sh`

All systems, species and defaults live in `config.sh`. Every default is also
overridable from the environment (`STRIDE=100 ./run_hb.sh 4M`).

**Add a system** — one line, `name="topology|trajectory"`, then add the name to
`SYSTEM_ORDER`:

```bash
declare -gA SYSTEMS=(
  [4M]="/path/md.gro|/path/md.xtc"
  [my_run]="/path/to/my.tpr|/path/to/my.xtc"     # <-- new
)
SYSTEM_ORDER=(4M 4M_500 4M_1000 my_run)
```

**Add a species** — `name="donor_heavy | acceptor | prand_file | extra_flags"`:

```bash
declare -gA SPECIES=(
  [WW]="resname SOL and name OW|resname SOL and name OW|P_rand_WW.dat|--symmetric --exclude-same-residue"
  [OuW]="resname SOL and name OW|resname UREA and name OU|P_rand_WH.dat|"
  [NuW]="resname UREA and name N1U N2U|resname SOL and name OW|P_rand_WW.dat|"
)
SPECIES_ORDER=(WW OuW NuW)
```

| key | meaning |
|-----|---------|
| `donor_heavy` | MDAnalysis selection for the donor **heavy** atom(s) |
| `acceptor` | MDAnalysis selection for the acceptor heavy atom(s) |
| `prand_file` | random-angle reference, looked up in `$PRAND_DIR` (default `..`) |
| `extra_flags` | passed verbatim to `hb_pmf_traj.py` — use `--symmetric --exclude-same-residue` for like–like pairs; empty for directional pairs |

Defaults (override via env): `STRIDE=200`, `SHELL_MIN=2.49`, `SHELL_MAX=3.25`,
`TEMPERATURE=300`, `PYTHON=python`, `PRAND_DIR=..`, `RESULTS_DIR=./results`.

---

## Generalizing to other mixtures & proteins

The registry-driven front-end (`hb_run.py`) automates the pipeline for **any**
binary osmolyte–water mixture (TMAO, DMSO, betaine, …), ternary osmolyte
mixtures, and protein-containing systems — you describe the *chemistry* once and
the *systems* once, and it enumerates the donor→acceptor species itself.

**1. Describe each molecule's H-bond chemistry** in `molecules.yaml`:

```yaml
molecules:
  tmao:
    resname: TMAO
    donors:    []        # acceptor-only  -> water donates TO it, never the reverse
    acceptors: [O]
  urea:
    resname: UREA
    donors:    [N1U, N2U]
    acceptors: [OU]
  water:
    resname: SOL
    donors:    [OW]
    acceptors: [OW]
    like_like_symmetric: true   # water self-pair is geometrically symmetric
```

**2. List systems and their components** in `systems.yaml`:

```yaml
systems:
  ternary_protein:
    topology:   /path/md.tpr     # .tpr/.psf (bonds) required for proteins
    trajectory: /path/md.xtc
    components: [water, urea, protein]
    exclude: molecule            # protein self-pairs excluded by chain, not residue
    neighbor_list: true          # capped search for large solvated systems
    shell: auto                  # per-pair first shell from the RDF first minimum
```

**3. Run** (dry-run first to see the plan):

```bash
./hb_run.py --dry-run ternary_protein   # prints the generated species + commands
./hb_run.py ternary_protein             # executes
```

### How species are generated

For a component list, `hb_run.py` emits:

| pair | result |
|------|--------|
| like–like, `like_like_symmetric` (water) | **one symmetric** species (both partners' H's) |
| same species, directional (urea) | **one self** species, donor→acceptor, intramolecular removed |
| cross species, each ordering | **one directional** species each → e.g. `waterD_ureaA` **and** `ureaD_waterA` |
| acceptor-only molecule (TMAO/DMSO/betaine) | only ever on the **acceptor** side |

So `components: [water, urea]` reproduces exactly the four validated species
(WW, Ou–W, Nu–W, U–U). `components: [water, tmao]` yields only `water→tmao` and
water–water — correct, since TMAO cannot donate.

### Two engine flags added for the general case

- `--exclude-same-molecule` — excludes intramolecular pairs by **bonded
  fragment** rather than residue. Essential for proteins, where one molecule
  spans many residues. Requires a topology with bonds (`.tpr`/`.psf`/`.pdb`+CONECT);
  errors clearly on a bare `.gro`. (For single-residue molecules it is identical
  to `--exclude-same-residue` — verified byte-for-byte.)
- `--neighbor-list` — replaces the full donor×acceptor distance matrix with a
  capped neighbour search. **Identical result** (verified byte-for-byte); much
  faster on large solvated systems, no benefit on small ones.

### What still needs your judgement (not auto-magic)

- **Chemistry is data you must get right once.** The atom names in
  `molecules.yaml` must match your force field; the acceptor-only vs
  donor/acceptor classification is real chemistry (TMAO/DMSO/betaine donate no
  H-bonds — do **not** set `--symmetric` for them).
- **First-shell window.** `shell: auto` reads the pair RDF (first peak → first
  minimum); it gave `(2.60, 3.13)` for water–water here (cf. the reference
  `[2.49, 3.25]`). Check it per pair, or set an explicit `shell:`.
- **P_rand(θ).** Unlike the radial `g(r)` (which cancels in dG), the *shape* of
  P_rand(θ) does **not** cancel. The default is analytic `sin θ`; set
  `prand_theta:` to a file to match a specific reference convention.
- **Protein numbers are unvalidated.** There is no reference protein grid to
  check against, so that tier is physically reasonable but carries no
  regression guarantee — benchmark before trusting it.

## Running on a cluster (SLURM)

The production machine here has **no scheduler**, so `run_hb.sh` is the primary
runner. On an HPC cluster, edit the `#SBATCH` lines and the module/venv block in
`submit.slurm`, then:

```bash
sbatch submit.slurm                          # all systems, one job
sbatch --export=ALL,STRIDE=100 submit.slurm  # override a default
# or enable "#SBATCH --array=0-2" + Option B for one system per array task
```

---

## Using the tools directly

**Build a grid** (arguments beyond the ones shown have sensible defaults — run
with `-h`):

```bash
python hb_pmf_traj.py \
    --top md.gro --traj md.xtc \
    --donor-heavy "resname SOL and name OW" \
    --acceptor    "resname SOL and name OW" \
    --symmetric --exclude-same-residue \
    --prand-theta ../P_rand_WW.dat \
    --stride 200 --temperature 300 \
    --shell 2.49 3.25 \
    -o WW_grid.dat
```

Key flags: `--symmetric` (use acceptor H's too), `--exclude-same-residue` (drop
intramolecular pairs, resolved on the *global* residue index so it is safe across
species), `-T/--temperature` (PMF column), `--rmax/--dr/--dtheta/--theta-max`
(binning), `--rdf` (optional `g(r)` for the PMF bulk correction only),
`--start/--stop/--stride` (frame range).

**Compute strengths** from one or many grids:

```bash
python hb_strength.py results/4M/*_1st_shell.dat -T 300 --all --csv out.csv
python hb_strength.py FULL_grid.dat --first-shell 2.49 3.25   # slice shell first
python hb_strength.py WW_grid_1st_shell.dat --plot pmf.png    # PMF landscape
```

---

## Output grid format

Six whitespace-separated columns, **r-major** (θ varies fastest):

```
r[Å]   θ[deg]   P(r,θ)   P_rand(r,θ)   PMF(r,θ)[kJ/mol]   P·PMF
```

`hb_strength.py` enforces r-major ordering on load and refuses to silently
transpose a θ-major file. Empty bins carry a finite (least-favourable) PMF so the
grid stays parseable by external/Fortran consumers; `dG_PMF = min(PMF)` is
unaffected.

---

## Validation

`validate.py` confirms the core reproduces the committed reference values within
`1e-3` kJ/mol:

| species | computed | reported |
|---------|---------:|---------:|
| W–W (water–water)        | −3.4941 | −3.4941 |
| U–W (urea–water, all)    | −1.3816 | −1.3816 |
| Ou–W (carbonyl O–water)  | −4.1218 | −4.1218 |
| Nu–W (amine N–water)     | −1.1155 | −1.1155 |

**End-to-end check** (grid rebuilt from `4M/md.xtc`, symmetric W–W, stride 500):
`dG = −3.467 kJ/mol` vs the committed Fortran grid `−3.494` — agreement to
≈0.03 kJ/mol, the residual being finite-frame sampling on the fine grid. Smaller
`STRIDE` tightens this.

---

## Notes & caveats

- **Temperature must match.** `hb_pmf_traj.py --temperature` sets the PMF column;
  `hb_strength.py -T` sets `dG` and reads that column for `dG_PMF`. `run_hb.sh`
  passes one `TEMPERATURE` to both so they always agree — keep them equal if you
  invoke the tools by hand.
- **RDF is optional.** `g(r)` cancels in `dG`; supply `--rdf` only if you want a
  bulk-referenced PMF column (`rdf_*_ow.dat` files live in `..`).
- **Symmetric mode counts each pair once** (`a_idx > donor`) and uses both
  partners' hydrogens — no double counting.
- **PBC & topology:** donor→H bonds are found within each molecule
  (`residue.atoms`) using the minimum-image convention, so wrapped molecules and
  `.gro` topologies with reused resids are handled correctly.

---

## References

1. A. Sapir, D. Harries, *Revisiting Hydrogen Bond Thermodynamics in Molecular
   Simulations*, **J. Chem. Theory Comput.** 2017, 13, 2851–2857.
   [doi:10.1021/acs.jctc.7b00238](https://doi.org/10.1021/acs.jctc.7b00238)
2. G. Olgenblum et al., **J. Chem. Theory Comput.** 2020, 16, 1249–1262.
3. Internal note, G. Olgenblum (2023), followed by this project.
