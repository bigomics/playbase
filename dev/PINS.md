# Pinned versions in the Omics Playground image chain

Every deliberate version pin across `playbase` and `omicsplayground`, why it exists,
and what has to happen before it can be dropped.

**Rule of thumb: do not add a pin unless the package is unobtainable otherwise.**
A pin that merely says "the version we happened to test" rots silently and gets
overridden downstream (see *Removed pins* below — most of the old ones did nothing).

## How dependency resolution actually works

Two different code paths, easy to confuse:

| Context | What drives installation |
|---|---|
| **Docker build** (`R/` not copied, `use.remotes=TRUE`) | `DESCRIPTION` `Imports:` + `Remotes:` via `remotes::install_deps('.')` |
| **Local dev** (`R/` present, `use.remotes=FALSE`) | `scan_packages()` in `dev/functions.R` (`remotes.url`, `add_github()`) |

So **`DESCRIPTION: Remotes:` is the operative pin list for the images.**
The lists in `dev/functions.R` only affect local developer installs.

`Remotes:` overrides *where* a package comes from — it applies to direct and
transitive dependencies alike, which is why a package can be pinned there
without appearing in `Imports:`.

## Load-bearing pins — DO NOT REMOVE

> **A `Remotes:` entry is not only a version pin — for a package that is not in
> `Imports:` it is the *only thing that installs the package at all*.** Delete it
> and the package silently disappears from the image; nothing errors at build
> time, the feature just breaks at runtime. This was confirmed the hard way on
> 2026-07-15 (see below). Always check `Imports:` **and** `grep -r 'pkg::' R/`
> before removing an entry.

These packages no longer exist upstream. Unpinning breaks the build.

| Package | Pin | Why | Droppable when |
|---|---|---|---|
| `KEGG.db` | `3.2.4` (BioC **3.11** tarball) | removed from Bioconductor after 3.11 | code stops using KEGG.db |
| `DeconRNASeq` | `1.44.0` (BioC **3.18** tarball) | removed from Bioconductor after 3.18 | code stops using it, or it returns to BioC |
| `grr` | `0.9.1` (CRAN Archive) | archived from CRAN | it returns to CRAN, or is vendored |

These two *look* like stale pins but are load-bearing **dependency declarations** —
neither is in `Imports:`, so the entry is what installs them:

| Package | Pin | Used by | Why still pinned |
|---|---|---|---|
| `infercnv` | `@infercnv-v1.3.3` | `R/pgx-cna.R` | not in `Imports:`; removing the pin removed the package entirely. Current BioC is **1.26.0** — a ~7-year API jump, so moving to it needs the CNA features actually tested, not just a green build. |
| `org.Pf.plasmo.db` | `3.14.0` (BioC 3.14) | `R/pgx-annot.R` | not in `Imports:`; same story. Present in current BioC annotation as 3.22.0. |

The clean fix for both is to declare them in `Imports:` and drop the pin, so they
resolve from the current repos — but that is a behaviour change that needs testing
of CNA analysis and Plasmodium annotation. Do it as its own task.

`DeconRNASeq` has a trap: `remotes` does **not** resolve dependencies of `url::`
remotes, so its dependency `limSolve` is installed explicitly in
`install_dependencies()`. If you ever repin DeconRNASeq, keep that line.

## Deliberate forks / tags — intentional, revisit occasionally

| Package | Pin | Why |
|---|---|---|
| `Seurat` | `satijalab/seurat@fix/v.5.3.1` | fork branch (`Dockerfile.update`) |
| `inspectdf` | `alastairrushworth/inspectdf@bugfix/cran-dplyr-failing` | upstream CRAN build broken (`Dockerfile.update`) |
| `visNetwork` | `bigomics/visNetwork` | bigomics fork (`Dockerfile.update`) |
| `RaMP` | `ncats/RaMP-DB@v3.0.2` | tag pin (`DESCRIPTION`) |
| `RefMet` | `metabolomicsworkbench/RefMet@v2` | tag pin (`DESCRIPTION`) |

Each should be re-checked when upstream releases: if the fix landed upstream, drop the fork.

## Platform / toolchain pins

| Thing | Pin | Where | Notes |
|---|---|---|---|
| Bioconductor | **3.22** | `Dockerfile.rbase`, `install_playbase.R`, `functions.R` | pinned in **3 places** — change all three together. Tied to R 4.5 (trixie); BiocManager refuses 3.18 on R 4.5. |
| CRAN snapshot | `__linux__/trixie/latest` | `dev/rspm.R`, `omicsplayground/dev/Rprofile` | the `<distro>` segment **must** match the base image, or PPM serves binaries built against the wrong distro. Needs `options(HTTPUserAgent=...)` in the same file or PPM silently serves source. |
| quarto CLI | **1.8.27** | `omicsplayground/docker/Dockerfile` **and** `docker/Dockerfile.update` | keep both in sync. Need >=1.8 for typst >=0.12. tinytex installs on top of the base one, so the base block cannot simply be deleted. |
| R | 4.5.0 | implied by `debian:trixie` | not pinned explicitly; comes from the distro. |

## Unpinned and floating — the real reproducibility gap

These track `HEAD`/`master` with no pin, so a rebuild silently picks up whatever
upstream became. There is **no `renv.lock`** in either repo.

`PCSF` · `playdata` · `EPIC` · `SuperCell` · `NNLM` · `azimuth` · `plaid` ·
`metaLINCS` · `grimon`

Deciding whether to pin these (or add a lockfile) is an open question — it is a
bigger risk than any stale version pin.

## Removed pins (2026-07-15) — kept here so nobody "restores" them

All of these were **already being overridden** later in the build; the shipped
image never contained the pinned version, so removing them changed nothing except
removing wasted work.

Each was verified by rebuild: the package is still installed, at the version the
image was *already* shipping. All are transitive dependencies of something else in
`Imports:`, which is why removing their entry does not remove them.

| Package | Was pinned to | Still installed as | Why removal is safe |
|---|---|---|---|
| `rms` | `6.8-0` (Archive) | **8.1-1**, as a PPM **binary** | build compiled 6.8-0 from source, then threw it away and installed 8.1-1 — pure waste |
| `rjson` | `0.2.21` (Archive) | **0.2.23** binary | pin never applied |
| `rliger` | `@v2.1.0` | **2.2.1** | in `Imports:`, resolves from CRAN |
| `BiocManager` | `1.30.23` | **1.30.27** | `install_playbase.R` runs `update.packages()`, which bumped it anyway. The pin also made the `BiocManager::valid()` build gate permanently unsatisfiable. |

`force.remotes.url` in `dev/functions.R` held `rjson`/`rms` and force-installed them
on every build regardless of need; it is gone.

### Why removing dead pins is worth doing

A `url::` pin always installs from a tarball, which means a **source build** — it
bypasses the Posit binary repo entirely. Dropping the `rms` pin alone let it arrive
as a prebuilt binary. Measured on the `playbase-pkg` layer:

| | before | after |
|---|---|---|
| `install_playbase.R` step | **8567 s** (143 min) | **2243 s** (37 min) |

So stale pins are not merely cosmetic — they cost build time by forcing source
compiles of packages that PPM already ships as binaries.
