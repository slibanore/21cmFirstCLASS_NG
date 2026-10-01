# 21cmFirstCLASS_NG

Local primordial non-Gaussianity in [`21cmFirstCLASS`](https://github.com/jordanflitter/21cmFirstCLASS).

**Authors:** Sarah Libanore, Karin Finish

This code adds local-type primordial non-Gaussianity, parametrised by
$f_{\rm NL}$, to `21cmFirstCLASS`. It enters in three independent places — the
initial density field, the conditional halo mass function, and the
unconditional halo mass function — so that each contribution to the 21-cm
signal can be isolated.

It is a fork of `21cmFirstCLASS` with the full commit history preserved. The
`fnl` branch of the upstream repository contains the same work in its
development form; this repository is the cleaned, documented version that
accompanies the paper.

> **Accompanying paper:** Finish, Libanore et al. (arXiv:2610.xxxxx)

---

## Please read this before comparing against `21cmFirstCLASS`

**This code runs in double precision.** Upstream `21cmFirstCLASS` uses `float`
for the density fields and single-precision FFTW (`fftwf_`); here those are
`double` and `fftw_`. A run with `F_NL = 0`, therefore, does not
reproduce published `21cmFirstCLASS` results bit for bit. The FFT precision
differs and the random number stream is consumed differently when building the
initial conditions. The physics is the same and the differences are small, but
they are not zero, so the Gaussian baseline for any comparison must be
generated with *this* code rather than taken from upstream.

The motivation is the non-Gaussian calculation itself. The correction to the
mass function involves a near-cancellation between two comparable terms, 
which amplifies rounding error by about an order of magnitude.

---

## What is new

### Parameters

All default to the Gaussian behaviour, so adding them changes nothing until you
set `F_NL`. Full descriptions are in the `CosmoParams` and `UserParams`
docstrings in `src/py21cmfast/inputs.py`.

**`CosmoParams`**

| Parameter | Default | Meaning |
|---|---|---|
| `F_NL` | `0.` | Amplitude of local non-Gaussianity, CMB (Planck) convention |
| `KCUT_FNL` | `1e-3` | Lower cut-off $k_{\rm cut}$ [1/Mpc] on every leg of the bispectrum |

`F_NL` follows the Planck convention,

$$\zeta = \zeta_G + \tfrac{3}{5} f_{\rm NL}\left(\zeta_G^2 - \langle\zeta_G^2\rangle\right)
\quad\Longleftrightarrow\quad
\Phi = \Phi_G + f_{\rm NL}\left(\Phi_G^2 - \langle\Phi_G^2\rangle\right)$$

with $\Phi = \tfrac{3}{5}\zeta$ the primordial potential, so values are
directly comparable to [Planck 2018 IX](https://arxiv.org/abs/1905.05697).

**`UserParams`** — where the non-Gaussianity acts:

| Flag | Default | Effect |
|---|---|---|
| `NON_GAUSS_IC` | `False` | Non-Gaussian initial density field |
| `NON_GAUSS_FCOLL_COND` | `False` | Conditional halo mass function (`dNdM_conditional`) |
| `NON_GAUSS_FCOLL_UNCOND` | `False` | Unconditional halo mass function (`dNion_General`) |

Which approximation is used:

| Flag | Default | Alternative |
|---|---|---|
| `USE_EDG_uncond_hmf` | Edgeworth, [arXiv:2009.01245](https://arxiv.org/abs/2009.01245) Eqs. 38–40 | saddlepoint |
| `USE_LD_cond_hmf` | Lidz / D'Aloisio closed form | two-scale saddlepoint |
| `NG_MODEL_APPROX` | high-barrier limit, [arXiv:1304.8049](https://arxiv.org/abs/1304.8049) | full, [arXiv:1206.3305](https://arxiv.org/abs/1206.3305) |

Numerical controls:

| Flag | Default | Effect |
|---|---|---|
| `EXTRA_DIM_FNL` | `1.5` | De-aliasing pad for squaring the potential (Orszag 3/2 rule) |
| `FORCE_MMAX` | `0.` | $\log_{10}(M_{\rm max}/M_\odot)$ cap on the mass integrals; 0 means off |
| `MAX_EPSILON_NG` | `0.` | Divergence guard for the saddlepoint scheme |
| `WRITE_CGF_DIAG` | `False` | Per-evaluation diagnostics to `files/` (slow; debugging only) |

---

## Installation

Same as upstream `21cmFirstCLASS` — see
[notebook 1](Tutorial/notebook_1.ipynb) — with one change: **double-precision
FFTW is required** (`libfftw3` and `libfftw3_omp`, not the `f` variants).

```bash
# Debian / Ubuntu
sudo apt-get install libfftw3-dev libgsl-dev

git clone https://github.com/slibanore/21cmFirstCLASS_NG.git
cd 21cmFirstCLASS_NG
pip install -e .
```

You also need `classy`, the CLASS Python wrapper, as for upstream.

---

## Getting started

```python
import py21cmfast as p21c

lightcone = p21c.run_lightcone(
    redshift=6.,
    random_seed=1,                    # fix this when comparing runs
    user_params={
        "BOX_LEN": 300., "HII_DIM": 100,
        "NON_GAUSS_IC": True,
        "NON_GAUSS_FCOLL_COND": True,
        "NON_GAUSS_FCOLL_UNCOND": True,
        "EXTRA_DIM_FNL": 1.5,
    },
    cosmo_params={"F_NL": 200., "KCUT_FNL": 1e-3},
)
```

Then work through the tutorials:

* **[`Tutorial fNL/notebook_fNL_collapse_fraction.ipynb`](Tutorial%20fNL/notebook_fNL_collapse_fraction.ipynb)**
  — how $f_{\rm NL}$ changes the halo mass function. 
* **[`Tutorial fNL/notebook_fNL_simulations.ipynb`](Tutorial%20fNL/notebook_fNL_simulations.ipynb)**
  — running simulations and comparing against the Gaussian baseline, via
  `Tutorial fNL/runsims_FNL.py`.

The four upstream notebooks in [`Tutorial/`](Tutorial) cover installation and
the non-$f_{\rm NL}$ features of `21cmFirstCLASS`; notebooks 1 and 4 carry small
modifications here.

---

## Things to be careful about

**The expansions are asymptotic, not convergent.** Both the Edgeworth and the
saddlepoint schemes break down at large $|f_{\rm NL}|$ and large $\nu$. The code
warns when it detects this and falls back to the Gaussian mass function. 

**Watch the clipping.** For $f_{\rm NL} < 0$ the Edgeworth correction can pass
through zero at the high-mass end, where it is clipped at zero. Clipping at zero
rather than at one keeps the correction continuous in mass, which matters
because it is integrated over $\mathrm{d}\ln M$ by an adaptive routine. 

**The two barriers differ.** `dNion_General` uses the Sheth-Tormen barrier
$\sqrt{a}\,\delta_c$, matching the Gaussian mass function it corrects;
`dNdM_conditional` uses the plain spherical-collapse $\delta_c$. Each is
consistent with its own Gaussian baseline.

**The mass function must have a barrier.** The correction is derived in a
Press-Schechter framework. Applying it to the Watson fits (`HMF = 2, 3`), which
are calibrated to simulations and have no barrier, is uncontrolled. The code
warns; prefer `HMF = 1`.

**Memory.** With `NON_GAUSS_IC = True` the de-aliasing grid is
`EXTRA_DIM_FNL × DIM` on a side and memory scales as its cube. 
Building the three-point tables needs further memory.

**`WRITE_CGF_DIAG` needs `files/` to exist** in the working directory, or the
records are discarded with a warning. It writes one record per integrand
evaluation, so use it on small boxes only.


---

## Acknowledging

If you use this code, please cite the accompanying paper (Finish, Libanore
et al., arXiv:2610.xxxxx) together with the `21cmFirstCLASS` papers:

* Jordan Flitter and Ely D. Kovetz, *"New tool for 21-cm cosmology. I. Probing ΛCDM and beyond"*, Phys. Rev. D **109** (2024) 043512 ([arXiv:2309.03942](https://arxiv.org/abs/2309.03942)).
* Jordan Flitter and Ely D. Kovetz, *"New tool for 21-cm cosmology. II. Investigating the effect of early linear fluctuations"*, Phys. Rev. D **109** (2024) 043513 ([arXiv:2309.03948](https://arxiv.org/abs/2309.03948)).

and the `21cmFAST` papers:

* Mesinger, Furlanetto and Cen, MNRAS **411** (2011) 955 ([arXiv:1003.3878](https://arxiv.org/abs/1003.3878)).
* Muñoz, Qin, Mesinger, Murray, Greig and Mason, MNRAS **511** (2022) 3657 ([arXiv:2110.13919](https://arxiv.org/abs/2110.13919)).

The non-Gaussian formalism implemented here is from:

* Sabti, Muñoz and Blas, *"New Roads to the Small-scale Universe: Measurements of the Clustering of Matter with the High-redshift UV Galaxy Luminosity Function"* ([arXiv:2009.01245](https://arxiv.org/abs/2009.01245)) — the Edgeworth expansion and the $\kappa_3$ calculation.
* LoVerde, Miller, Shandera and Verde ([arXiv:0711.4126](https://arxiv.org/abs/0711.4126)) — the original Edgeworth derivation.
* D'Aloisio, Zhang, Jeong and Shapiro ([arXiv:1206.3305](https://arxiv.org/abs/1206.3305)) and Lidz, Baxter, Adshead and Dodelson ([arXiv:1304.8049](https://arxiv.org/abs/1304.8049)) — the conditional mass function.

`21cmFirstCLASS` integrates several other open source codes; please also cite
the relevant ones listed in the
[upstream README](https://github.com/jordanflitter/21cmFirstCLASS#acknowledging)
— `CLASS`, `HYREC-2`, `powerbox`, `21cmSense`, `AxiCLASS` and `dmeff-CLASS`.

## License

As upstream — see [LICENSE](LICENSE).
