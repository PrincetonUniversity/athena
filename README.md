athena
======
<!-- Jenkins Status Badge in Markdown (with view), unprotected, flat style -->
<!-- In general, need to be on Princeton VPN, logged into Princeton CAS, with ViewStatus access to Jenkins instance to click on unprotected Build Status Badge, but server is configured to whitelist GitHub -->
<!-- [![Jenkins Build Status](https://jenkins.princeton.edu/buildStatus/icon?job=athena/PrincetonUniversity_athena_jenkins_master)](https://jenkins.princeton.edu/job/athena/job/PrincetonUniversity_athena_jenkins_master/) -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11660592.svg)](https://doi.org/10.5281/zenodo.11660592)
[![codecov](https://codecov.io/gh/PrincetonUniversity/athena/branch/master/graph/badge.svg?token=ZzniY084kP)](https://codecov.io/gh/PrincetonUniversity/athena)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.0-4baaaa.svg)](code_of_conduct.md)

<!--[![Public GitHub  issues](https://img.shields.io/github/issues/PrincetonUniversity/athena-public-version.svg)](https://github.com/PrincetonUniversity/athena-public-version/issues)
[![Public GitHub pull requests](https://img.shields.io/github/issues-pr/PrincetonUniversity/athena-public-version.svg)](https://github.com/PrincetonUniversity/athena-public-version/pulls) -->

<p align="center">
	  <img width="345" height="345" src="https://user-images.githubusercontent.com/1410981/115276281-759d8580-a108-11eb-9fc9-833480b97f95.png">
</p>

Athena++ radiation GRMHD code and adaptive mesh refinement (AMR) framework

Please read [our contributing guidelines](./CONTRIBUTING.md) for details on how to participate.

## Citation
To cite Athena++ in your publication, please use the following BibTeX to refer to the code's [method paper](https://ui.adsabs.harvard.edu/abs/2020ApJS..249....4S/abstract):
```
@article{Stone2020,
	doi = {10.3847/1538-4365/ab929b},
	url = {https://doi.org/10.3847%2F1538-4365%2Fab929b},
	year = 2020,
	month = jun,
	publisher = {American Astronomical Society},
	volume = {249},
	number = {1},
	pages = {4},
	author = {James M. Stone and Kengo Tomida and Christopher J. White and Kyle G. Felker},
	title = {The Athena$\mathplus$$\mathplus$ Adaptive Mesh Refinement Framework: Design and Magnetohydrodynamic Solvers},
	journal = {The Astrophysical Journal Supplement Series},
}
```
Additionally, you can add a reference to `https://github.com/PrincetonUniversity/athena` in a footnote.

Finally, we have minted DOIs for each released version of Athena++ on Zenodo. This practice encourages computational reproducibility, since you can specify exactly which version of the code was used to produce the results in your publication. `10.5281/zenodo.4455879` is the DOI which cites _all_ versions of the code; it will always resolve to the latest release. Click on the Zenodo badge above to get access to BibTeX, etc. info related to these DOIs, e.g.:

```
@software{athena,
  author       = {Athena++ development team},
  title        = {PrincetonUniversity/athena: Athena++ v24.0},
  month        = jun,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {24.0},
  doi          = {10.5281/zenodo.11660592},
  url          = {https://doi.org/10.5281/zenodo.11660592}
}
```
