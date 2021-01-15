# Cool dwarfs in multiple systems
## TBD
  
> This repository contains the code written for the work <a href="#" target="_blank">Cifuentes et al. 2021</a>.

[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![Publication](https://img.shields.io/badge/Published%3F-yes-brightgreen.svg)](https://www.aanda.org/articles/aa/abs/2020/10/aa38295-20/aa38295-20.html)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Ask Me Anything !](https://img.shields.io/badge/Ask%20me-anything-1abc9c.svg)](https://GitHub.com/ccifuentesr)

## Table of Contents

- [Installation](#installation)
- [Structure](#structure)
- [Support](#support)
- [License](#license)
- [Suggested Resources](#resources)
- [The CARMENES input catalogue of M dwarfs series](#volumes)

---

## Installation

> The files are self-contained, self-consistent, homogenoeusly formatted, fairly self-explanatory.

- The code is provided as `*.py` files meant to be run individually.
- They may be run as Python Notebooks. The symbol `# %%` starts a cell that can be run separately.
- Most of the `*.py` files require of additional data contained in the folders stored in the repository.
- Cloning or downloading the complete repository is strongly recommended (see below).
- The installation of some basic libraries is a prerequisite: `numpy`, `scipy`, `astropy`, `matplotlib` and `pyperclip`. Other modules are included in the Python distribution and do not need additional installation (e.g., `csv`).

### Clone

- Clone this repo to your local machine using `git clone https://github.com/ccifuentesr/cif21-multiplicity`, or
- Download this repo as a .zip and run the scripts in your local machine.

## Structure

### Directories

- Directory ./: Includes all the code files classified as detailed below and the master table (`Name.csv`). [XX files, YY MB]
- Directory ./Nnn: Includes ... [XX files, YY MB]

The total size is ZZ MB.

### Files

All files are named `cif20.xxx_yyy_zzz.py`, where `xxx` defines the kind of output that it produces, `yyy` gives additional information about the output, and `zzz` enumerates the main variables involved. For example, the script `cif20.plot_literature_Mabs_SpT.py` produces an absolute magnitude vs. spectral type plot, and compares the values with those of the literature. The complete list of files and their description can be found below sorted by output, in alphabetical order.

| File | Description | Input<sup id="a1">[1](#f1)</sup>| Output<sup id="a2">[2](#f2)</sup> | 
| --- | --- | --- | --- | 
| Xxx.py | ... | ... | ... |

1. <small id="f1"> ... </small> [↩](#a1) 
2. <small id="f2"> ... </small> [↩](#a2)

**Notes:** 
- If necessary.

### The master table

`Xxx.csv` is the full version of Table A.3 (summary table) in <a href="https://arxiv.org/abs/2007.15077" target="_blank">Cifuentes et al. 2020</a>, available at the CDS via anonymous ftp to cdsarc.u-strasbg.fr(130.79.128.5) or via [*link to CDS in progress*]. The complete version available here contains 2483 rows and 175 columns. It is stored in the root directory and can be manipulated separately with tabular data managemente software such as <a href="http://www.star.bris.ac.uk/~mbt/topcat/" target="_blank">TOPCAT</a>.

**Version history**

|	Version	|	Date	| Comments |
|	---	|	---	| --- |
| 01 | February 2021 | Current version |

Row-by-row description of `Xxx.csv`.

|	ID	|	Name	|	Units	|	Description	|	Annotations
|	---	|	---	|	---	|	---	|	---
|		|	ID_System	|	-	|	System identifier	|	
|		|	ID_Star	|	-	|	Star identifier 	|	Non-resolved spectroscopic binaries are identified with a single ID.
|		|	Component	|	-	|	Component identifier from this work	|	Alphabetic designation of the components (A is the primary). The use of lowercases designates instances where later observations have revealed a closer companion, such as spectroscopic pairs (e.g. Aab). Componentes are named based on its apparent magnitude in G.
|		|	Name	|	-	|	Discovery or most common name 	|	
|		|	RA_J2015	|	hms	|	Right ascension (J2016.0 epoch)	|	In equinox J2000, from Gaia EDR3.
|		|	DE_J2015	|	dms	|	Declination (J2016.0 epoch)	|	In equinox J2000, from Gaia EDR3.
|		|	Koenigstuhl	|	-	|	Star identifier (JHHMMm+DDdAAA)	|	Königstuhl designation (KO N, se the works by J. A. Caballero).
|		|	Karmn	|	-	|	Star identifier (JHHMMm+DDdAAA)	|	Carmencita identifier (JHHMMm+DDdAAA, Cortés-Contreras et al. 2016).
|		|	SpT	|	-	|	Spectral type 	|	When estimated, lowercase is used (e.g. m3 V).
|		|	SpTnum	|	-	|	Spectral type in numerical format 	|	SpTnum = -2.0 for K5V ; -1.0 for K7V ; 0.0 for M0.0V ; 0.5 for M0.5V ; ... ; 10.0 for L0.0 ; 10.5 for L0.5 ; etc.
|		|	SpT_ref	|	-	|	Reference for the spectral type 	|	See references after this table.
|		|	Discoverer	|	-	|	Reference for the discoverer	|	See note (X)
|		|	WDS_name	|	-	|	WDS name (based on J2000 position)	|	See Component annotation.
|		|	WDS_disc	|	-	|	Discoverer Code (1 to 4 letters) and Number	|	Originally 3 letters represent the discoverer; an additional 'A' denotes an appendix, 'B' a second appendix, e.g. in the lists of F. Struve: STF1004 is the 1004th system in the main list, STFA 11 is the 11th system of the first appendix, STFB 12 is the 12th system of the second appendix, etc (from WDS).
|		|	WDS_comp	|	-	|	Component identifier from WDS	|	Components when more than 2
|		|	Teff	|	K	|	Effective temperature from VOSA	|	BT-Settl CIFIST grid of synthetic spectra uncertainties are 50 K for Teff (25 K for Teff <= 2400 K).
|		|	logg	|	-	|	Surface gravity from VOSA	|	BT-Settl CIFIST grid of synthetic spectra uncertainties are 0.5 dex for logg.
|		|	Lbol	|	solLum	|	Bolometric luminosity from VOSA 	|	BT-Settl CIFIST grid of synthetic spectra metallicity is fixed to solar.
|		|	Lberr	|	solLum	|	Bolometric luminosity error from VOSA 	|	BT-Settl CIFIST grid of synthetic spectra metallicity is fixed to solar.
|		|	Mass_Lbol	|	solMass	|	Stellar mass computed from Lbol	|	Using the Mass vs. Lbol relation from Pecaut et al. 2013 (hotter than F7 V) and Cifuentes et al. 2020 (cooler than K5 V).
|		|	Mass_MG	|	solMass	|	Stellar mass computed from MG	|	Using the Mass vs. MG relation from Pecaut et al. 2013 (hotter than F7 V) and Cifuentes et al. 2020 (cooler than K5 V).
|		|	Mass_A	|	solMass	|	Stellar mass adopted for the A component	|	
|		|	Mass_B	|	solMass	|	Stellar mass adopted for the B component	|	
|		|	Mass_C	|	solMass	|	Stellar mass adopted for the C component	|	
|		|	ra_A	|	deg	|	Right ascension (J2016.0 epoch) for A	|	
|		|	ra_B	|	deg	|	Right ascension (J2016.0 epoch) for B	|	
|		|	ra_C	|	deg	|	Right ascension (J2016.0 epoch) for C	|	
|		|	dec_A	|	deg	|	Declination (J2016.0 epoch) for A	|	
|		|	dec_B	|	deg	|	Declination (J2016.0 epoch) for A	|	
|		|	dec_C	|	deg	|	Declination (J2016.0 epoch) for A	|	
|		|	pmra_A	|	mas a-1	|	Proper motion in Right Ascension for A	|	
|		|	pmra_B	|	mas a-1	|	Proper motion in Right Ascension for B	|	
|		|	pmra_C	|	mas a-1	|	Proper motion in Right Ascension for C	|	
|		|	pmdec_A	|	mas a-1	|	Proper motion in Declination for A	|	
|		|	pmdec_B	|	mas a-1	|	Proper motion in Declination for B	|	
|		|	pmdec_C	|	mas a-1	|	Proper motion in Declination for C	|	
|		|	d_A	|	pc	|	Distance for A component	|	
|		|	d_B	|	pc	|	Distance for B component	|	
|		|	d_C	|	pc	|	Distance for C component	|	
|		|	rho_AB	|	arcsec	|	Angular separation between A and B components	|	
|		|	rho_AC	|	arcsec	|	Angular separation between A and C components	|	
|		|	theta_AB	|	deg	|	Position angle between A and B components	|	
|		|	theta_AC	|	deg	|	Position angle between A and C components	|	
|		|	muratio_AB	|	-	|	mu ratio between A and B	|	As described by Montes et al. 2018.
|		|	muratio_AC	|	-	|	mu ratio between A and C	|	As described by Montes et al. 2018.
|		|	deltaPA_AB	|	deg	|	Difference of positional angle between A and B	|	As described by Montes et al. 2018.
|		|	deltaPA_AC	|	deg	|	Difference of positional angle between A and C	|	As described by Montes et al. 2018.
|		|	deltad_AB	|	-	|	Distance ratio between A and B	|	Described as ∆d/d_A.
|		|	deltad_AC	|	-	|	Distance ratio between A and C	|	Described as ∆d/d_A.
|		|	s_AB	|	au	|	Projected physical separation between A and B	|	
|		|	s_AC	|	au	|	Projected physical separation between A and C	|	
|		|	Ug_AB	|	1E33 J	|	Binding energy between A and B	|	
|		|	Ug_AC	|	1E33 J	|	Binding energy between A and C	|	
|		|	Porb_AB	|	a	|	Orbital period of A and B	|	
|		|	Porb_AC	|	a	|	Orbital period of A and C	|	
|		|	gaia_id	|	-	|	Gaia identification number	|	
|		|	ra	|	deg	|	Right ascension (J2016.0 epoch)	|	
|		|	ra_error	|	deg	|	Right ascension error (J2016.0 epoch)	|	
|		|	dec	|	deg	|	Declination (J2016.0 epoch)	|	
|		|	dec_error	|	deg	|	Declination error (J2016.0 epoch)	|	
|		|	parallax	|	mas	|	Parallax 	|	
|		|	parallax_error	|	mas	|	Parallax error 	|	
|		|	d_pc	|	pac	|	Distance 	|	
|		|	ed_pc	|	pac	|	Distance error 	|	
|		|	parallax_ref	|	-	|	Reference for the parallax 	|	
|		|	pm	|	mas a-1	|	Total proper motion 	|	
|		|	pmra	|	mas a-1	|	Proper motion in Right Ascension	|	
|		|	pmra_error	|	mas a-1	|	Proper motion error in Right Ascension	|	
|		|	pmdec	|	mas a-1	|	Proper motion in Declination	|	
|		|	pmdec_error	|	mas a-1	|	Proper motion error in Declination	|	
|		|	pm_ref	|	mas a-1	|	Reference for the proper motion	|	
|		|	ruwe	|	-	|	Re-normalised Unit Weight Error	|	See https://www.cosmos.esa.int/web/gaia/dr2-known-issues
|		|	phot_g_mean_mag	|	mag	|	*Gaia* EDR3 *G* magnitude	|	
|		|	phot_g_mean_mag_error	|	mag	|	*Gaia* EDR3 *G* magnitude error	|	
|		|	phot_bp_mean_mag	|	mag	|	*Gaia* EDR3 *BP* magnitude	|	
|		|	phot_bp_mean_mag_error	|	mag	|	*Gaia* EDR3 *BP* magnitude error	|	
|		|	phot_rp_mean_mag	|	mag	|	*Gaia* EDR3 *RP* magnitude	|	
|		|	phot_rp_mean_mag_error	|	mag	|	*Gaia* EDR3 *RP* magnitude error	|	
|		|	phot_bp_rp_excess_factor	|	mag	|	*BP*-*RP* excess factor	|	
|		|	radial_velocity	|	km s-1	|	Radial velocity	|	
|		|	radial_velocity_error	|	km s-1	|	Radial velocity error	|	
|		|	radial_velocity_ref	|	-	|	Radial velocity reference	|	
|		|	2MASS_id	|	-	|	2MASS identification number	|	
|		|	RAJ2000	|	deg	|	Right ascension (J2000 epoch)	|	
|		|	DEJ2000	|	deg	|	Declination (J2000 epoch)	|	
|		|	Jmag	|	mag	|	2MASS *J* magnitude	|	
|		|	e_Jmag	|	mag	|	2MASS *J* magnitude error	|	
|		|	Hmag	|	mag	|	2MASS *H* magnitude	|	
|		|	e_Hmag	|	mag	|	2MASS *H* magnitude error	|	
|		|	Kmag	|	mag	|	2MASS *Ks* magnitude	|	
|		|	e_Kmag	|	mag	|	2MASS *Ks* magnitude error	|	
|		|	Qfl	|	-	|	*JHKs* Photometric quality flag [ABCUXZ]	|	Three character flag, one character per band, that provides a summary of the net quality of the default photometry in each band (from 2MASS).
|		|	AllWISE	|	-	|	WISE identificator	|	WISE All-Sky Release Catalog name, based on J2000 position
|		|	W1mag	|	mag	|	WISE *W1* magnitude	|	
|		|	e_W1mag	|	mag	|	WISE *W1* magnitude error	|	
|		|	W2mag	|	mag	|	WISE *W2* magnitude	|	
|		|	e_W2mag	|	mag	|	WISE *W2* magnitude error	|	
|		|	W3mag	|	mag	|	WISE *W3* magnitude	|	
|		|	e_W3mag	|	mag	|	WISE *W3* magnitude error	|	
|		|	W4mag	|	mag	|	WISE *W4* magnitude	|	
|		|	e_W4mag	|	mag	|	WISE *W4* magnitude error	|	
|		|	qph	|	-	|	*W1W2W3W4* Photometric quality flag [ABCUXZ]	|	Four character flag, one character per band [W1/W2/W3/W4], that provides a shorthand summary of the quality of the profile-fit photometry measurement in each band, as derived from the measurement signal-to-noise ratio.

**Notes:** 
- Uncertainties in *Gaia* DR2 photometry are computed from the associated bolometric flux and its error.
- In the description, 'B' (not *B*) should be read 'the closest component in multiple systems'.

---

## Support

Reach out to me at <a href="mailto:ccifuentes@cab.inta-csic.es">`ccifuentes@cab.inta-csic.es`</a>.

---

## License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)

- **[MIT license](http://opensource.org/licenses/mit-license.php)**

---

## Suggested Resources

- <a href="https://www.python.org/dev/peps/pep-0008/" target="_blank">Style Guide for Python Code (PEP 8)</a>
- <a href="https://carmenes.caha.es" target="_blank">CARMENES Website</a>
