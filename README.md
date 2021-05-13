# One<sup>[*](#one)</sup> is the loneliest number: multiplicity in cool dwarfs
## C. Cifuentes<sup>1</sup>, J. A. Caballero<sup>1</sup> and S. Agustí<sup>2</sup>

1. <small id="f1">Centro de Astrobiología (CSIC-INTA), Madrid, Spain</small>
2. <small id="f2">Lyceè Français, Madrid, Spain</small>

<a href="https://zenodo.org/badge/latestdoi/329311083"><img src="https://zenodo.org/badge/329311083.svg" alt="DOI"></a>
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![Publication](https://img.shields.io/badge/Published%3F-no-orange.svg)]()
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Ask Me Anything !](https://img.shields.io/badge/Ask%20me-anything-1abc9c.svg)](https://GitHub.com/ccifuentesr)
[![Music](https://img.shields.io/badge/Music%3F-yes-green.svg)](https://www.youtube.com/watch?v=DYzY7-V5vxY)

> This repository contains the data, the code and complementary science for the publication <a href="#" target="_blank">Cifuentes et al. 2021</a> (submitted to RNAAS).

---

## Table of Contents

- [Summary](#summary)
- [Parameters](#parameters)
- [Structure](#structure)
- [Charts](#charts)
- [Results](#results)
- [References](#references)
- [Code usage](#usage)
- [Support](#support)
- [Suggested Resources](#resources)
- [One](#one)

---

## Summary

> Stars in multiple systems offer a unique opportunity to learn about stellar formation and evolution. As they settle down into stable configurations, multiple systems arrange in a variety of hierarchies and separations between the components. We examine 11 known and 11 newly discovered multiple systems including at least one M dwarf with the latest astrometric data from Gaia Early Data Release 3 (EDR3). We find that the individual components of systems at very wide separations are often multiple systems themselves.

Multiple star systems are particularly important in astrophysics.
Measuring their orbits reveals very relevant information about them, only accessible theoretically otherwise.
For instance, the masses of the components of these systems can be determined directly, applying fundamental laws of gravity.
Then, radii and densities can be estimated from the masses, and mass-luminosity relations can be drawn from further observation.

Stars that are close together in the sky can be revealed as visual or as physical binaries.
By analysing their distances, the relative movement and even their composition, we can determine if the pair is physically linked.
In this case the two objects, bounded by gravity, follow elliptical orbits around their common barycenter, and we refer to it as a binary system.
Very often, especially when the stars are distant, they appear as a single point of light to the unaided eye, but then are revealed as multiple by other means.
Indirect techniques, such as spectroscopy, astrometry or photometry, can reveal systems of this kind.
For example, the so called spectroscopic binaries are visible under a detailed study of signatures in their spectra, via the Doppler effect.

Systems of two or more stars are called multiple star systems.
When dinamically stable, they are organised hierarchically and consist of nested orbits that can be treated individually as a Two-body problems. 
For instance, in physical triple star systems, two stars form a close binary system, and the third orbits this pair at a distance much larger than that of the binary orbit.
Most multiple star systems are found to be triple.
On the other hand, trapezia systems are usually very young and unstable.
They interact strongly and chaotically in an N-body problem manner. 
The competition for stable orbits is resolved with the fragmentation into multiple hierarchical systems.

---

## Parameters

The following parameters are computed in this work. The subindex <img src="https://render.githubusercontent.com/render/math?math=B"> designates any other component different than the primary (i.e. B, C, D).

| Parameter | Unit | Formula | Annotations | 
| --- | --- | --- | --- | 
| Distance  | pc | <img src="https://render.githubusercontent.com/render/math?math=d = 1000/\varpi"> | <img src="https://render.githubusercontent.com/render/math?math=\varpi"> is the trigonometric parallax in milliarcseconds. |
| Physical separation  | au | <img src="https://render.githubusercontent.com/render/math?math=s = \rho d"> | |
| Apparent separation (e.g. A and B)  | arcsec | <img src="https://render.githubusercontent.com/render/math?math=\rho = 3600 \sqrt{(\alpha_A-\alpha_B)^2 %2B (\delta_A-\delta_B)^2}"> | |
| Total proper motion | km s-1 | <img src="https://render.githubusercontent.com/render/math?math=\mu = \sqrt{\mu_\alpha \cos{\delta}^2 %2B \mu_\delta^2}"> | |
| Orbital period | a |  <img src="https://render.githubusercontent.com/render/math?math=P_{\orb} = 2\pi \sqrt{\frac{a^3}{\mu}}"> | <img src="https://render.githubusercontent.com/render/math?math=a \sim s"> and <img src="https://render.githubusercontent.com/render/math?math=\mu = GM"> is the standard gravitational parameter, where <img src="https://render.githubusercontent.com/render/math?math=M"> is the mass of the more massive body.  |
| Binding energy | J |  <img src="https://render.githubusercontent.com/render/math?math=U_g = -M_A M_B/r"> | We approximate <img src="https://render.githubusercontent.com/render/math?math=r \sim s">, i.e. <img src="https://render.githubusercontent.com/render/math?math=-U_g^* = M_A M_B/s"> |

### Criteria for physical parity

The following parameters are used to discriminate between physical and visual pairs. 

| Parameter | Unit | Formula | Annotations | 
| --- | --- | --- | --- | 
| <img src="https://render.githubusercontent.com/render/math?math=\mu"> ratio | - |  <img src="https://render.githubusercontent.com/render/math?math=(\mu {\ratio})^2 = \frac{(\mu_\alpha \cos{\delta}_A - \mu_\alpha \cos{\delta}_B)^2 %2B (\mu_\delta_A - \mu_\delta_B)^2}{(\mu_\alpha \cos{\delta})^2 %2B (\mu_\delta_A)^2} "> | |
| Proper motion position angle difference (PA) | deg |  <img src="https://render.githubusercontent.com/render/math?math=\Delta PA = \lvert PA_A-PA_B \rvert"> | <img src="https://render.githubusercontent.com/render/math?math=PA_i"> is the angle between <img src="https://render.githubusercontent.com/render/math?math=\mu_\alpha \cos{\delta}_i"> and <img src="https://render.githubusercontent.com/render/math?math=\mu_{\delta,i}">.|
| <img src="https://render.githubusercontent.com/render/math?math=\Delta d"> | - | <img src="https://render.githubusercontent.com/render/math?math=\Delta d = \lvert d_A-d_B \rvert/d_A"> | |

**Notes:** 
- <img src="https://render.githubusercontent.com/render/math?math=\mu"> ratio and <img src="https://render.githubusercontent.com/render/math?math=\Delta PA"> defined by <a href="https://ui.adsabs.harvard.edu/abs/2018MNRAS.479.1332M/abstract" target="_blank">Montes et al. 2018</a> (Mon18).
- We define positive physical parity if:
	- <img src="https://render.githubusercontent.com/render/math?math=\mu"> ratio < 0.15, and
	- <img src="https://render.githubusercontent.com/render/math?math=\Delta PA"> < 15 deg, and
	- <img src="https://render.githubusercontent.com/render/math?math=\Delta d"> < 0.10.

### Bolometric luminosities

All luminosities were computed using the Virtual Observatory SED Analyzer (<a href="http://svo2.cab.inta-csic.es/theory/vosa/" target="_blank">VOSA</a>), using the BT-Settl CIFIST grid of models ([Fe/H] = 0.0). Close binaries (<img src="https://render.githubusercontent.com/render/math?math=\rho"> < 5 arcsec) are dismissed in this calculation.

### Masses

By priority:

**Later or equal than K5 V:**
1. Via bolometric luminosity using Cif20.
2. Via spectral type using Cif20.
3. Via absolute magnitud in G using Cif20 

**Earlier than K5 V:**
1. Via bolometric luminosity using Pec13.
2. Via spectral type using Pec13.
3. Via absolute magnitud in G using Pec13.

### Spectral types

By priority:

**Later or equal than K5 V:**
1. Via bolometric luminosity using Cif20.
2. Via absolute magnitud in G using Cif20.
	
**Earlier than K5 V:**
1. Via bolometric luminosity using Pec13.
2. Via absolute magnitud in G using Pec13.
3. Via absolute magnitud in J using Pec13.
4. Via interpolation with known spectral type using absolute magnitud in J (only if 1-3 are not available).

---

## Charts

The file `cif21.charts.pdf` is a 22-page pdf document that displays the 22 systems investigated in this work.
Each page contains a single system, with the following information:

- WDS name and discover(s) code(s), if the system is known, or 'New', if the system is reported for the first time in this work.
- Common name of the resolved components.
- Finding chart with a visual representation of the proper motion of each resolved component.
- Description of the system by pairs, including a quantitative analysis of physical connection, based on the criteria described above. 
- Astro-photometric description of the resolved components.   
- Brief one-paragraph description of the system.

As a visual aid, values that comply with the criteria for physical binding (see section [Parameters](#parameters)) and good astrometric quality data (i.e. RUWE < 1.4), are displayed in green. In case of poor or not complying, they are displayed in orange or red, respectively.

<p float="center">
 <img src="https://github.com/ccifuentesr/cif21-multiplicity/blob/main/chart_example.png" width="50%" />
</p>

---

## Structure

### Directories

- Directory ./: Stores all the files detailed below, including the master table (`cif21.multiplicity.csv`). [19 files, 54 MB]

### Files

The complete list of files and their description goes as follows:

| File | Type | Description | 
| --- | --- | --- | 
| `cif21.multiplicity.csv` | Main table | The master table as described below |
| `Mamajek_Pec13.csv` | Input  | Tabular data from <a href="https://ui.adsabs.harvard.edu/abs/2013ApJS..208....9P/abstract" target="_blank">Pecaut & Mamajek (2013)</a> |
| `KO6AB.csv` | Input | Positions during 11 epochs (1953-2015) for the system KO6 AB |
| `cif21.charts.pdf` | Code | Obtains the separation as a function of time spanning several epochs (e.g. `KO6AB.csv`) |
| `cif21.charts.py` | Code | Preliminar LaTeX chart for describing multiple systems |
| `cif21.dphot.py` | Code |  Photometric distance formulas using _G_ and _J_ magnitudes (Table 5 in <a href="https://ui.adsabs.harvard.edu/abs/2020A%26A...642A.115C/abstract" target="_blank">Cifuentes et al. 2020</a>) |
| `cif21.MR.py` | Code | Radii and masses from Stefan-Boltzmann and <a href="https://ui.adsabs.harvard.edu/abs/2019A%26A...625A..68S/abstract" target="_blank">Schweitzer et al. (2019)</a> |
| `cif21.params.py` | Code | Computes parameters to decide on the multiplicity of a system (see section [Parameters](#parameters)) |
| `cif21.plots.py` | Code | Produces the figures shown in section [Results](#results) |
| `cif21.rho_epochs.py` | Code | Obtains the separation as a function of time spanning several epochs (e.g. `KO6AB.csv`, see section [Results](#results)) |
| `*.png` | Images | Image files included in section [Results](#results) |


### The master table

`cif21.multiplicity.csv` contains **53** rows and **108** columns. It is stored in the root directory and can be manipulated separately with tabular data management software such as <a href="http://www.star.bris.ac.uk/~mbt/topcat/" target="_blank">TOPCAT</a>.

**Row-by-row description of `cif21.multiplicity.csv`.**

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
|		|	Discoverer	|	-	|	Reference for the discoverer	|	See references after this table.
|		|	WDS_name	|	-	|	WDS name (based on J2000 position)	|	See Component annotation.
|		|	WDS_disc	|	-	|	Discoverer Code (1 to 4 letters) and Number	|	Originally 3 letters represent the discoverer; an additional 'A' denotes an appendix, 'B' a second appendix, e.g. in the lists of F. Struve: STF1004 is the 1004th system in the main list, STFA 11 is the 11th system of the first appendix, STFB 12 is the 12th system of the second appendix, etc (from WDS).
|		|	WDS_comp	|	-	|	Component identifier from WDS	|	Components when more than 2
|		|	Teff	|	K	|	Effective temperature from VOSA	|	BT-Settl CIFIST grid of synthetic spectra uncertainties are 50 K for Teff (25 K for Teff <= 2400 K).
|		|	logg	|	-	|	Surface gravity from VOSA	|	BT-Settl CIFIST grid of synthetic spectra uncertainties are 0.5 dex for logg.
|		|	Lbol	|	solLum	|	Bolometric luminosity from VOSA 	|	BT-Settl CIFIST grid of synthetic spectra metallicity is fixed to solar.
|		|	Lberr	|	solLum	|	Bolometric luminosity error from VOSA 	|	BT-Settl CIFIST grid of synthetic spectra metallicity is fixed to solar.
|		|	Mass_Lbol	|	solMass	|	Stellar mass computed from Lbol	|	Using the Mass vs. Lbol relation from <a href="https://ui.adsabs.harvard.edu/abs/2013ApJS..208....9P/abstract" target="_blank">Pecaut et al. 2013</a> (hotter than F7 V) and <a href="https://ui.adsabs.harvard.edu/abs/2020A%26A...642A.115C/abstract" target="_blank">Cifuentes et al. 2020</a> (cooler than K5 V).
|		|	Mass_MG	|	solMass	|	Stellar mass computed from MG	|	Using the Mass vs. MG relation from <a href="https://ui.adsabs.harvard.edu/abs/2013ApJS..208....9P/abstract" target="_blank">Pecaut et al. 2013</a> (hotter than F7 V) and <a href="https://ui.adsabs.harvard.edu/abs/2020A%26A...642A.115C/abstract" target="_blank">Cifuentes et al. 2020</a> (cooler than K5 V).
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
|		|	muratio_AB	|	-	|	mu ratio between A and B	|	As described by Mon18.
|		|	muratio_AC	|	-	|	mu ratio between A and C	|	As described by Mon18.
|		|	deltaPA_AB	|	deg	|	Difference of positional angle between A and B	|	As described by Mon18.
|		|	deltaPA_AC	|	deg	|	Difference of positional angle between A and C	|	As described by Mon18.
|		|	deltad_AB	|	-	|	Distance ratio between A and B	|	Described as <img src="https://render.githubusercontent.com/render/math?math=d_A-d_B/d_A"> .
|		|	deltad_AC	|	-	|	Distance ratio between A and C	|	Described as <img src="https://render.githubusercontent.com/render/math?math=\Delta d/d_A">.
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
- Uncertainties in *Gaia* EDR3 photometry are computed from the associated bolometric flux and its error.
- In the description, 'B' (not *B*) should be read 'the closest component in multiple systems'. This apparent misnomer turns out to be convenient to designate a complete system in a single row of the table, with independence of the number of components.

---

## Results

<p float="center">
  <img src="https://github.com/ccifuentesr/cif21-multiplicity/blob/main/cif21_deltaPA_muratio.png" width="49%" />
  <img src="https://github.com/ccifuentesr/cif21-multiplicity/blob/main/cif21_d2_d1.png" width="49%" /> 
Left: <img src="https://render.githubusercontent.com/render/math?math=\mu"> ratio against proper motion position angle difference. Right: Distances comparison between primary and additional components. Crosses denote sources that do not satisfy the criteria for physical parity as described in this document.
</p>
<p float="center">
  <img src="https://github.com/ccifuentesr/cif21-multiplicity/blob/main/cif21_Ug_M.png" width="49%" />
  <img src="https://github.com/ccifuentesr/cif21-multiplicity/blob/main/cif21_MG_GJ.png" width="49%" /> 
Left: Binding energy against total mass. Right: Absolute magnitude <img src="https://render.githubusercontent.com/render/math?math=M_G"> against <img src="https://render.githubusercontent.com/render/math?math=G-J"> colour diagram. HD 134494 is a star that we re-classify as sub-giant, PYC J07311+4556 is a young candidate to the AB Doradus group, and LSPM J1633+0311S is a white dwarf.
</p>

### The case of Königstuhl 6 AB

LP 209-28 and LP 209-27 (KO6 AB) is a pair proposed as a binary system by <a href="https://ui.adsabs.harvard.edu/abs/2012Obs...132..252C/abstract" target="_blank">Caballero et al. 2012</a>. *Gaia* EDR3 introduces a notable dissimilarity in distance between components, accompanied by a good single-star model fitting (i.e. RUWE < 1.4). This means that we do not expect additional multiplicity in either component, and the difference in parallactic determination would not be the cause of close, unresolved companions. Therefore, we propose this system as a visual pair.

<p float="center">
  <img src="https://github.com/ccifuentesr/cif21-multiplicity/blob/main/KO6AB.png" width="80%" />
</p>

Additionally, measuring the projected separation, <img src="https://render.githubusercontent.com/render/math?math=\rho">, as a function of time for 10 epochs of observation available spanning 62 years (from March 16, 1953 to August 5, 2015), we find that the separation between the components increases in 0.248 arcsec, which implies 0.004 arcsec per year.

---

## References

The following are the bibliographic references in alphabetical order for the data used in this work.

### Spectral type

<a href="#"> </a>

| Reference | Bibcode |
| --- | --- |
| AF15 | <a href="#">2015A&A...577A.128A</a> |
| Bar14 | <a href="#">2014ApJ...794..143B</a> |
| Cab10 | <a href="#">2010A&A...520A..91C</a> |
| Can93 | <a href="#">1993yCat.3135....0C</a> |
| Cab07b | <a href="#">2007ApJ...667..520C</a> |
| Cif20 | <a href="#">2020A&A...642A.115C</a> |
| Clo02 | <a href="#">2002ApJ...567L..53C</a> |
| Cru03 | <a href="#">2003AJ....126.2421C</a> |
| Dhi10 | <a href="#">2010AJ....139.2566D</a> |
| Gau12 | <a href="#">2012MNRAS.427.2457G</a> |
| Gra03 | <a href="#">2003AJ....126.2048G</a> |
| Gra06 | <a href="#">2006AJ....132..161G</a> |
| Hou88 | <a href="#">1988MSS...C04....0H</a> |
| Hou99 | <a href="#">1999MSS...C05....0H</a> |
| Joh86 | <a href="#">1986ApJ...310..354J</a> |
| Joy49 | <a href="#">1949ApJ...109..231J</a> |
| Law08 | <a href="#">2008MNRAS.384..150L</a> |
| Lep13 | <a href="#">2013AJ....145..102L</a> |
| Mon18 | <a href="#">2018MNRAS.479.1332M</a> |
| New14 | <a href="#">2014AJ....147...20N</a> |
| Pec13 | <a href="#">2013ApJS..208....9P</a> |
| Rei08 | <a href="#">2008AJ....136.1290R</a> |
| Ria06 | <a href="#">2006AJ....132..866R</a> |
| Ste86 | <a href="#">1986AJ.....92..139S</a> |

### Discoverer

| Reference | Bibcode |
| --- | --- |
| Cab07a | <a href="#">2007A&A...462L..61C</a> (KO 1) |
| Cab07b | <a href="#">2007ApJ...667..520C</a> (KO 2 & KO 3) |
| Cab12a | <a href="#">2012Obs...132....1C</a> (KO 4) |
| Cab12b | <a href="#">2012Obs...132..176C</a> (KO 5) |
| Cab12c | <a href="#">2012Obs...132..252C</a> (KO 6) |
| Dhi10 | <a href="#">2010AJ....139.2566D</a> |
| Dun29 | <a href="#">1829MmRAS&3..257D</a> |
| Gau12 | <a href="#">2012MNRAS.427.2457G</a> |
| Gic61 | <a href="#">1961LowOB...5...61G</a> |
| Her26 | <a href="#">1826MmRAS...2..459H</a> |
| Hor12 | <a href="#">2012AJ....144..165H</a> |
| Jan06 | <a href="#">2006A&A...453..609J</a> |
| Jan06 | <a href="#">2006A&A...453..609J</a> |
| Jan14 | <a href="#">2014ApJ...789..102J</a> |
| Kna15 | <a href="#">2015JDSO...11..384K</a> |
| Law08 | <a href="#">2008MNRAS.384..150L</a> |
| Lep01 | <a href="#">2001AJ....122.3407L</a> |
| Lep13 | <a href="#">2013AJ....145..102L</a> |
| Ric12 | <a href="#">2012JDSO....8..160R</a> |
| Tok79 | <a href="#">1979SvAL....5..229T</a> |
| Vys56 | <a href="#">1956AJ.....61..201V</a> |

### Parallaxes and proper motions

| Reference | Bibcode |
| --- | --- |
| Dit14 | <a href="#">2014ApJ...784..156D</a> |
| Gaia2 | <a href="#">2018yCat.1345....0G</a> |
| GaiaEDR3 | <a href="#">2021A&A...649A...1G</a> |
| Lep05 | <a href="#">2005AJ....129.1483L</a> |

### Radial velocity

| Reference | Bibcode |
| --- | --- |
| Bur15 | <a href="#">2015ApJS..220...18B</a> |
| Gaia2 | <a href="#">2018yCat.1345....0G</a> |
| New14 | <a href="#">2014AJ....147...20N</a> |
| Shk10 | <a href="#">2010ApJ...716.1522S</a> |
| Ter15 | <a href="#">2015ApJS..220...16T</a> |

---
     
## Code usage

> The files are self-contained, self-consistent, homogenoeusly formatted, fairly self-explanatory.

- The code is provided as `*.py` files meant to be run individually.
- They may be run as Python Notebooks. The symbol `# %%` starts a cell that can be run separately.
- Cloning or downloading the complete repository is strongly recommended (see below).
- The installation of some basic libraries is a prerequisite: `numpy`, `scipy`, `astropy`, `matplotlib` and `pyperclip`. Other modules are included in the Python distribution and do not need additional installation (e.g. `csv`).

### Clone

- Clone this repo to your local machine using `git clone https://github.com/ccifuentesr/cif21-multiplicity`, or
- Download this repo as a .zip and run the scripts in your local machine.

---

## Support

Reach out to me at <a href="mailto:ccifuentes@cab.inta-csic.es">`ccifuentes@cab.inta-csic.es`</a>.

---

## Suggested Resources

- <a href="https://www.python.org/dev/peps/pep-0008/" target="_blank">Style Guide for Python Code (PEP 8)</a>
- <a href="https://carmenes.caha.es" target="_blank">CARMENES Website</a>

---

## <sup>*</sup>One

_One is the loneliest number that you'll ever do  
Two can be as bad as one  
It’s the loneliest number since the number one_  

<a href="https://www.youtube.com/watch?v=DYzY7-V5vxY" target="_blank">One</a>, Harry Nilsson, _Aerial Ballet_ (1968)

> One is a song written and recorded by Harry Nilsson and made famous by Three Dog Night. It is known for its opening line "_One is the loneliest number that you'll ever do_". Nilsson wrote the song after calling someone and getting a busy signal. He stayed on the line listening to the "beep, beep, beep, beep..." tone, writing the song. The busy signal became the opening notes.
> Three Dog Night played <a href="https://open.spotify.com/track/0TGKBG5wK1ZGSACf6uso3H?si=Lsc8MKNcSOiEyZv3XA961Q" target="_blank">recording</a> in the key of F minor, and it was released as the second single from Three Dog Night's eponymous first album. It became their first of seven gold records over the next five years and reached number five on the U.S. Billboard Hot 100 in 1969 and number four in Canada. (From <a href="https://en.wikipedia.org/wiki/One_(Harry_Nilsson_song)" target="_blank">Wikipedia</a>).
