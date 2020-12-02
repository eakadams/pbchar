# pbchar
Apertif primary beam characterization

To do PB-correction, source-finding and cross-matching for all data in DR1:
`python run_beams.py --PB 191002 --ncores 24`
Runs on happili-05; only drift scans 190912 and 191002 currently implemented

To collect cross-matched, PB-corrected source information:
`python combine_beams.py /tank/adams/pbchar/191002`

To make plots: <br>
`python run_pbchar.py filebase` <br>
To see options for running, including filtering on matches, turning on plots, etc: <br>
`python run_pbchar.py` -h

To make first, simple plots (in a python interpreter):

```
import pbchar as pbchar
B00 = pbchar.PB('files/191002_00.csv')
B00.go()
```


On the do-do list:
- Make specification of drift scan to run_beams more robust (look for normal stored files, not resized)
- Compare peak to total fluxes (more for understanding flux issues)
- Do intra-Apertif comparisons using MDS fields
- Add a Bayesian / MCMC analysis to get intrinsic scatter in flux scale
- Rewrite beam.Beam to use astropy/python exclusively, as opposed to miriad (mainly because of BPA=0 header issues, also will hopefully be more efficient)
