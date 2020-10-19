# pbchar
Apertif primary beam characterization

To do PB-correction, source-finding and cross-matching for all data in DR1:
`python run_beams.py --PB 191002 --ncores 24`
Runs on happili-05; only drift scans 190912 and 191002 currently implemented

To collect cross-matched, PB-corrected source information:
`python collect_beams.py /tank/adams/pbchar/191002`

To make first, simple plots (in a python interpreter):

```
import pbchar as pbchar
B00 = pbchar.PB('files/191002_00.csv')
B00.go()
```


On the do-do list:
- Write a wrapper function for PB characterization to run all compound beams for a given set of PB images
- Make sure run_beams (and beam.Beam) can handle PB images from Gaussian Process method in addition to drift scans
- Rewrite beam.Beam to use astropy/python exclusively, as opposed to miriad (mainly because of BPA=0 header issues, also will hopefully be more efficient)
