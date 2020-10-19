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
