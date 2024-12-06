## Reference
* `euclid_obseq.pkl` contains the observation sequence for Euclid in pandas DataFrame format. It makes use of MultiIndex system to access the 5 kinds of observation for each pointing. See (Scaramella et al.)[https://arxiv.org/pdf/2108.01201] Fig. 8 for an illustration of the observation sequence (the timings used here are slightly different). 
It looks like:
```python
                         date  exptime         ra        dec          pa        lon        lat  pa_ecliptic        saa  filter   date_euclid  obs_id  pointing_id  patch_id  dither_id
visit obs_kind
0     VIS_LONG   61346.858513    566.0  65.025888 -36.679767  249.647989  51.745415 -56.852863   267.537147  89.376641     VIS  59155.858513       1            1         0          0
      VIS_SHORT  61346.862077     95.0  65.025888 -36.679767  249.647989  51.745415 -56.852863   267.537147  89.376641     VIS  59155.862077       1            1         0          0
      NISP_J     61346.865526     87.2  65.025888 -36.679767  249.647989  51.745415 -56.852863   267.537147  89.376641  NISP_J  59155.865526       1            1         0          0
      NISP_H     61346.866068     87.2  65.025888 -36.679767  249.647989  51.745415 -56.852863   267.537147  89.376641  NISP_H  59155.866068       1            1         0          0
      NISP_Y     61346.866610     87.2  65.025888 -36.679767  249.647989  51.745415 -56.852863   267.537147  89.376641  NISP_Y  59155.866610       1            1         0          0
1     VIS_LONG   61346.869173    566.0  65.032330 -36.644968  249.644155  51.773928 -56.821327   267.513288  89.384369     VIS  59155.869173       1            2         0          1
      VIS_SHORT  61346.872737     95.0  65.032330 -36.644968  249.644155  51.773928 -56.821327   267.513288  89.384369     VIS  59155.872737       1            2         0          1
      NISP_J     61346.876186     87.2  65.032330 -36.644968  249.644155  51.773928 -56.821327   267.513288  89.384369  NISP_J  59155.876186       1            2         0          1
      NISP_H     61346.876728     87.2  65.032330 -36.644968  249.644155  51.773928 -56.821327   267.513288  89.384369  NISP_H  59155.876728       1            2         0          1
      NISP_Y     61346.877270     87.2  65.032330 -36.644968  249.644155  51.773928 -56.821327   267.513288  89.384369  NISP_Y  59155.877270       1            2         0          1
...
```
To access the information of the long VIS observation for the pointing 42 you will do:
```python
import pandas as pd

obs_seq = pd.read_pickle("euclid_obseq.pkl")
info = obs_seq.loc[42, "VIS_LONG"]
ra, dec = info["ra"], info["dec"]
```
Note: The `date` parameter correspond to the Euclid planned observation shifted by 2 years to overlap with the observation sequence of the Roman/Rubin simulation. The original Euclid date is stored in `date_euclid`
Note 2: This observation sequence has been made from the `rsd2024a`: https://euclid.roe.ac.uk/dmsf/files/20839/view
Note 3: The observation timing are made from Figure 4.4 of MOCDC_v4.2: https://euclid.roe.ac.uk/dmsf/files/20821/view
