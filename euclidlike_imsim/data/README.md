## Reference
* `euclid_obseq.pkl` contains the observation sequence for Euclid in pandas DataFrame format. It makes use of MultiIndex system to access the 5 kinds of observation for each pointing. See (Scaramella et al.)[https://arxiv.org/pdf/2108.01201] Fig. 8 for an illustration of the observation sequence (the timing used here a slightly different). 
It looks like:
```python
                         date  exptime      ra     dec      pa     saa  filter   date_euclid  obs_id  pointing_id  patch_id  dither_id  line_obseq
visit obs_kind
0     VIS_LONG   61085.020984    566.0  51.745 -56.853 -92.463  89.377     VIS  60354.020984       1            1         1          0          77
      VIS_SHORT  61085.024549     95.0  51.745 -56.853 -92.463  89.377     VIS  60354.024549       1            1         1          0          77
      NISP_J     61085.027998     87.2  51.745 -56.853 -92.463  89.377  NISP_J  60354.027998       1            1         1          0          77
      NISP_H     61085.028539     87.2  51.745 -56.853 -92.463  89.377  NISP_H  60354.028539       1            1         1          0          77
      NISP_Y     61085.029081     87.2  51.745 -56.853 -92.463  89.377  NISP_Y  60354.029081       1            1         1          0          77
1     VIS_LONG   61085.033322    566.0  51.774 -56.821 -92.487  89.384     VIS  60354.033322       1            2         1          1          79
      VIS_SHORT  61085.036887     95.0  51.774 -56.821 -92.487  89.384     VIS  60354.036887       1            2         1          1          79
      NISP_J     61085.040336     87.2  51.774 -56.821 -92.487  89.384  NISP_J  60354.040336       1            2         1          1          79
      NISP_H     61085.040877     87.2  51.774 -56.821 -92.487  89.384  NISP_H  60354.040877       1            2         1          1          79
      NISP_Y     61085.041419     87.2  51.774 -56.821 -92.487  89.384  NISP_Y  60354.041419       1            2         1          1          79
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
