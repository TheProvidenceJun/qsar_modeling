<div style="background-color:lightgreen; border-style: dotted;border-color: yellow;text-align:center;color:blue;font-size:20px;border-radius: 50px;padding: 20px 20px;">
<strong> pyQSAR   (2025.11) <br>
    <br>
 chemoinformatics Lab., Soongsil University, Seoul, Korea<br>
    </strong></div>

## Developers
Sinyoung Kim, Sunghyun Moon, Minhyeong Kim, and Kwang-Hwi Cho (chokh@ssu.ac.kr)<br>

## Revision History
|Year| Date | Description |
|---|:---|:---|
|2017| December 27 | Created |
|2020| April 25    | migrated to python3 |
|2021| November 16 | Add Genetic Algorithm |
|2022| Feburary 28 | Add Clustering Algorithm (K-Means, S.O.M) |
|2023| November 16 | Filtering for Mordred, SDF generation, data split (Random or Block) |


## Tested Python Package 
| Package Name | version | Description |
|:--|:--|:--|
|  python | 3.9.16| |
| bokeh |2.4.3| |
| ipywidgets |8.0.2| |
| IPython | 7.31.1 | |
| matplotlib |3.3.4| | 
| numpy  | 1.24.2 | |
|  pandas | 1.5.3 | |
| PIL | 8.2.0| |
| py3Dmol | 2.0.3||
| plotly| 5.6.0| |
| pubchempy | 1.0.4 | |
| rdkit | 2022.09.5 | 2022.02 Rdkit doesn't works with version 3.8 or higher|
| scipy |1.10.1 | | |
| seaborn |0.11.1 | |
|  sklearn | 1.2.2| | |



```python
#from pyqsar import data_setting as ds
from pyqsar import data_tools as dt
from pyqsar import model_tools as mt
from pyqsar import draw_mol
import pandas as pd
import numpy as np
```


```python
dt.VCHECK()
```

    pyqsar: 2023-11-20
    Numpy : 1.24.2
    Scipy : 1.10.1
    pandas : 1.5.3
    RDkit : 2022.09.5
    sklearn : 1.2.2
    bokeh : 2.4.3
    py3Dmol : 2.0.3
    IPython : 7.31.1
    ipywidgets : 8.0.2
    pubchempy : 1.0.4



```python
## Set the maximum number of rows and columns you want to see in a data frame
#pd.set_option('display.max_rows',None)
#pd.set_option('display.max_columns',None)
```


```python
dt.man('pyQsar.PNG')
```


    
![png](output_6_0.png)
    



```python
# OLD : needs modification
#dt.man('workflow.png')
```

<h1>1. Data processing</h1>
<h2>1.1. Raw Data Read</h2>
<ul>
<li> the sample data(qdb.csv) is taken from 
<a href="https://qsardb.org/repository/handle/10967/229">qsarDB</a> </li>
<li> Training set(194) and Test set(98) are combined to make qdb.csv(292).</li>
<li> Then, descriptors are calculated using mordred.</li>
</ul>


```python
orgFormat = pd.read_csv("qdb.csv") 
orgFormat
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Name</th>
      <th>Cid</th>
      <th>Smiles</th>
      <th>logLC50</th>
      <th>ABC</th>
      <th>ABCGG</th>
      <th>nAcid</th>
      <th>nBase</th>
      <th>SpAbs_A</th>
      <th>SpMax_A</th>
      <th>...</th>
      <th>SRW10</th>
      <th>TSRW10</th>
      <th>MW</th>
      <th>AMW</th>
      <th>WPath</th>
      <th>WPol</th>
      <th>Zagreb1</th>
      <th>Zagreb2</th>
      <th>mZagreb1</th>
      <th>mZagreb2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>m001</td>
      <td>11859</td>
      <td>C1=C(C=C(C(=C1Cl)Cl)Cl)O</td>
      <td>-5.46</td>
      <td>7.427747</td>
      <td>7.188678</td>
      <td>0</td>
      <td>0</td>
      <td>11.643052</td>
      <td>2.307250</td>
      <td>...</td>
      <td>9.078065</td>
      <td>39.748909</td>
      <td>195.924948</td>
      <td>15.071150</td>
      <td>110</td>
      <td>13</td>
      <td>48.0</td>
      <td>54.0</td>
      <td>4.944444</td>
      <td>2.222222</td>
    </tr>
    <tr>
      <th>1</th>
      <td>m002</td>
      <td>11899</td>
      <td>C1=CC(=C(C=C1Cl)Cl)[N+](=O)[O-]</td>
      <td>-4.66</td>
      <td>8.134854</td>
      <td>7.770338</td>
      <td>0</td>
      <td>0</td>
      <td>12.675204</td>
      <td>2.302776</td>
      <td>...</td>
      <td>9.094144</td>
      <td>41.023148</td>
      <td>190.954084</td>
      <td>13.639577</td>
      <td>150</td>
      <td>14</td>
      <td>52.0</td>
      <td>58.0</td>
      <td>5.194444</td>
      <td>2.472222</td>
    </tr>
    <tr>
      <th>2</th>
      <td>m003</td>
      <td>2723704</td>
      <td>CNC(=S)N</td>
      <td>-3.98</td>
      <td>3.047207</td>
      <td>3.305183</td>
      <td>0</td>
      <td>0</td>
      <td>5.226252</td>
      <td>1.847759</td>
      <td>...</td>
      <td>6.834109</td>
      <td>27.254130</td>
      <td>90.025169</td>
      <td>8.184106</td>
      <td>18</td>
      <td>2</td>
      <td>16.0</td>
      <td>14.0</td>
      <td>3.361111</td>
      <td>1.333333</td>
    </tr>
    <tr>
      <th>3</th>
      <td>m004</td>
      <td>3032338</td>
      <td>CCNC(=S)N</td>
      <td>-4.00</td>
      <td>3.754314</td>
      <td>4.057055</td>
      <td>0</td>
      <td>0</td>
      <td>6.155367</td>
      <td>1.902113</td>
      <td>...</td>
      <td>7.131699</td>
      <td>29.439488</td>
      <td>104.040819</td>
      <td>7.431487</td>
      <td>32</td>
      <td>3</td>
      <td>20.0</td>
      <td>18.0</td>
      <td>3.611111</td>
      <td>1.583333</td>
    </tr>
    <tr>
      <th>4</th>
      <td>m005</td>
      <td>2346</td>
      <td>C1=CC=C(C=C1)CN=C=S</td>
      <td>-6.54</td>
      <td>7.071068</td>
      <td>6.547760</td>
      <td>0</td>
      <td>0</td>
      <td>12.932143</td>
      <td>2.154341</td>
      <td>...</td>
      <td>8.438366</td>
      <td>38.130322</td>
      <td>149.029920</td>
      <td>8.766466</td>
      <td>133</td>
      <td>9</td>
      <td>42.0</td>
      <td>44.0</td>
      <td>3.111111</td>
      <td>2.500000</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>287</th>
      <td>m288</td>
      <td>15292</td>
      <td>CC1=CC(=CC(=C1O)Cl)Cl</td>
      <td>-5.65</td>
      <td>7.427747</td>
      <td>7.188678</td>
      <td>0</td>
      <td>0</td>
      <td>11.643052</td>
      <td>2.307250</td>
      <td>...</td>
      <td>9.078065</td>
      <td>39.748909</td>
      <td>175.979570</td>
      <td>10.998723</td>
      <td>110</td>
      <td>13</td>
      <td>48.0</td>
      <td>54.0</td>
      <td>4.944444</td>
      <td>2.222222</td>
    </tr>
    <tr>
      <th>288</th>
      <td>m289</td>
      <td>15767</td>
      <td>COC1=C(C(=C(C(=C1Cl)Cl)Cl)Cl)Cl</td>
      <td>-7.01</td>
      <td>9.496696</td>
      <td>9.522267</td>
      <td>0</td>
      <td>0</td>
      <td>15.688761</td>
      <td>2.426720</td>
      <td>...</td>
      <td>9.602855</td>
      <td>44.616134</td>
      <td>277.862653</td>
      <td>17.366416</td>
      <td>220</td>
      <td>23</td>
      <td>64.0</td>
      <td>77.0</td>
      <td>6.916667</td>
      <td>3.000000</td>
    </tr>
    <tr>
      <th>289</th>
      <td>m290</td>
      <td>1551919</td>
      <td>CCCCNC(=S)N</td>
      <td>-3.85</td>
      <td>5.168527</td>
      <td>5.361851</td>
      <td>0</td>
      <td>0</td>
      <td>8.762573</td>
      <td>1.949856</td>
      <td>...</td>
      <td>7.475906</td>
      <td>33.090360</td>
      <td>132.072119</td>
      <td>6.603606</td>
      <td>79</td>
      <td>5</td>
      <td>28.0</td>
      <td>26.0</td>
      <td>4.111111</td>
      <td>2.083333</td>
    </tr>
    <tr>
      <th>290</th>
      <td>m291</td>
      <td>2566</td>
      <td>CC1(CC2=C(O1)C(=CC=C2)OC(=O)NC)C</td>
      <td>-6.52</td>
      <td>12.367162</td>
      <td>11.107564</td>
      <td>0</td>
      <td>0</td>
      <td>19.625463</td>
      <td>2.470049</td>
      <td>...</td>
      <td>9.709599</td>
      <td>62.929068</td>
      <td>221.105193</td>
      <td>7.132426</td>
      <td>425</td>
      <td>22</td>
      <td>84.0</td>
      <td>97.0</td>
      <td>6.256944</td>
      <td>3.472222</td>
    </tr>
    <tr>
      <th>291</th>
      <td>m292</td>
      <td>13905</td>
      <td>CCNC1=NC(=NC(=N1)SC)NCC</td>
      <td>-3.63</td>
      <td>9.899495</td>
      <td>9.242468</td>
      <td>0</td>
      <td>0</td>
      <td>17.184731</td>
      <td>2.298128</td>
      <td>...</td>
      <td>9.110410</td>
      <td>44.568778</td>
      <td>213.104816</td>
      <td>7.348442</td>
      <td>318</td>
      <td>17</td>
      <td>62.0</td>
      <td>68.0</td>
      <td>5.333333</td>
      <td>3.500000</td>
    </tr>
  </tbody>
</table>
<p>292 rows × 1830 columns</p>
</div>



<h2>1.2. format change</h2>
<ul>
<li> pyQsar needs the data with the following orders </li>
<li> "ID",  "EP",  "SMI", ~Descriptors~~~~ </li>
<li> "ID" = index or name, "EP" = End Point, "SMI" = smiles code</li>
</ul>


```python
changed = orgFormat.rename(columns={'Smiles': 'SMI', 'Cid': 'ID', "logLC50" : "EP"})
newFormat = changed[["ID", "EP", "SMI"]+ list(changed.columns)[4:]]
newFormat.to_csv("qdb1.csv", index=False)
newFormat
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ID</th>
      <th>EP</th>
      <th>SMI</th>
      <th>ABC</th>
      <th>ABCGG</th>
      <th>nAcid</th>
      <th>nBase</th>
      <th>SpAbs_A</th>
      <th>SpMax_A</th>
      <th>SpDiam_A</th>
      <th>...</th>
      <th>SRW10</th>
      <th>TSRW10</th>
      <th>MW</th>
      <th>AMW</th>
      <th>WPath</th>
      <th>WPol</th>
      <th>Zagreb1</th>
      <th>Zagreb2</th>
      <th>mZagreb1</th>
      <th>mZagreb2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>11859</td>
      <td>-5.46</td>
      <td>C1=C(C=C(C(=C1Cl)Cl)Cl)O</td>
      <td>7.427747</td>
      <td>7.188678</td>
      <td>0</td>
      <td>0</td>
      <td>11.643052</td>
      <td>2.307250</td>
      <td>4.614501</td>
      <td>...</td>
      <td>9.078065</td>
      <td>39.748909</td>
      <td>195.924948</td>
      <td>15.071150</td>
      <td>110</td>
      <td>13</td>
      <td>48.0</td>
      <td>54.0</td>
      <td>4.944444</td>
      <td>2.222222</td>
    </tr>
    <tr>
      <th>1</th>
      <td>11899</td>
      <td>-4.66</td>
      <td>C1=CC(=C(C=C1Cl)Cl)[N+](=O)[O-]</td>
      <td>8.134854</td>
      <td>7.770338</td>
      <td>0</td>
      <td>0</td>
      <td>12.675204</td>
      <td>2.302776</td>
      <td>4.605551</td>
      <td>...</td>
      <td>9.094144</td>
      <td>41.023148</td>
      <td>190.954084</td>
      <td>13.639577</td>
      <td>150</td>
      <td>14</td>
      <td>52.0</td>
      <td>58.0</td>
      <td>5.194444</td>
      <td>2.472222</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2723704</td>
      <td>-3.98</td>
      <td>CNC(=S)N</td>
      <td>3.047207</td>
      <td>3.305183</td>
      <td>0</td>
      <td>0</td>
      <td>5.226252</td>
      <td>1.847759</td>
      <td>3.695518</td>
      <td>...</td>
      <td>6.834109</td>
      <td>27.254130</td>
      <td>90.025169</td>
      <td>8.184106</td>
      <td>18</td>
      <td>2</td>
      <td>16.0</td>
      <td>14.0</td>
      <td>3.361111</td>
      <td>1.333333</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3032338</td>
      <td>-4.00</td>
      <td>CCNC(=S)N</td>
      <td>3.754314</td>
      <td>4.057055</td>
      <td>0</td>
      <td>0</td>
      <td>6.155367</td>
      <td>1.902113</td>
      <td>3.804226</td>
      <td>...</td>
      <td>7.131699</td>
      <td>29.439488</td>
      <td>104.040819</td>
      <td>7.431487</td>
      <td>32</td>
      <td>3</td>
      <td>20.0</td>
      <td>18.0</td>
      <td>3.611111</td>
      <td>1.583333</td>
    </tr>
    <tr>
      <th>4</th>
      <td>2346</td>
      <td>-6.54</td>
      <td>C1=CC=C(C=C1)CN=C=S</td>
      <td>7.071068</td>
      <td>6.547760</td>
      <td>0</td>
      <td>0</td>
      <td>12.932143</td>
      <td>2.154341</td>
      <td>4.308683</td>
      <td>...</td>
      <td>8.438366</td>
      <td>38.130322</td>
      <td>149.029920</td>
      <td>8.766466</td>
      <td>133</td>
      <td>9</td>
      <td>42.0</td>
      <td>44.0</td>
      <td>3.111111</td>
      <td>2.500000</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>287</th>
      <td>15292</td>
      <td>-5.65</td>
      <td>CC1=CC(=CC(=C1O)Cl)Cl</td>
      <td>7.427747</td>
      <td>7.188678</td>
      <td>0</td>
      <td>0</td>
      <td>11.643052</td>
      <td>2.307250</td>
      <td>4.614501</td>
      <td>...</td>
      <td>9.078065</td>
      <td>39.748909</td>
      <td>175.979570</td>
      <td>10.998723</td>
      <td>110</td>
      <td>13</td>
      <td>48.0</td>
      <td>54.0</td>
      <td>4.944444</td>
      <td>2.222222</td>
    </tr>
    <tr>
      <th>288</th>
      <td>15767</td>
      <td>-7.01</td>
      <td>COC1=C(C(=C(C(=C1Cl)Cl)Cl)Cl)Cl</td>
      <td>9.496696</td>
      <td>9.522267</td>
      <td>0</td>
      <td>0</td>
      <td>15.688761</td>
      <td>2.426720</td>
      <td>4.853440</td>
      <td>...</td>
      <td>9.602855</td>
      <td>44.616134</td>
      <td>277.862653</td>
      <td>17.366416</td>
      <td>220</td>
      <td>23</td>
      <td>64.0</td>
      <td>77.0</td>
      <td>6.916667</td>
      <td>3.000000</td>
    </tr>
    <tr>
      <th>289</th>
      <td>1551919</td>
      <td>-3.85</td>
      <td>CCCCNC(=S)N</td>
      <td>5.168527</td>
      <td>5.361851</td>
      <td>0</td>
      <td>0</td>
      <td>8.762573</td>
      <td>1.949856</td>
      <td>3.899712</td>
      <td>...</td>
      <td>7.475906</td>
      <td>33.090360</td>
      <td>132.072119</td>
      <td>6.603606</td>
      <td>79</td>
      <td>5</td>
      <td>28.0</td>
      <td>26.0</td>
      <td>4.111111</td>
      <td>2.083333</td>
    </tr>
    <tr>
      <th>290</th>
      <td>2566</td>
      <td>-6.52</td>
      <td>CC1(CC2=C(O1)C(=CC=C2)OC(=O)NC)C</td>
      <td>12.367162</td>
      <td>11.107564</td>
      <td>0</td>
      <td>0</td>
      <td>19.625463</td>
      <td>2.470049</td>
      <td>4.778732</td>
      <td>...</td>
      <td>9.709599</td>
      <td>62.929068</td>
      <td>221.105193</td>
      <td>7.132426</td>
      <td>425</td>
      <td>22</td>
      <td>84.0</td>
      <td>97.0</td>
      <td>6.256944</td>
      <td>3.472222</td>
    </tr>
    <tr>
      <th>291</th>
      <td>13905</td>
      <td>-3.63</td>
      <td>CCNC1=NC(=NC(=N1)SC)NCC</td>
      <td>9.899495</td>
      <td>9.242468</td>
      <td>0</td>
      <td>0</td>
      <td>17.184731</td>
      <td>2.298128</td>
      <td>4.596257</td>
      <td>...</td>
      <td>9.110410</td>
      <td>44.568778</td>
      <td>213.104816</td>
      <td>7.348442</td>
      <td>318</td>
      <td>17</td>
      <td>62.0</td>
      <td>68.0</td>
      <td>5.333333</td>
      <td>3.500000</td>
    </tr>
  </tbody>
</table>
<p>292 rows × 1829 columns</p>
</div>



<h2>1.3. Generating SDF</h2>
<ul>
<li> Generates SDF from SMILES for future use </li>
</ul>


```python
dt.SDFgeneration()
```

    [1mPath :/home/jun/pyqsar/pyqsar3r25[0m
    Index   File Name
       0    qdb.csv
       1    qdb1.csv


    Enter .csv file index to load :  1


    
    
              ID    EP                               SMI
    0      11859 -5.46          C1=C(C=C(C(=C1Cl)Cl)Cl)O
    1      11899 -4.66   C1=CC(=C(C=C1Cl)Cl)[N+](=O)[O-]
    2    2723704 -3.98                          CNC(=S)N
    3    3032338 -4.00                         CCNC(=S)N
    4       2346 -6.54               C1=CC=C(C=C1)CN=C=S
    ..       ...   ...                               ...
    287    15292 -5.65             CC1=CC(=CC(=C1O)Cl)Cl
    288    15767 -7.01   COC1=C(C(=C(C(=C1Cl)Cl)Cl)Cl)Cl
    289  1551919 -3.85                       CCCCNC(=S)N
    290     2566 -6.52  CC1(CC2=C(O1)C(=CC=C2)OC(=O)NC)C
    291    13905 -3.63           CCNC1=NC(=NC(=N1)SC)NCC
    
    [292 rows x 3 columns]


    /home/jun/pyqsar/pyqsar3r25/pyqsar/data_tools.py:600: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      smi_df['Molecule'] = smi_df['SMI'].apply(Chem.MolFromSmiles)


<h2>1.4.  Data Filtering</h2>
<ul>
<li> eliminating non-numerical data </li>
</ul>


```python
X_data = newFormat.iloc[:,3:]
X_data
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ABC</th>
      <th>ABCGG</th>
      <th>nAcid</th>
      <th>nBase</th>
      <th>SpAbs_A</th>
      <th>SpMax_A</th>
      <th>SpDiam_A</th>
      <th>SpAD_A</th>
      <th>SpMAD_A</th>
      <th>LogEE_A</th>
      <th>...</th>
      <th>SRW10</th>
      <th>TSRW10</th>
      <th>MW</th>
      <th>AMW</th>
      <th>WPath</th>
      <th>WPol</th>
      <th>Zagreb1</th>
      <th>Zagreb2</th>
      <th>mZagreb1</th>
      <th>mZagreb2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>7.427747</td>
      <td>7.188678</td>
      <td>0</td>
      <td>0</td>
      <td>11.643052</td>
      <td>2.307250</td>
      <td>4.614501</td>
      <td>11.643052</td>
      <td>1.164305</td>
      <td>3.206363</td>
      <td>...</td>
      <td>9.078065</td>
      <td>39.748909</td>
      <td>195.924948</td>
      <td>15.071150</td>
      <td>110</td>
      <td>13</td>
      <td>48.0</td>
      <td>54.0</td>
      <td>4.944444</td>
      <td>2.222222</td>
    </tr>
    <tr>
      <th>1</th>
      <td>8.134854</td>
      <td>7.770338</td>
      <td>0</td>
      <td>0</td>
      <td>12.675204</td>
      <td>2.302776</td>
      <td>4.605551</td>
      <td>12.675204</td>
      <td>1.152291</td>
      <td>3.294669</td>
      <td>...</td>
      <td>9.094144</td>
      <td>41.023148</td>
      <td>190.954084</td>
      <td>13.639577</td>
      <td>150</td>
      <td>14</td>
      <td>52.0</td>
      <td>58.0</td>
      <td>5.194444</td>
      <td>2.472222</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3.047207</td>
      <td>3.305183</td>
      <td>0</td>
      <td>0</td>
      <td>5.226252</td>
      <td>1.847759</td>
      <td>3.695518</td>
      <td>5.226252</td>
      <td>1.045250</td>
      <td>2.408576</td>
      <td>...</td>
      <td>6.834109</td>
      <td>27.254130</td>
      <td>90.025169</td>
      <td>8.184106</td>
      <td>18</td>
      <td>2</td>
      <td>16.0</td>
      <td>14.0</td>
      <td>3.361111</td>
      <td>1.333333</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3.754314</td>
      <td>4.057055</td>
      <td>0</td>
      <td>0</td>
      <td>6.155367</td>
      <td>1.902113</td>
      <td>3.804226</td>
      <td>6.155367</td>
      <td>1.025895</td>
      <td>2.595100</td>
      <td>...</td>
      <td>7.131699</td>
      <td>29.439488</td>
      <td>104.040819</td>
      <td>7.431487</td>
      <td>32</td>
      <td>3</td>
      <td>20.0</td>
      <td>18.0</td>
      <td>3.611111</td>
      <td>1.583333</td>
    </tr>
    <tr>
      <th>4</th>
      <td>7.071068</td>
      <td>6.547760</td>
      <td>0</td>
      <td>0</td>
      <td>12.932143</td>
      <td>2.154341</td>
      <td>4.308683</td>
      <td>12.932143</td>
      <td>1.293214</td>
      <td>3.179653</td>
      <td>...</td>
      <td>8.438366</td>
      <td>38.130322</td>
      <td>149.029920</td>
      <td>8.766466</td>
      <td>133</td>
      <td>9</td>
      <td>42.0</td>
      <td>44.0</td>
      <td>3.111111</td>
      <td>2.500000</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>287</th>
      <td>7.427747</td>
      <td>7.188678</td>
      <td>0</td>
      <td>0</td>
      <td>11.643052</td>
      <td>2.307250</td>
      <td>4.614501</td>
      <td>11.643052</td>
      <td>1.164305</td>
      <td>3.206363</td>
      <td>...</td>
      <td>9.078065</td>
      <td>39.748909</td>
      <td>175.979570</td>
      <td>10.998723</td>
      <td>110</td>
      <td>13</td>
      <td>48.0</td>
      <td>54.0</td>
      <td>4.944444</td>
      <td>2.222222</td>
    </tr>
    <tr>
      <th>288</th>
      <td>9.496696</td>
      <td>9.522267</td>
      <td>0</td>
      <td>0</td>
      <td>15.688761</td>
      <td>2.426720</td>
      <td>4.853440</td>
      <td>15.688761</td>
      <td>1.206828</td>
      <td>3.465917</td>
      <td>...</td>
      <td>9.602855</td>
      <td>44.616134</td>
      <td>277.862653</td>
      <td>17.366416</td>
      <td>220</td>
      <td>23</td>
      <td>64.0</td>
      <td>77.0</td>
      <td>6.916667</td>
      <td>3.000000</td>
    </tr>
    <tr>
      <th>289</th>
      <td>5.168527</td>
      <td>5.361851</td>
      <td>0</td>
      <td>0</td>
      <td>8.762573</td>
      <td>1.949856</td>
      <td>3.899712</td>
      <td>8.762573</td>
      <td>1.095322</td>
      <td>2.887985</td>
      <td>...</td>
      <td>7.475906</td>
      <td>33.090360</td>
      <td>132.072119</td>
      <td>6.603606</td>
      <td>79</td>
      <td>5</td>
      <td>28.0</td>
      <td>26.0</td>
      <td>4.111111</td>
      <td>2.083333</td>
    </tr>
    <tr>
      <th>290</th>
      <td>12.367162</td>
      <td>11.107564</td>
      <td>0</td>
      <td>0</td>
      <td>19.625463</td>
      <td>2.470049</td>
      <td>4.778732</td>
      <td>19.625463</td>
      <td>1.226591</td>
      <td>3.704782</td>
      <td>...</td>
      <td>9.709599</td>
      <td>62.929068</td>
      <td>221.105193</td>
      <td>7.132426</td>
      <td>425</td>
      <td>22</td>
      <td>84.0</td>
      <td>97.0</td>
      <td>6.256944</td>
      <td>3.472222</td>
    </tr>
    <tr>
      <th>291</th>
      <td>9.899495</td>
      <td>9.242468</td>
      <td>0</td>
      <td>0</td>
      <td>17.184731</td>
      <td>2.298128</td>
      <td>4.596257</td>
      <td>17.184731</td>
      <td>1.227481</td>
      <td>3.514645</td>
      <td>...</td>
      <td>9.110410</td>
      <td>44.568778</td>
      <td>213.104816</td>
      <td>7.348442</td>
      <td>318</td>
      <td>17</td>
      <td>62.0</td>
      <td>68.0</td>
      <td>5.333333</td>
      <td>3.500000</td>
    </tr>
  </tbody>
</table>
<p>292 rows × 1826 columns</p>
</div>




```python
X_data, sel_desc=dt.NonNumricFilter(X_data)
```

    Start :  (292, 1826)
    Filterd : (292, 952)



```python
filtered = newFormat[["ID","EP","SMI"]+sel_desc]
filtered.to_csv("qdb1_F.csv", index=False)
filtered
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ID</th>
      <th>EP</th>
      <th>SMI</th>
      <th>ABC</th>
      <th>ABCGG</th>
      <th>nAcid</th>
      <th>nBase</th>
      <th>SpAbs_A</th>
      <th>SpMax_A</th>
      <th>SpDiam_A</th>
      <th>...</th>
      <th>SRW10</th>
      <th>TSRW10</th>
      <th>MW</th>
      <th>AMW</th>
      <th>WPath</th>
      <th>WPol</th>
      <th>Zagreb1</th>
      <th>Zagreb2</th>
      <th>mZagreb1</th>
      <th>mZagreb2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>11859</td>
      <td>-5.46</td>
      <td>C1=C(C=C(C(=C1Cl)Cl)Cl)O</td>
      <td>7.427747</td>
      <td>7.188678</td>
      <td>0</td>
      <td>0</td>
      <td>11.643052</td>
      <td>2.307250</td>
      <td>4.614501</td>
      <td>...</td>
      <td>9.078065</td>
      <td>39.748909</td>
      <td>195.924948</td>
      <td>15.071150</td>
      <td>110</td>
      <td>13</td>
      <td>48.0</td>
      <td>54.0</td>
      <td>4.944444</td>
      <td>2.222222</td>
    </tr>
    <tr>
      <th>1</th>
      <td>11899</td>
      <td>-4.66</td>
      <td>C1=CC(=C(C=C1Cl)Cl)[N+](=O)[O-]</td>
      <td>8.134854</td>
      <td>7.770338</td>
      <td>0</td>
      <td>0</td>
      <td>12.675204</td>
      <td>2.302776</td>
      <td>4.605551</td>
      <td>...</td>
      <td>9.094144</td>
      <td>41.023148</td>
      <td>190.954084</td>
      <td>13.639577</td>
      <td>150</td>
      <td>14</td>
      <td>52.0</td>
      <td>58.0</td>
      <td>5.194444</td>
      <td>2.472222</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2723704</td>
      <td>-3.98</td>
      <td>CNC(=S)N</td>
      <td>3.047207</td>
      <td>3.305183</td>
      <td>0</td>
      <td>0</td>
      <td>5.226252</td>
      <td>1.847759</td>
      <td>3.695518</td>
      <td>...</td>
      <td>6.834109</td>
      <td>27.254130</td>
      <td>90.025169</td>
      <td>8.184106</td>
      <td>18</td>
      <td>2</td>
      <td>16.0</td>
      <td>14.0</td>
      <td>3.361111</td>
      <td>1.333333</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3032338</td>
      <td>-4.00</td>
      <td>CCNC(=S)N</td>
      <td>3.754314</td>
      <td>4.057055</td>
      <td>0</td>
      <td>0</td>
      <td>6.155367</td>
      <td>1.902113</td>
      <td>3.804226</td>
      <td>...</td>
      <td>7.131699</td>
      <td>29.439488</td>
      <td>104.040819</td>
      <td>7.431487</td>
      <td>32</td>
      <td>3</td>
      <td>20.0</td>
      <td>18.0</td>
      <td>3.611111</td>
      <td>1.583333</td>
    </tr>
    <tr>
      <th>4</th>
      <td>2346</td>
      <td>-6.54</td>
      <td>C1=CC=C(C=C1)CN=C=S</td>
      <td>7.071068</td>
      <td>6.547760</td>
      <td>0</td>
      <td>0</td>
      <td>12.932143</td>
      <td>2.154341</td>
      <td>4.308683</td>
      <td>...</td>
      <td>8.438366</td>
      <td>38.130322</td>
      <td>149.029920</td>
      <td>8.766466</td>
      <td>133</td>
      <td>9</td>
      <td>42.0</td>
      <td>44.0</td>
      <td>3.111111</td>
      <td>2.500000</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>287</th>
      <td>15292</td>
      <td>-5.65</td>
      <td>CC1=CC(=CC(=C1O)Cl)Cl</td>
      <td>7.427747</td>
      <td>7.188678</td>
      <td>0</td>
      <td>0</td>
      <td>11.643052</td>
      <td>2.307250</td>
      <td>4.614501</td>
      <td>...</td>
      <td>9.078065</td>
      <td>39.748909</td>
      <td>175.979570</td>
      <td>10.998723</td>
      <td>110</td>
      <td>13</td>
      <td>48.0</td>
      <td>54.0</td>
      <td>4.944444</td>
      <td>2.222222</td>
    </tr>
    <tr>
      <th>288</th>
      <td>15767</td>
      <td>-7.01</td>
      <td>COC1=C(C(=C(C(=C1Cl)Cl)Cl)Cl)Cl</td>
      <td>9.496696</td>
      <td>9.522267</td>
      <td>0</td>
      <td>0</td>
      <td>15.688761</td>
      <td>2.426720</td>
      <td>4.853440</td>
      <td>...</td>
      <td>9.602855</td>
      <td>44.616134</td>
      <td>277.862653</td>
      <td>17.366416</td>
      <td>220</td>
      <td>23</td>
      <td>64.0</td>
      <td>77.0</td>
      <td>6.916667</td>
      <td>3.000000</td>
    </tr>
    <tr>
      <th>289</th>
      <td>1551919</td>
      <td>-3.85</td>
      <td>CCCCNC(=S)N</td>
      <td>5.168527</td>
      <td>5.361851</td>
      <td>0</td>
      <td>0</td>
      <td>8.762573</td>
      <td>1.949856</td>
      <td>3.899712</td>
      <td>...</td>
      <td>7.475906</td>
      <td>33.090360</td>
      <td>132.072119</td>
      <td>6.603606</td>
      <td>79</td>
      <td>5</td>
      <td>28.0</td>
      <td>26.0</td>
      <td>4.111111</td>
      <td>2.083333</td>
    </tr>
    <tr>
      <th>290</th>
      <td>2566</td>
      <td>-6.52</td>
      <td>CC1(CC2=C(O1)C(=CC=C2)OC(=O)NC)C</td>
      <td>12.367162</td>
      <td>11.107564</td>
      <td>0</td>
      <td>0</td>
      <td>19.625463</td>
      <td>2.470049</td>
      <td>4.778732</td>
      <td>...</td>
      <td>9.709599</td>
      <td>62.929068</td>
      <td>221.105193</td>
      <td>7.132426</td>
      <td>425</td>
      <td>22</td>
      <td>84.0</td>
      <td>97.0</td>
      <td>6.256944</td>
      <td>3.472222</td>
    </tr>
    <tr>
      <th>291</th>
      <td>13905</td>
      <td>-3.63</td>
      <td>CCNC1=NC(=NC(=N1)SC)NCC</td>
      <td>9.899495</td>
      <td>9.242468</td>
      <td>0</td>
      <td>0</td>
      <td>17.184731</td>
      <td>2.298128</td>
      <td>4.596257</td>
      <td>...</td>
      <td>9.110410</td>
      <td>44.568778</td>
      <td>213.104816</td>
      <td>7.348442</td>
      <td>318</td>
      <td>17</td>
      <td>62.0</td>
      <td>68.0</td>
      <td>5.333333</td>
      <td>3.500000</td>
    </tr>
  </tbody>
</table>
<p>292 rows × 955 columns</p>
</div>



## Check if there are descriptors with one value


```python
dt.uCheck(filtered)
```


    
![png](output_19_0.png)
    


<h2>1.5.  Data scaling</h2>
<h3> ScalingTools(scale, ext='.csv') </h3>
<ul>
<li> scale : Scaling Method, 'minmax','standard' </li>
<li> ext : call csv format file with ext extension' </li>
</ul>


```python
st = dt.ScalingTools(scale='standard')
```

    [1mPath :/home/jun/pyqsar/pyqsar3r25[0m
    Index   File Name
       0    qdb.csv
       1    qdb1.csv
       2    qdb1_F.csv


    Enter .csv file index to load :  2


    
    


<h3> train_scaler() </h3>
<ul>
<li> * fit and transform data </li>
<li> * save scaler file to use in model information </li>
</ul>    


```python
st.train_scaler()
```

    Scaler File : qdb1_F_standard.info file saved





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ID</th>
      <th>EP</th>
      <th>SMI</th>
      <th>ABC</th>
      <th>ABCGG</th>
      <th>nAcid</th>
      <th>nBase</th>
      <th>SpAbs_A</th>
      <th>SpMax_A</th>
      <th>SpDiam_A</th>
      <th>...</th>
      <th>SRW10</th>
      <th>TSRW10</th>
      <th>MW</th>
      <th>AMW</th>
      <th>WPath</th>
      <th>WPol</th>
      <th>Zagreb1</th>
      <th>Zagreb2</th>
      <th>mZagreb1</th>
      <th>mZagreb2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>11859</td>
      <td>-5.46</td>
      <td>C1=C(C=C(C(=C1Cl)Cl)Cl)O</td>
      <td>-0.151062</td>
      <td>-0.093948</td>
      <td>-0.117851</td>
      <td>-0.162791</td>
      <td>-0.257506</td>
      <td>0.465983</td>
      <td>0.541285</td>
      <td>...</td>
      <td>0.420833</td>
      <td>-0.102373</td>
      <td>0.159083</td>
      <td>0.948995</td>
      <td>-0.374862</td>
      <td>-0.050123</td>
      <td>-0.111334</td>
      <td>-0.087220</td>
      <td>0.080589</td>
      <td>-0.324793</td>
    </tr>
    <tr>
      <th>1</th>
      <td>11899</td>
      <td>-4.66</td>
      <td>C1=CC(=C(C=C1Cl)Cl)[N+](=O)[O-]</td>
      <td>0.008807</td>
      <td>0.063539</td>
      <td>-0.117851</td>
      <td>-0.162791</td>
      <td>-0.112802</td>
      <td>0.449525</td>
      <td>0.524199</td>
      <td>...</td>
      <td>0.432253</td>
      <td>-0.000520</td>
      <td>0.102523</td>
      <td>0.676300</td>
      <td>-0.284085</td>
      <td>0.044915</td>
      <td>0.016656</td>
      <td>0.015890</td>
      <td>0.196860</td>
      <td>-0.113150</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2723704</td>
      <td>-3.98</td>
      <td>CNC(=S)N</td>
      <td>-1.141451</td>
      <td>-1.145427</td>
      <td>-0.117851</td>
      <td>-0.162791</td>
      <td>-1.157113</td>
      <td>-1.223984</td>
      <td>-1.213185</td>
      <td>...</td>
      <td>-1.173080</td>
      <td>-1.101110</td>
      <td>-1.045867</td>
      <td>-0.362894</td>
      <td>-0.583651</td>
      <td>-1.095539</td>
      <td>-1.135254</td>
      <td>-1.118326</td>
      <td>-0.655795</td>
      <td>-1.077300</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3032338</td>
      <td>-4.00</td>
      <td>CCNC(=S)N</td>
      <td>-0.981582</td>
      <td>-0.941853</td>
      <td>-0.117851</td>
      <td>-0.162791</td>
      <td>-1.026855</td>
      <td>-1.024075</td>
      <td>-1.005646</td>
      <td>...</td>
      <td>-0.961698</td>
      <td>-0.926429</td>
      <td>-0.886394</td>
      <td>-0.506258</td>
      <td>-0.551879</td>
      <td>-1.000502</td>
      <td>-1.007264</td>
      <td>-1.015215</td>
      <td>-0.539523</td>
      <td>-0.865657</td>
    </tr>
    <tr>
      <th>4</th>
      <td>2346</td>
      <td>-6.54</td>
      <td>C1=CC=C(C=C1)CN=C=S</td>
      <td>-0.231703</td>
      <td>-0.267480</td>
      <td>-0.117851</td>
      <td>-0.162791</td>
      <td>-0.076781</td>
      <td>-0.096403</td>
      <td>-0.042566</td>
      <td>...</td>
      <td>-0.033554</td>
      <td>-0.231750</td>
      <td>-0.374499</td>
      <td>-0.251963</td>
      <td>-0.322665</td>
      <td>-0.430274</td>
      <td>-0.303319</td>
      <td>-0.344997</td>
      <td>-0.772066</td>
      <td>-0.089634</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>287</th>
      <td>15292</td>
      <td>-5.65</td>
      <td>CC1=CC(=CC(=C1O)Cl)Cl</td>
      <td>-0.151062</td>
      <td>-0.093948</td>
      <td>-0.117851</td>
      <td>-0.162791</td>
      <td>-0.257506</td>
      <td>0.465983</td>
      <td>0.541285</td>
      <td>...</td>
      <td>0.420833</td>
      <td>-0.102373</td>
      <td>-0.067860</td>
      <td>0.173252</td>
      <td>-0.374862</td>
      <td>-0.050123</td>
      <td>-0.111334</td>
      <td>-0.087220</td>
      <td>0.080589</td>
      <td>-0.324793</td>
    </tr>
    <tr>
      <th>288</th>
      <td>15767</td>
      <td>-7.01</td>
      <td>COC1=C(C(=C(C(=C1Cl)Cl)Cl)Cl)Cl</td>
      <td>0.316703</td>
      <td>0.537885</td>
      <td>-0.117851</td>
      <td>-0.162791</td>
      <td>0.309685</td>
      <td>0.905381</td>
      <td>0.997454</td>
      <td>...</td>
      <td>0.793598</td>
      <td>0.286676</td>
      <td>1.091388</td>
      <td>1.386213</td>
      <td>-0.125224</td>
      <td>0.900256</td>
      <td>0.400626</td>
      <td>0.505665</td>
      <td>0.997839</td>
      <td>0.333651</td>
    </tr>
    <tr>
      <th>289</th>
      <td>1551919</td>
      <td>-3.85</td>
      <td>CCCCNC(=S)N</td>
      <td>-0.661845</td>
      <td>-0.588572</td>
      <td>-0.117851</td>
      <td>-0.162791</td>
      <td>-0.661337</td>
      <td>-0.848482</td>
      <td>-0.823350</td>
      <td>...</td>
      <td>-0.717203</td>
      <td>-0.634606</td>
      <td>-0.567448</td>
      <td>-0.663959</td>
      <td>-0.445215</td>
      <td>-0.810426</td>
      <td>-0.751284</td>
      <td>-0.808994</td>
      <td>-0.306981</td>
      <td>-0.442372</td>
    </tr>
    <tr>
      <th>290</th>
      <td>2566</td>
      <td>-6.52</td>
      <td>CC1(CC2=C(O1)C(=CC=C2)OC(=O)NC)C</td>
      <td>0.965682</td>
      <td>0.967113</td>
      <td>-0.117851</td>
      <td>-0.162791</td>
      <td>0.861594</td>
      <td>1.064743</td>
      <td>0.854826</td>
      <td>...</td>
      <td>0.869420</td>
      <td>1.750471</td>
      <td>0.445589</td>
      <td>-0.563225</td>
      <td>0.340012</td>
      <td>0.805218</td>
      <td>1.040576</td>
      <td>1.021218</td>
      <td>0.691012</td>
      <td>0.733421</td>
    </tr>
    <tr>
      <th>291</th>
      <td>13905</td>
      <td>-3.63</td>
      <td>CCNC1=NC(=NC(=N1)SC)NCC</td>
      <td>0.407771</td>
      <td>0.462127</td>
      <td>-0.117851</td>
      <td>-0.162791</td>
      <td>0.519414</td>
      <td>0.432433</td>
      <td>0.506455</td>
      <td>...</td>
      <td>0.443807</td>
      <td>0.282890</td>
      <td>0.354559</td>
      <td>-0.522077</td>
      <td>0.097182</td>
      <td>0.330029</td>
      <td>0.336631</td>
      <td>0.273667</td>
      <td>0.261455</td>
      <td>0.756937</td>
    </tr>
  </tbody>
</table>
<p>292 rows × 955 columns</p>
</div>



<h2>1.6.  Outlier Check - You can Skip this part!</h2>
<h3> 1.6.1. check_outlier(standard) - You can Skip! </h3>
<ul>
<li> * Check outlier value in data </li>
<li> * See with histogram, and delete if you want </li>
</ul>    


```python
data = st.check_outlier(9)
```


    interactive(children=(Dropdown(description='Descriptor', options=('nBase', 'VR1_A', 'VR2_A', 'nBridgehead', 'n…


<h3> 1.6.2. dt.outplot(figsize=(2000,2000)) - You can skip! </h3>

* Open File which ext is .out, show histogram of outlier molecule and scatter plot 
* Can Choose File by input index number 
* figsize= (default) (2000,2000) : (weight, height) 
* Return 
   output, dictionary type data   
 ** To import image about outlier index plot, requires the kaleido package : pip install -U kaleido 



```python
st.outplot()
```

    /home/jun/pyqsar/pyqsar3r25
    File Index   File Name
    0            qdb1_F_standard_9.out


    Enter File Index Num of Above List :  0


    qdb1_F_standard_9.out is selected



    interactive(children=(Dropdown(description='ID:', options=('11712', '12358480', '1269845', '12901', '13218779'…


<h2> 1.7. Save Scaled Data </h2>


```python
st.save()
```

    File name to save scaling as
    (default) qdb1_F_s.csv
    
    - 


    qdb1_F_s.csv is saved
    
    


<h2>1.8. Data Split</h2>
<h3> split data into training and test sets </h3>
<h3> data_split('.csv, test_num, Block) </h3>
    * test_num : number of test molecules, the rest of them is training date <br>
    * Random : True : select test_num molecues randomly<br>
               False : take test_num molecues from the end of the data 
                


```python
dt.data_split('.csv',test_num=98, Random=True)  # false는 뒤에서 가져온다 
```

    [1mPath :/home/jun/pyqsar/pyqsar3r25[0m
    Index   File Name
       0    qdb.csv
       1    qdb1.csv
       2    qdb1_F.csv
       3    qdb1_F_s.csv


    Enter .csv file index to load :  3


    
    
    drop_index : ['NaasN', 'SaasN', 'n6AHRing', 'n8FRing', 'n8FHRing', 'n8FARing', 'n8FAHRing']
    try   1
    ['NaasN', 'SaasN', 'n6AHRing', 'n8FRing', 'n8FHRing', 'n8FARing', 'n8FAHRing']
    --------------------
    drop_index : ['NddssS', 'SddssS', 'n6AHRing']
    try   2
    ['NddssS', 'SddssS', 'n6AHRing']
    --------------------
    drop_index : ['nI', 'C2SP1', 'NsI', 'SsI']
    try   3
    ['nI', 'C2SP1', 'NsI', 'SsI']
    --------------------
    drop_index : ['nI', 'C2SP1', 'NaaO', 'NsI', 'SaaO', 'SsI', 'n6AHRing', 'n7FRing', 'n7FARing', 'n9FARing', 'n9FAHRing']
    try   4
    ['nI', 'C2SP1', 'NaaO', 'NsI', 'SaaO', 'SsI', 'n6AHRing', 'n7FRing', 'n7FARing', 'n9FARing', 'n9FAHRing']
    --------------------
    drop_index : ['nI', 'C2SP1', 'NsI', 'SsI', 'n7FRing', 'n7FARing']
    try   5
    ['nI', 'C2SP1', 'NsI', 'SsI', 'n7FRing', 'n7FARing']
    --------------------
    drop_index : ['n7FRing', 'n8FRing', 'n8FHRing', 'n7FARing', 'n8FARing', 'n8FAHRing']
    try   6
    ['n7FRing', 'n8FRing', 'n8FHRing', 'n7FARing', 'n8FARing', 'n8FAHRing']
    --------------------
    drop_index : ['nI', 'C2SP1', 'NsI', 'SsI', 'n11FRing', 'n11FHRing', 'n11FARing', 'n11FAHRing']
    try   7
    ['nI', 'C2SP1', 'NsI', 'SsI', 'n11FRing', 'n11FHRing', 'n11FARing', 'n11FAHRing']
    --------------------
    drop_index : ['NaaO', 'NdssS', 'SaaO', 'SdssS', 'n7Ring', 'n7HRing', 'n7ARing', 'n7AHRing', 'n12FRing', 'n10FHRing', 'n12FHRing', 'n9FaRing', 'n9FaHRing', 'n10FaHRing', 'n12FARing', 'n12FAHRing']
    try   8
    ['NaaO', 'NdssS', 'SaaO', 'SdssS', 'n7Ring', 'n7HRing', 'n7ARing', 'n7AHRing', 'n12FRing', 'n10FHRing', 'n12FHRing', 'n9FaRing', 'n9FaHRing', 'n10FaHRing', 'n12FARing', 'n12FAHRing']
    --------------------
    drop_index : ['NdssS', 'SdssS', 'n7Ring', 'n7HRing', 'n7ARing', 'n7AHRing', 'n11FRing', 'n12FRing', 'n11FHRing', 'n12FHRing', 'nG12FaHRing', 'n11FARing', 'n12FARing', 'n11FAHRing', 'n12FAHRing']
    try   9
    ['NdssS', 'SdssS', 'n7Ring', 'n7HRing', 'n7ARing', 'n7AHRing', 'n11FRing', 'n12FRing', 'n11FHRing', 'n12FHRing', 'nG12FaHRing', 'n11FARing', 'n12FARing', 'n11FAHRing', 'n12FAHRing']
    --------------------
    drop_index : ['n7FRing', 'n8FRing', 'n11FRing', 'n8FHRing', 'n11FHRing', 'n7FARing', 'n8FARing', 'n11FARing', 'n8FAHRing', 'n11FAHRing']
    try  10
    ['n7FRing', 'n8FRing', 'n11FRing', 'n8FHRing', 'n11FHRing', 'n7FARing', 'n8FARing', 'n11FARing', 'n8FAHRing', 'n11FAHRing']
    train :  (194, 945)
    qdb1_F_s.train  is saved.
    test :  (98, 955)


    /home/jun/pyqsar/pyqsar3r25/pyqsar/data_tools.py:520: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      train_frame_scaled.insert(0, "ID", _id.iloc[train_index].values)
    /home/jun/pyqsar/pyqsar3r25/pyqsar/data_tools.py:521: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      train_frame_scaled.insert(1, "EP", _ep.iloc[train_index].values)
    /home/jun/pyqsar/pyqsar3r25/pyqsar/data_tools.py:522: PerformanceWarning: DataFrame is highly fragmented.  This is usually the result of calling `frame.insert` many times, which has poor performance.  Consider joining all columns at once using pd.concat(axis=1) instead. To get a de-fragmented frame, use `newframe = frame.copy()`
      train_frame_scaled.insert(2, "SMI", _smi.iloc[train_index].values)


    qdb1_F_s.test   is saved.


<h1> 2. Clustering</h1>
<h2> 2.1. Data Load </h2>
 * Split data by [descriptor, EP] set from train/test files 
   


```python
# You can restart from here
from pyqsar import data_tools as dt
from pyqsar import model_tools as mt
from pyqsar import draw_mol
import pandas as pd
import numpy as np
```


```python
X_train, y_train = mt.split_xy('.train')
```

    [1mPath :/home/jun/pyqsar/pyqsar3r25[0m
    Index   File Name
       0    qdb1_F_s.train
       1    qdb1_F_s_B.train


    Enter .train file index to load :  0


    
    


<h1> 2.2-1. KMeans</h1>
<h2> mt.FeatureCluster_KMeans(train_Data, K value=500, init='k-means++', algorithm='auto') : Apply KMeans Class </h2>

    * qstart = saved file name of clust info
    * K value = # cluster (default) 500
    * init = mthod of finding centroid 
        Possible : (default)'k-means++', 'random' = select centroid randomly in rows
    * algorithm
        Possible : (default)'auto', 'elkan', 'full'
            full : classical EM(expectaion maximization) algorithm => appropriate sparse data
            elkan : appropriate dense data
            auto : automatically select according to the data, usually choose elkan
    clust = fit KMeans(train_data, K) Model 
    set_cluster = get cluster value as array
    cluster_dist = clust evaluation by histogram
    model_evaluation = clust evaluation by silhouette score
         0 : cluster overlapped
         1 : clustering well
        -1 : clustering bad  


```python
dt.man('kmeans.png')
```


    
![png](output_36_0.png)
    



```python
# KMeans Clustering 
clust = mt.FeatureCluster_KMeans(X_train,400,init='k-means++',algorithm='auto')
clust_info = clust.set_cluster()
```

    
     [1;46mCluster[0m 0 ['ATSC5s'] 
     [1;46mCluster[0m 1 ['ATS4d', 'ATS4se', 'ATS4pe', 'ATS4are', 'ATS4i'] 
     [1;46mCluster[0m 2 ['ATS3Z', 'ATS4Z', 'ATS3m', 'ATS4m'] 
     [1;46mCluster[0m 3 ['ATSC8s'] 
     [1;46mCluster[0m 4 ['VR1_A', 'VR2_A', 'C4SP3', 'Xch-4d', 'Xch-5d', 'Xch-4dv', 'Xch-5dv'] 
     [1;46mCluster[0m 5 ['nO', 'SMR_VSA1', 'MID_O'] 
     [1;46mCluster[0m 6 ['piPC2', 'piPC3', 'piPC4', 'piPC5', 'piPC6', 'TpiPC10'] 
     [1;46mCluster[0m 7 ['VE2_A', 'VE2_DzZ', 'VE2_Dzm', 'VE2_Dzv', 'VE2_Dzse', 'VE2_Dzpe', 'VE2_Dzare', 'VE2_Dzp', 'VE2_Dzi', 'VE2_D'] 
     [1;46mCluster[0m 8 ['AATS1Z', 'AATS1m', 'MZ', 'Mm', 'AMW'] 
     [1;46mCluster[0m 9 ['SpAbs_DzZ', 'SpMax_DzZ', 'SpDiam_DzZ', 'SpAD_DzZ', 'LogEE_DzZ', 'SpAbs_Dzm', 'SpMax_Dzm', 'SpDiam_Dzm', 'SpAD_Dzm', 'LogEE_Dzm', 'SpAbs_D', 'SpMax_D', 'SpDiam_D', 'SpAD_D', 'LogEE_D', 'ECIndex', 'WPath'] 
     [1;46mCluster[0m 10 ['nBridgehead', 'C3SP3', 'Xch-7dv'] 
     [1;46mCluster[0m 11 ['LogEE_A', 'VE1_A', 'VE3_A', 'VE3_DzZ', 'VR3_DzZ', 'VE3_Dzm', 'VR3_Dzm', 'VE3_Dzv', 'VR3_Dzv', 'VE3_Dzse', 'VR3_Dzse', 'VE3_Dzpe', 'VR3_Dzpe', 'VE3_Dzare', 'VR3_Dzare', 'VE3_Dzp', 'VR3_Dzp', 'VE3_Dzi', 'VR3_Dzi', 'VE3_D', 'VR3_D', 'piPC1', 'VAdjMat', 'MWC02', 'SRW02'] 
     [1;46mCluster[0m 12 ['MATS1se', 'MATS1pe', 'MATS1are'] 
     [1;46mCluster[0m 13 ['AATS0i', 'AATS1i', 'Mi'] 
     [1;46mCluster[0m 14 ['ETA_epsilon_2', 'ETA_epsilon_5'] 
     [1;46mCluster[0m 15 ['AATSC2Z', 'AATSC2m'] 
     [1;46mCluster[0m 16 ['n10FRing', 'n10FaRing'] 
     [1;46mCluster[0m 17 ['nN', 'MID_N'] 
     [1;46mCluster[0m 18 ['CIC2'] 
     [1;46mCluster[0m 19 ['MIC2', 'MIC3', 'MIC4', 'MIC5'] 
     [1;46mCluster[0m 20 ['SMR_VSA5', 'SlogP_VSA5'] 
     [1;46mCluster[0m 21 ['NtsC', 'StsC'] 
     [1;46mCluster[0m 22 ['naHRing'] 
     [1;46mCluster[0m 23 ['C2SP2', 'SMR_VSA7'] 
     [1;46mCluster[0m 24 ['NdssS', 'n7Ring', 'n7HRing', 'n7ARing', 'n7AHRing', 'n12FRing', 'n12FHRing', 'n12FARing', 'n12FAHRing'] 
     [1;46mCluster[0m 25 ['SM1_Dzv', 'SM1_Dzp'] 
     [1;46mCluster[0m 26 ['NaaO', 'SaaO'] 
     [1;46mCluster[0m 27 ['ATSC6Z', 'ATSC6m'] 
     [1;46mCluster[0m 28 ['ATSC4se', 'ATSC4pe', 'ATSC4are'] 
     [1;46mCluster[0m 29 ['NaasN', 'SaasN'] 
     [1;46mCluster[0m 30 ['ATSC0p'] 
     [1;46mCluster[0m 31 ['ATSC5v', 'ATSC5p'] 
     [1;46mCluster[0m 32 ['AATS0s', 'AATS2s', 'AATSC0s'] 
     [1;46mCluster[0m 33 ['NaaaC', 'SaaaC', 'nFaRing'] 
     [1;46mCluster[0m 34 ['IC3', 'IC4', 'IC5'] 
     [1;46mCluster[0m 35 ['ATSC3pe', 'ATSC3are'] 
     [1;46mCluster[0m 36 ['NdsN', 'SdsN'] 
     [1;46mCluster[0m 37 ['ATSC2v', 'ATSC2p'] 
     [1;46mCluster[0m 38 ['n9FRing', 'n9FHRing'] 
     [1;46mCluster[0m 39 ['ATS0v', 'ATS1v', 'ATS4v', 'ATS0p', 'ATS1p', 'ATS2p', 'ATS3p', 'ATS4p', 'Xp-1dv', 'Sv', 'Sp', 'ETA_alpha', 'VMcGowan', 'apol', 'SMR'] 
     [1;46mCluster[0m 40 ['ATSC7se', 'ATSC7pe', 'ATSC7are'] 
     [1;46mCluster[0m 41 ['MATS2s'] 
     [1;46mCluster[0m 42 ['NsNH2', 'SsNH2', 'VSA_EState4'] 
     [1;46mCluster[0m 43 ['NsOH', 'SsOH'] 
     [1;46mCluster[0m 44 ['BCUTdv-1l'] 
     [1;46mCluster[0m 45 ['AATS1pe', 'AATS0are', 'AATS1are', 'AATS2are', 'Mare'] 
     [1;46mCluster[0m 46 ['GATS1Z', 'GATS1m'] 
     [1;46mCluster[0m 47 ['ATSC5se', 'ATSC5pe', 'ATSC5are'] 
     [1;46mCluster[0m 48 ['SIC3', 'SIC4', 'SIC5', 'BIC3', 'BIC4', 'BIC5'] 
     [1;46mCluster[0m 49 ['SIC1', 'BIC1'] 
     [1;46mCluster[0m 50 ['SM1_DzZ', 'SM1_Dzm'] 
     [1;46mCluster[0m 51 ['ATSC6dv'] 
     [1;46mCluster[0m 52 ['nF', 'NsF', 'SsF'] 
     [1;46mCluster[0m 53 ['ATSC8v'] 
     [1;46mCluster[0m 54 ['ATSC3v'] 
     [1;46mCluster[0m 55 ['AXp-0dv', 'AXp-1dv'] 
     [1;46mCluster[0m 56 ['AATSC1i', 'MATS1i'] 
     [1;46mCluster[0m 57 ['nP', 'NdsssP'] 
     [1;46mCluster[0m 58 ['GATS2Z', 'GATS2m'] 
     [1;46mCluster[0m 59 ['BCUTpe-1l', 'BCUTare-1l'] 
     [1;46mCluster[0m 60 ['NdsCH', 'SdsCH'] 
     [1;46mCluster[0m 61 ['TopoShapeIndex', 'PetitjeanIndex'] 
     [1;46mCluster[0m 62 ['NsSH', 'SsSH'] 
     [1;46mCluster[0m 63 ['n9FaRing', 'n9FaHRing'] 
     [1;46mCluster[0m 64 ['nX', 'MID_X'] 
     [1;46mCluster[0m 65 ['nBr', 'NsBr', 'SsBr', 'ETA_dPsi_B'] 
     [1;46mCluster[0m 66 ['NaaNH', 'SaaNH'] 
     [1;46mCluster[0m 67 ['ATSC6se', 'ATSC6pe', 'ATSC6are'] 
     [1;46mCluster[0m 68 ['MATS1Z', 'MATS1m'] 
     [1;46mCluster[0m 69 ['GATS2se', 'GATS2pe', 'GATS2are'] 
     [1;46mCluster[0m 70 ['n5aRing'] 
     [1;46mCluster[0m 71 ['ATS6d', 'ATS6v', 'ATS7v', 'ATS6p', 'ATS7p'] 
     [1;46mCluster[0m 72 ['ATSC0Z', 'ATSC0m'] 
     [1;46mCluster[0m 73 ['Xch-3d', 'Xch-3dv', 'n3Ring', 'n3ARing', 'SRW03'] 
     [1;46mCluster[0m 74 ['ATSC2se', 'ATSC2pe', 'ATSC2are'] 
     [1;46mCluster[0m 75 ['GATS1se', 'GATS1pe', 'GATS1are'] 
     [1;46mCluster[0m 76 ['AATS0v', 'AATS1v', 'AATS2v', 'Mv'] 
     [1;46mCluster[0m 77 ['NdCH2', 'SdCH2'] 
     [1;46mCluster[0m 78 ['nI', 'C2SP1', 'NsI', 'SsI'] 
     [1;46mCluster[0m 79 ['n9FARing', 'n9FAHRing'] 
     [1;46mCluster[0m 80 ['SssssC'] 
     [1;46mCluster[0m 81 ['nARing'] 
     [1;46mCluster[0m 82 ['SM1_Dzse', 'SM1_Dzpe'] 
     [1;46mCluster[0m 83 ['NaaN', 'SaaN'] 
     [1;46mCluster[0m 84 ['AETA_eta_B', 'AETA_eta_BR', 'JGI1'] 
     [1;46mCluster[0m 85 ['SRW05', 'SRW07', 'SRW09'] 
     [1;46mCluster[0m 86 ['Xc-6d', 'Xc-5dv', 'Xc-6dv', 'NssssC'] 
     [1;46mCluster[0m 87 ['ATSC4p'] 
     [1;46mCluster[0m 88 ['C1SP2'] 
     [1;46mCluster[0m 89 ['SpMax_A', 'SpDiam_A', 'MWC03', 'MWC04', 'MWC05', 'MWC06', 'MWC07', 'MWC08', 'MWC09', 'MWC10', 'SRW04', 'SRW06', 'SRW08', 'SRW10'] 
     [1;46mCluster[0m 90 ['SlogP_VSA2'] 
     [1;46mCluster[0m 91 ['NssNH', 'SssNH'] 
     [1;46mCluster[0m 92 ['AATSC1dv', 'MATS1dv'] 
     [1;46mCluster[0m 93 ['SpAbs_A', 'SpAD_A', 'nHeavyAtom', 'ATS0d', 'VR1_DzZ', 'VR1_Dzm', 'VR1_Dzv', 'VR1_Dzse', 'VR1_Dzpe', 'VR1_Dzare', 'VR1_Dzp', 'VR1_Dzi', 'nBondsO', 'Xp-0d', 'Xp-1d', 'Xp-2d', 'VR1_D', 'ETA_beta_s', 'ETA_eta_R', 'ETA_eta_RL', 'MID', 'MWC01'] 
     [1;46mCluster[0m 94 ['MPC8', 'MPC9', 'MPC10'] 
     [1;46mCluster[0m 95 ['ATS1Z', 'ATS1m', 'Xp-0dv', 'SZ', 'Sm', 'ZMIC0', 'LabuteASA', 'MW'] 
     [1;46mCluster[0m 96 ['ATSC2c'] 
     [1;46mCluster[0m 97 ['Xp-5d', 'Xp-6d', 'Xp-7d', 'MPC4', 'MPC5'] 
     [1;46mCluster[0m 98 ['AETA_eta_RL'] 
     [1;46mCluster[0m 99 ['BCUTi-1l'] 
     [1;46mCluster[0m 100 ['AETA_beta', 'AETA_beta_ns', 'AETA_eta_FL', 'AETA_dBeta'] 
     [1;46mCluster[0m 101 ['BCUTZ-1l', 'BCUTm-1l', 'BCUTse-1l'] 
     [1;46mCluster[0m 102 ['AATS0dv', 'AATS1dv', 'AATS2dv'] 
     [1;46mCluster[0m 103 ['ATSC2i', 'AATSC2i', 'MATS2i'] 
     [1;46mCluster[0m 104 ['NssS', 'SssS'] 
     [1;46mCluster[0m 105 ['ATSC0c'] 
     [1;46mCluster[0m 106 ['NaaS', 'SaaS'] 
     [1;46mCluster[0m 107 ['n10FHRing', 'n10FaHRing'] 
     [1;46mCluster[0m 108 ['NsssN', 'SsssN'] 
     [1;46mCluster[0m 109 ['AATSC2se', 'AATSC2pe', 'AATSC2are'] 
     [1;46mCluster[0m 110 ['NddC', 'SddC'] 
     [1;46mCluster[0m 111 ['BCUTv-1l', 'BCUTp-1l'] 
     [1;46mCluster[0m 112 ['nBase'] 
     [1;46mCluster[0m 113 ['ATSC7Z', 'ATSC7m'] 
     [1;46mCluster[0m 114 ['ATSC1se', 'ATSC1pe', 'ATSC1are'] 
     [1;46mCluster[0m 115 ['ATS0s', 'ATSC0s', 'EState_VSA10'] 
     [1;46mCluster[0m 116 ['SlogP_VSA7'] 
     [1;46mCluster[0m 117 ['NdS', 'SdS'] 
     [1;46mCluster[0m 118 ['EState_VSA8'] 
     [1;46mCluster[0m 119 ['nG12FRing', 'nG12FHRing'] 
     [1;46mCluster[0m 120 ['EState_VSA7'] 
     [1;46mCluster[0m 121 ['GGI9', 'GGI10'] 
     [1;46mCluster[0m 122 ['SlogP_VSA11'] 
     [1;46mCluster[0m 123 ['BCUTare-1h'] 
     [1;46mCluster[0m 124 ['Xp-5dv', 'Xp-6dv', 'Xp-7dv'] 
     [1;46mCluster[0m 125 ['NssCH2', 'SssCH2'] 
     [1;46mCluster[0m 126 ['ATSC5Z', 'ATSC5m'] 
     [1;46mCluster[0m 127 ['ATSC1s'] 
     [1;46mCluster[0m 128 ['nAtom', 'ATS0se', 'ATS1se', 'ATS2se', 'ATS3se', 'ATS0pe', 'ATS1pe', 'ATS2pe', 'ATS3pe', 'ATS0are', 'ATS1are', 'ATS2are', 'ATS3are', 'ATS0i', 'ATS1i', 'ATS2i', 'ATS3i', 'ATSC0v', 'nBonds', 'nBondsKS', 'Sse', 'Spe', 'Sare', 'Si', 'TIC0'] 
     [1;46mCluster[0m 129 ['EState_VSA5'] 
     [1;46mCluster[0m 130 ['PEOE_VSA5'] 
     [1;46mCluster[0m 131 ['AMID_C'] 
     [1;46mCluster[0m 132 ['ETA_dEpsilon_D'] 
     [1;46mCluster[0m 133 ['VSA_EState5'] 
     [1;46mCluster[0m 134 ['C1SP3'] 
     [1;46mCluster[0m 135 ['ZMIC3', 'ZMIC4', 'ZMIC5'] 
     [1;46mCluster[0m 136 ['ATSC4d'] 
     [1;46mCluster[0m 137 ['ATSC3Z', 'ATSC3m'] 
     [1;46mCluster[0m 138 ['ATSC1p', 'AATSC1p'] 
     [1;46mCluster[0m 139 ['CIC3', 'CIC4', 'CIC5'] 
     [1;46mCluster[0m 140 ['ATSC6p'] 
     [1;46mCluster[0m 141 ['SsssCH'] 
     [1;46mCluster[0m 142 ['ATSC7c'] 
     [1;46mCluster[0m 143 ['PEOE_VSA8'] 
     [1;46mCluster[0m 144 ['SlogP_VSA3'] 
     [1;46mCluster[0m 145 ['ATSC2Z', 'ATSC2m'] 
     [1;46mCluster[0m 146 ['PEOE_VSA4'] 
     [1;46mCluster[0m 147 ['GATS2v', 'GATS2p'] 
     [1;46mCluster[0m 148 ['AATSC1c'] 
     [1;46mCluster[0m 149 ['n6ARing'] 
     [1;46mCluster[0m 150 ['AATSC1Z', 'AATSC1m'] 
     [1;46mCluster[0m 151 ['ATSC8d'] 
     [1;46mCluster[0m 152 ['ATSC5c'] 
     [1;46mCluster[0m 153 ['SdssS'] 
     [1;46mCluster[0m 154 ['ATSC3dv'] 
     [1;46mCluster[0m 155 ['EState_VSA3'] 
     [1;46mCluster[0m 156 ['ATS1s', 'ATS2s', 'ATS3s', 'ATS4s', 'ATS5s', 'ATSC0dv'] 
     [1;46mCluster[0m 157 ['Xc-3dv'] 
     [1;46mCluster[0m 158 ['NdO', 'SdO', 'VSA_EState2'] 
     [1;46mCluster[0m 159 ['ATSC2s', 'AATSC2s'] 
     [1;46mCluster[0m 160 ['nAromAtom', 'nAromBond', 'nBondsA', 'nBondsM', 'n6Ring', 'naRing', 'n6aRing'] 
     [1;46mCluster[0m 161 ['ATSC2d'] 
     [1;46mCluster[0m 162 ['NdssC'] 
     [1;46mCluster[0m 163 ['JGI2'] 
     [1;46mCluster[0m 164 ['AATSC0v'] 
     [1;46mCluster[0m 165 ['GATS1v', 'GATS1p'] 
     [1;46mCluster[0m 166 ['AATS0se', 'AATS1se', 'AATS2se', 'AATS0pe', 'AATS2pe', 'Mse', 'Mpe', 'ETA_epsilon_1', 'ETA_epsilon_4', 'ETA_dEpsilon_A'] 
     [1;46mCluster[0m 167 ['SIC0', 'BIC0'] 
     [1;46mCluster[0m 168 ['EState_VSA4'] 
     [1;46mCluster[0m 169 ['AATSC0p'] 
     [1;46mCluster[0m 170 ['ATSC3c'] 
     [1;46mCluster[0m 171 ['nBondsT', 'NtN', 'StN', 'SMR_VSA2'] 
     [1;46mCluster[0m 172 ['SdssC'] 
     [1;46mCluster[0m 173 ['AXp-0d', 'ETA_shape_p'] 
     [1;46mCluster[0m 174 ['nG12FaHRing'] 
     [1;46mCluster[0m 175 ['JGI5'] 
     [1;46mCluster[0m 176 ['PEOE_VSA11'] 
     [1;46mCluster[0m 177 ['ATSC8i'] 
     [1;46mCluster[0m 178 ['AATSC0Z', 'AATSC0m'] 
     [1;46mCluster[0m 179 ['MIC0', 'MIC1'] 
     [1;46mCluster[0m 180 ['ATSC8Z', 'ATSC8m'] 
     [1;46mCluster[0m 181 ['GATS2i'] 
     [1;46mCluster[0m 182 ['Xpc-4d', 'GGI2'] 
     [1;46mCluster[0m 183 ['RPCG'] 
     [1;46mCluster[0m 184 ['BCUTs-1l'] 
     [1;46mCluster[0m 185 ['ATSC1v'] 
     [1;46mCluster[0m 186 ['ATSC7p'] 
     [1;46mCluster[0m 187 ['nFARing', 'nFAHRing'] 
     [1;46mCluster[0m 188 ['AATSC0i'] 
     [1;46mCluster[0m 189 ['BCUTd-1l'] 
     [1;46mCluster[0m 190 ['ATSC3d'] 
     [1;46mCluster[0m 191 ['piPC9'] 
     [1;46mCluster[0m 192 ['ETA_dPsi_A'] 
     [1;46mCluster[0m 193 ['SaasC'] 
     [1;46mCluster[0m 194 ['NddsN'] 
     [1;46mCluster[0m 195 ['NsssCH'] 
     [1;46mCluster[0m 196 ['JGI7'] 
     [1;46mCluster[0m 197 ['PEOE_VSA3'] 
     [1;46mCluster[0m 198 ['PEOE_VSA6'] 
     [1;46mCluster[0m 199 ['ATSC0se', 'ATSC0pe', 'ATSC0are'] 
     [1;46mCluster[0m 200 ['ATS8Z', 'ATS8m', 'GGI7'] 
     [1;46mCluster[0m 201 ['SIC2', 'BIC2'] 
     [1;46mCluster[0m 202 ['nBondsKD', 'ETA_beta', 'ETA_beta_ns', 'ETA_eta_F', 'AETA_eta_F', 'ETA_eta_FL'] 
     [1;46mCluster[0m 203 ['NsCH3', 'SsCH3'] 
     [1;46mCluster[0m 204 ['GATS1c'] 
     [1;46mCluster[0m 205 ['AATSC2v', 'AATSC2p', 'MATS2v', 'MATS2p'] 
     [1;46mCluster[0m 206 ['SddssS'] 
     [1;46mCluster[0m 207 ['ETA_shape_x'] 
     [1;46mCluster[0m 208 ['JGI4'] 
     [1;46mCluster[0m 209 ['JGI8'] 
     [1;46mCluster[0m 210 ['MATS2Z', 'MATS2m'] 
     [1;46mCluster[0m 211 ['BalabanJ'] 
     [1;46mCluster[0m 212 ['RNCG'] 
     [1;46mCluster[0m 213 ['PEOE_VSA12'] 
     [1;46mCluster[0m 214 ['SMR_VSA3'] 
     [1;46mCluster[0m 215 ['n5AHRing'] 
     [1;46mCluster[0m 216 ['FilterItLogS'] 
     [1;46mCluster[0m 217 ['EState_VSA2'] 
     [1;46mCluster[0m 218 ['ATSC2dv', 'AATSC2dv', 'MATS2dv'] 
     [1;46mCluster[0m 219 ['VSA_EState9'] 
     [1;46mCluster[0m 220 ['BCUTZ-1h', 'BCUTm-1h'] 
     [1;46mCluster[0m 221 ['Xc-3d', 'ETA_eta_B', 'ETA_eta_BR', 'GGI1', 'mZagreb1'] 
     [1;46mCluster[0m 222 ['n6HRing', 'n6aHRing'] 
     [1;46mCluster[0m 223 ['nFHRing'] 
     [1;46mCluster[0m 224 ['AATSC1d'] 
     [1;46mCluster[0m 225 ['MATS1s'] 
     [1;46mCluster[0m 226 ['ATSC4i'] 
     [1;46mCluster[0m 227 ['SlogP_VSA10'] 
     [1;46mCluster[0m 228 ['nRot'] 
     [1;46mCluster[0m 229 ['ATSC5d'] 
     [1;46mCluster[0m 230 ['ATS0dv', 'ATS1dv', 'ATS2dv', 'ATS3dv', 'ATS4dv', 'ATS5dv', 'GGI5'] 
     [1;46mCluster[0m 231 ['JGI3'] 
     [1;46mCluster[0m 232 ['PEOE_VSA9'] 
     [1;46mCluster[0m 233 ['BCUTi-1h'] 
     [1;46mCluster[0m 234 ['ETA_beta_ns_d', 'AETA_beta_ns_d'] 
     [1;46mCluster[0m 235 ['ATSC4Z', 'ATSC4m'] 
     [1;46mCluster[0m 236 ['PEOE_VSA1'] 
     [1;46mCluster[0m 237 ['PEOE_VSA10'] 
     [1;46mCluster[0m 238 ['nHBAcc', 'TopoPSA'] 
     [1;46mCluster[0m 239 ['AATSC0c'] 
     [1;46mCluster[0m 240 ['ETA_dAlpha_B'] 
     [1;46mCluster[0m 241 ['ATSC1Z', 'ATSC1m'] 
     [1;46mCluster[0m 242 ['AATS0d', 'AATS1d', 'AATS2d'] 
     [1;46mCluster[0m 243 ['SlogP_VSA8'] 
     [1;46mCluster[0m 244 ['ATSC7v'] 
     [1;46mCluster[0m 245 ['nG12FARing', 'nG12FAHRing'] 
     [1;46mCluster[0m 246 ['AMID_X'] 
     [1;46mCluster[0m 247 ['GATS2c'] 
     [1;46mCluster[0m 248 ['ATSC8dv', 'ATSC8se', 'ATSC8pe', 'ATSC8are'] 
     [1;46mCluster[0m 249 ['ATS5Z', 'ATS6Z', 'ATS5m', 'ATS6m'] 
     [1;46mCluster[0m 250 ['EState_VSA1'] 
     [1;46mCluster[0m 251 ['ATSC6c'] 
     [1;46mCluster[0m 252 ['EState_VSA6'] 
     [1;46mCluster[0m 253 ['NssO', 'SssO'] 
     [1;46mCluster[0m 254 ['ATSC4c'] 
     [1;46mCluster[0m 255 ['SdsssP'] 
     [1;46mCluster[0m 256 ['AATS2Z', 'AATS2m'] 
     [1;46mCluster[0m 257 ['ATSC3s'] 
     [1;46mCluster[0m 258 ['VSA_EState7'] 
     [1;46mCluster[0m 259 ['Xch-6d', 'Xch-7d', 'Xch-6dv'] 
     [1;46mCluster[0m 260 ['SddsN'] 
     [1;46mCluster[0m 261 ['NddssS'] 
     [1;46mCluster[0m 262 ['PEOE_VSA7'] 
     [1;46mCluster[0m 263 ['ATS5d', 'ATS5v', 'ATS5se', 'ATS5pe', 'ATS5are', 'ATS5p', 'ATS5i'] 
     [1;46mCluster[0m 264 ['C3SP2'] 
     [1;46mCluster[0m 265 ['JGI6'] 
     [1;46mCluster[0m 266 ['SlogP_VSA1'] 
     [1;46mCluster[0m 267 ['ATSC7i'] 
     [1;46mCluster[0m 268 ['nH', 'nBondsS', 'bpol'] 
     [1;46mCluster[0m 269 ['ATS2Z', 'ATS2m'] 
     [1;46mCluster[0m 270 ['Xpc-4dv', 'Xpc-5dv', 'Xpc-6dv'] 
     [1;46mCluster[0m 271 ['ATSC7d'] 
     [1;46mCluster[0m 272 ['SMR_VSA9'] 
     [1;46mCluster[0m 273 ['ATSC6i'] 
     [1;46mCluster[0m 274 ['n3HRing', 'n3AHRing'] 
     [1;46mCluster[0m 275 ['ATSC6d'] 
     [1;46mCluster[0m 276 ['ATS0Z', 'ATS0m'] 
     [1;46mCluster[0m 277 ['AATS2i', 'GATS1i'] 
     [1;46mCluster[0m 278 ['ATSC1i'] 
     [1;46mCluster[0m 279 ['ATSC3i'] 
     [1;46mCluster[0m 280 ['ETA_epsilon_3', 'AMID', 'nRing'] 
     [1;46mCluster[0m 281 ['AMID_N'] 
     [1;46mCluster[0m 282 ['ATSC4dv'] 
     [1;46mCluster[0m 283 ['GATS2s'] 
     [1;46mCluster[0m 284 ['SMR_VSA4'] 
     [1;46mCluster[0m 285 ['ATSC4s'] 
     [1;46mCluster[0m 286 ['n5HRing'] 
     [1;46mCluster[0m 287 ['BCUTd-1h'] 
     [1;46mCluster[0m 288 ['AETA_eta'] 
     [1;46mCluster[0m 289 ['GATS1s'] 
     [1;46mCluster[0m 290 ['GATS2dv'] 
     [1;46mCluster[0m 291 ['VSA_EState3'] 
     [1;46mCluster[0m 292 ['ETA_shape_y'] 
     [1;46mCluster[0m 293 ['SLogP'] 
     [1;46mCluster[0m 294 ['SpMAD_DzZ', 'SpMAD_Dzm'] 
     [1;46mCluster[0m 295 ['JGI9'] 
     [1;46mCluster[0m 296 ['BCUTv-1h', 'BCUTp-1h'] 
     [1;46mCluster[0m 297 ['MATS1v', 'MATS1p'] 
     [1;46mCluster[0m 298 ['nBondsD'] 
     [1;46mCluster[0m 299 ['AATS0p', 'AATS1p', 'AATS2p', 'Mp'] 
     [1;46mCluster[0m 300 ['SlogP_VSA4'] 
     [1;46mCluster[0m 301 ['IC1', 'IC2'] 
     [1;46mCluster[0m 302 ['MATS2c'] 
     [1;46mCluster[0m 303 ['BCUTc-1h'] 
     [1;46mCluster[0m 304 ['NaasC'] 
     [1;46mCluster[0m 305 ['C1SP1'] 
     [1;46mCluster[0m 306 ['SpMAD_Dzv', 'SpMAD_Dzse', 'SpMAD_Dzpe', 'SpMAD_Dzare', 'SpMAD_Dzi'] 
     [1;46mCluster[0m 307 ['AATSC1se', 'AATSC1pe', 'AATSC1are'] 
     [1;46mCluster[0m 308 ['n6AHRing'] 
     [1;46mCluster[0m 309 ['Xp-2dv', 'Xp-3dv', 'Xp-4dv'] 
     [1;46mCluster[0m 310 ['ATSC1dv'] 
     [1;46mCluster[0m 311 ['VR3_A', 'VE1_DzZ', 'VR2_DzZ', 'VE1_Dzm', 'VR2_Dzm', 'VE1_Dzv', 'VR2_Dzv', 'VE1_Dzse', 'VR2_Dzse', 'VE1_Dzpe', 'VR2_Dzpe', 'VE1_Dzare', 'VR2_Dzare', 'VE1_Dzp', 'VR2_Dzp', 'VE1_Dzi', 'VR2_Dzi', 'VE1_D', 'VR2_D', 'AETA_eta_R', 'TMWC10'] 
     [1;46mCluster[0m 312 ['nAHRing'] 
     [1;46mCluster[0m 313 ['RotRatio'] 
     [1;46mCluster[0m 314 ['AATSC0dv', 'BCUTdv-1h'] 
     [1;46mCluster[0m 315 ['ATSC5i'] 
     [1;46mCluster[0m 316 ['ETA_psi_1'] 
     [1;46mCluster[0m 317 ['NaaCH', 'SaaCH', 'SlogP_VSA6', 'VSA_EState6'] 
     [1;46mCluster[0m 318 ['ATSC7dv'] 
     [1;46mCluster[0m 319 ['AATSC2d'] 
     [1;46mCluster[0m 320 ['ATS7dv', 'ATS8dv', 'GGI8'] 
     [1;46mCluster[0m 321 ['Xc-4d', 'Xc-4dv'] 
     [1;46mCluster[0m 322 ['MATS1c'] 
     [1;46mCluster[0m 323 ['ATSC1c'] 
     [1;46mCluster[0m 324 ['ATSC3p'] 
     [1;46mCluster[0m 325 ['ETA_dEpsilon_B'] 
     [1;46mCluster[0m 326 ['VSA_EState1'] 
     [1;46mCluster[0m 327 ['nFaHRing'] 
     [1;46mCluster[0m 328 ['GATS1dv'] 
     [1;46mCluster[0m 329 ['AATS1s'] 
     [1;46mCluster[0m 330 ['ZMIC1', 'ZMIC2'] 
     [1;46mCluster[0m 331 ['SpMAD_A'] 
     [1;46mCluster[0m 332 ['AETA_beta_s'] 
     [1;46mCluster[0m 333 ['BCUTc-1l'] 
     [1;46mCluster[0m 334 ['AMID_h'] 
     [1;46mCluster[0m 335 ['PEOE_VSA13'] 
     [1;46mCluster[0m 336 ['BCUTse-1h', 'BCUTpe-1h'] 
     [1;46mCluster[0m 337 ['nHBDon'] 
     [1;46mCluster[0m 338 ['C2SP3'] 
     [1;46mCluster[0m 339 ['MATS2se', 'MATS2pe', 'MATS2are'] 
     [1;46mCluster[0m 340 ['fragCpx', 'MPC6', 'MPC7', 'TMPC10'] 
     [1;46mCluster[0m 341 ['ATSC8c'] 
     [1;46mCluster[0m 342 ['nFRing'] 
     [1;46mCluster[0m 343 ['SMR_VSA6'] 
     [1;46mCluster[0m 344 ['TopoPSA(NO)'] 
     [1;46mCluster[0m 345 ['ETA_dEpsilon_C'] 
     [1;46mCluster[0m 346 ['AATSC1s'] 
     [1;46mCluster[0m 347 ['ATS7d', 'ATS8d', 'ATS8s', 'ATS8v', 'ATS8se', 'ATS8pe', 'ATS8are', 'ATS8p', 'ATS8i'] 
     [1;46mCluster[0m 348 ['AATSC0are', 'AMID_O'] 
     [1;46mCluster[0m 349 ['ATS1d', 'ATS2d', 'ATS3d', 'ATS2v', 'ATS3v', 'ATSC0d', 'BertzCT', 'Xp-3d', 'Xp-4d', 'MPC2', 'MPC3', 'TSRW10', 'WPol', 'Zagreb1', 'Zagreb2'] 
     [1;46mCluster[0m 350 ['IC0'] 
     [1;46mCluster[0m 351 ['AATSC0d'] 
     [1;46mCluster[0m 352 ['nC', 'MID_C'] 
     [1;46mCluster[0m 353 ['ATSC3se'] 
     [1;46mCluster[0m 354 ['FCSP3'] 
     [1;46mCluster[0m 355 ['TIC1', 'TIC2', 'TIC3', 'TIC4', 'TIC5'] 
     [1;46mCluster[0m 356 ['AATSC0se', 'AATSC0pe'] 
     [1;46mCluster[0m 357 ['nHRing'] 
     [1;46mCluster[0m 358 ['n5Ring', 'n5ARing'] 
     [1;46mCluster[0m 359 ['SpAbs_Dzv', 'SpMax_Dzv', 'SpDiam_Dzv', 'SpAD_Dzv', 'LogEE_Dzv', 'SpAbs_Dzp', 'SpMax_Dzp', 'SpDiam_Dzp', 'SpAD_Dzp', 'SpMAD_Dzp', 'LogEE_Dzp'] 
     [1;46mCluster[0m 360 ['nG12FaRing'] 
     [1;46mCluster[0m 361 ['AXp-1d'] 
     [1;46mCluster[0m 362 ['ATSC6s'] 
     [1;46mCluster[0m 363 ['nS'] 
     [1;46mCluster[0m 364 ['ATS6se', 'ATS6pe', 'ATS6are', 'ATS6i'] 
     [1;46mCluster[0m 365 ['AATS0Z', 'AATS0m'] 
     [1;46mCluster[0m 366 ['ATSC5dv'] 
     [1;46mCluster[0m 367 ['ETA_dBeta'] 
     [1;46mCluster[0m 368 ['ATSC7s'] 
     [1;46mCluster[0m 369 ['VSA_EState8'] 
     [1;46mCluster[0m 370 ['n5aHRing'] 
     [1;46mCluster[0m 371 ['CIC0'] 
     [1;46mCluster[0m 372 ['AATSC2c'] 
     [1;46mCluster[0m 373 ['ETA_eta', 'ETA_eta_L'] 
     [1;46mCluster[0m 374 ['AETA_eta_L'] 
     [1;46mCluster[0m 375 ['ATS6dv', 'ATS6s'] 
     [1;46mCluster[0m 376 ['JGT10'] 
     [1;46mCluster[0m 377 ['CIC1'] 
     [1;46mCluster[0m 378 ['piPC7', 'piPC8'] 
     [1;46mCluster[0m 379 ['AATSC1v'] 
     [1;46mCluster[0m 380 ['AETA_alpha', 'ETA_dAlpha_A'] 
     [1;46mCluster[0m 381 ['ATSC4v'] 
     [1;46mCluster[0m 382 ['Xc-5d', 'Xpc-5d', 'Xpc-6d', 'GGI3'] 
     [1;46mCluster[0m 383 ['ATS7s', 'ATS7se', 'ATS7pe', 'ATS7are', 'ATS7i'] 
     [1;46mCluster[0m 384 ['JGI10'] 
     [1;46mCluster[0m 385 ['nCl', 'NsCl', 'SsCl', 'EState_VSA9'] 
     [1;46mCluster[0m 386 ['ATSC8p'] 
     [1;46mCluster[0m 387 ['fMF'] 
     [1;46mCluster[0m 388 ['ATS7Z', 'ATS7m', 'GGI6'] 
     [1;46mCluster[0m 389 ['ATSC1d'] 
     [1;46mCluster[0m 390 ['piPC10'] 
     [1;46mCluster[0m 391 ['nHetero', 'MID_h'] 
     [1;46mCluster[0m 392 ['SpAbs_Dzse', 'SpMax_Dzse', 'SpDiam_Dzse', 'SpAD_Dzse', 'LogEE_Dzse', 'SpAbs_Dzpe', 'SpMax_Dzpe', 'SpDiam_Dzpe', 'SpAD_Dzpe', 'LogEE_Dzpe', 'SpAbs_Dzare', 'SpMax_Dzare', 'SpDiam_Dzare', 'SpAD_Dzare', 'LogEE_Dzare', 'SpAbs_Dzi', 'SpMax_Dzi', 'SpDiam_Dzi', 'SpAD_Dzi', 'LogEE_Dzi'] 
     [1;46mCluster[0m 393 ['SpMAD_D', 'Kier1', 'Diameter', 'Radius', 'mZagreb2'] 
     [1;46mCluster[0m 394 ['BCUTs-1h'] 
     [1;46mCluster[0m 395 ['PEOE_VSA2'] 
     [1;46mCluster[0m 396 ['ATSC6v'] 
     [1;46mCluster[0m 397 ['ATSC0i'] 
     [1;46mCluster[0m 398 ['GGI4'] 
     [1;46mCluster[0m 399 ['SM1_Dzare', 'SM1_Dzi'] 
    
    Cluster info file [1mqdb1_F_s_kmeans.cluster[0m file saved



```python
# Evaluating clustered results
clust.cluster_dist()
```


    
![png](output_38_0.png)
    



    interactive(children=(ToggleButtons(description='Bins Index', options=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), value=0)…



```python
clust.model_evaluation()
## Silhouette Score에 계산에 해당하는 Nearest Cluster 확인하는 방법 찾아볼 것. (아직 진행하지 못함)
```

    silhouette average score : 0.584343
    silhouette socre group by 
     Cluster
    0      0.000000
    1      0.794959
    2      0.711150
    3      0.000000
    4      0.604672
             ...   
    395    0.000000
    396    0.000000
    397    0.000000
    398    0.000000
    399    0.609406
    Name: silhouette_coeff, Length: 400, dtype: float64
    
    
    overlapped Cluseter Num :  [0, 3, 18, 22, 30, 41, 44, 51, 53, 54, 70, 80, 81, 87, 88, 90, 96, 98, 99, 105, 112, 116, 118, 120, 122, 123, 127, 129, 130, 131, 132, 133, 134, 136, 140, 141, 142, 143, 144, 146, 148, 149, 151, 152, 153, 154, 155, 157, 161, 162, 163, 164, 168, 169, 170, 172, 174, 175, 176, 177, 181, 183, 184, 185, 186, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 204, 206, 207, 208, 209, 211, 212, 213, 214, 215, 216, 217, 219, 223, 224, 225, 226, 227, 228, 229, 231, 232, 233, 236, 237, 239, 240, 243, 244, 246, 247, 250, 251, 252, 254, 255, 257, 258, 260, 261, 262, 264, 265, 266, 267, 271, 272, 273, 275, 278, 279, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 295, 298, 300, 302, 303, 304, 305, 308, 310, 312, 313, 315, 316, 318, 319, 322, 323, 324, 325, 326, 327, 328, 329, 331, 332, 333, 334, 335, 337, 338, 341, 342, 343, 344, 345, 346, 350, 351, 353, 354, 357, 360, 361, 362, 363, 366, 367, 368, 369, 370, 371, 372, 374, 376, 377, 379, 381, 384, 386, 387, 389, 390, 394, 395, 396, 397, 398]
    wrong/incorrect Cluster Num :  []
    
    
    kmeans.eval file saved



    
![png](output_39_1.png)
    


<h2> 2.2-2. S.O.M (Caution : Takes long time)</h2>
<h3> mt.FeatureCluster_Minisom(train_data, error, random_seed, train,neighborhood_function, topology) : Apply S.O.M Class </h3>

    * data = Train Data
    * Error = (default) qe, te
    * random_seed = (default) 0
    * train = time to init, (default) 1000
    * neighborhood_function = function weight neightborhookd of a position in mpa
        Possible : (default) 'gaussian', 'mexican_hat', 'bubble', 'triangle'
    * topology = (default) 'hexagonal', 'rectangular'
        - Error te can't use hexagonal topology
    
    Parameter_Combine(min_size,max_size) = expect #cluster n x n in range(min_size to max_size)
    
    set_cluster(map_n,sigma,lr,init) = get cluster value as array
        * map_n : map size
        * sigma
        * lr : initial learning_rate
        * init : metohd of initai weight
            Possible : 'pca', 'random'


```python
dt.man('som.png')
```


    
![png](output_42_0.png)
    



```python
# S.O.M의 경우 args에 따라서 cluster의 개수가 다양하기 때문에, 해당 부분을 미리 확인해보는 단계
clust = mt.FeatureCluster_Minisom(X_train)
clust.parameter_combine(18,19)
```

    시작시간: 2026-05-07 10:42:26.405380 
    


    /home/jun/pyqsar/pyqsar3r25/pyqsar/minisom.py:477: RuntimeWarning: invalid value encountered in sqrt
      return sqrt(-2 * cross_term + input_data_sq + weights_flat_sq.T)


    
    종료시간: 2026-05-07 10:43:20.050035 
    총 소요시간: 0 days 00:00:53.644655





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>map_size</th>
      <th>sigma</th>
      <th>learning_rate</th>
      <th>init_method</th>
      <th>qe</th>
      <th>n_cluster</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>18x18</td>
      <td>0.1</td>
      <td>0.9</td>
      <td>random_init</td>
      <td>3.243813</td>
      <td>317</td>
    </tr>
    <tr>
      <th>1</th>
      <td>18x18</td>
      <td>0.2</td>
      <td>0.9</td>
      <td>random_init</td>
      <td>3.243814</td>
      <td>317</td>
    </tr>
    <tr>
      <th>2</th>
      <td>18x18</td>
      <td>0.3</td>
      <td>0.9</td>
      <td>random_init</td>
      <td>3.244328</td>
      <td>315</td>
    </tr>
    <tr>
      <th>3</th>
      <td>18x18</td>
      <td>0.4</td>
      <td>0.9</td>
      <td>random_init</td>
      <td>3.256003</td>
      <td>303</td>
    </tr>
    <tr>
      <th>4</th>
      <td>18x18</td>
      <td>0.1</td>
      <td>0.8</td>
      <td>random_init</td>
      <td>3.266577</td>
      <td>317</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>175</th>
      <td>18x18</td>
      <td>0.3</td>
      <td>0.1</td>
      <td>random_init</td>
      <td>8.901528</td>
      <td>26</td>
    </tr>
    <tr>
      <th>176</th>
      <td>18x18</td>
      <td>0.2</td>
      <td>0.1</td>
      <td>random_init</td>
      <td>8.901602</td>
      <td>26</td>
    </tr>
    <tr>
      <th>177</th>
      <td>18x18</td>
      <td>0.1</td>
      <td>0.1</td>
      <td>random_init</td>
      <td>8.901602</td>
      <td>26</td>
    </tr>
    <tr>
      <th>178</th>
      <td>18x18</td>
      <td>0.4</td>
      <td>0.1</td>
      <td>random_init</td>
      <td>8.918685</td>
      <td>26</td>
    </tr>
    <tr>
      <th>179</th>
      <td>18x18</td>
      <td>0.5</td>
      <td>0.1</td>
      <td>random_init</td>
      <td>8.924816</td>
      <td>27</td>
    </tr>
  </tbody>
</table>
<p>180 rows × 6 columns</p>
</div>




```python
# SOM Clustering 
clust_info = clust.set_cluster(map_n=18,sigma=0.1,lr=0.9,init='random')
```

    /home/jun/pyqsar/pyqsar3r25/pyqsar/minisom.py:477: RuntimeWarning: invalid value encountered in sqrt
      return sqrt(-2 * cross_term + input_data_sq + weights_flat_sq.T)


    Quantization Error : 3.243813384760734
    
     [1;46mCluster[0m 0 ['JGI3'] 
     [1;46mCluster[0m 1 ['Xp-2d', 'Xp-3d', 'Xp-4d', 'MPC2', 'Zagreb1', 'Zagreb2'] 
     [1;46mCluster[0m 2 ['nO', 'ATSC0c', 'SMR_VSA1', 'MID_O', 'AMID_O'] 
     [1;46mCluster[0m 3 ['ATS8s'] 
     [1;46mCluster[0m 4 ['VE3_Dzp', 'VE3_D'] 
     [1;46mCluster[0m 5 ['SpAbs_D', 'SpMax_D', 'SpDiam_D', 'SpAD_D', 'LogEE_D'] 
     [1;46mCluster[0m 6 ['NaaNH'] 
     [1;46mCluster[0m 7 ['MIC5'] 
     [1;46mCluster[0m 8 ['SdssC', 'VSA_EState5'] 
     [1;46mCluster[0m 9 ['n9FRing', 'n9FHRing'] 
     [1;46mCluster[0m 10 ['ATS3se', 'ATS3pe', 'ATS3are'] 
     [1;46mCluster[0m 11 ['GATS1s', 'GATS1se', 'GATS1pe', 'GATS1are'] 
     [1;46mCluster[0m 12 ['AATSC2dv', 'MATS2dv'] 
     [1;46mCluster[0m 13 ['NssNH', 'SssNH'] 
     [1;46mCluster[0m 14 ['ATSC1d'] 
     [1;46mCluster[0m 15 ['ATS4se', 'ATS4pe', 'ATS4are'] 
     [1;46mCluster[0m 16 ['MIC1', 'MIC2'] 
     [1;46mCluster[0m 17 ['ATS2v', 'ATS3v'] 
     [1;46mCluster[0m 18 ['C2SP3', 'NssCH2', 'SssCH2', 'VSA_EState7'] 
     [1;46mCluster[0m 19 ['AETA_eta_B', 'JGI1'] 
     [1;46mCluster[0m 20 ['SpAbs_Dzp', 'SpMax_Dzp', 'SpAD_Dzp', 'LogEE_Dzp'] 
     [1;46mCluster[0m 21 ['BCUTs-1l', 'BCUTi-1l'] 
     [1;46mCluster[0m 22 ['ATSC4c', 'ATSC4dv', 'ATSC4s', 'ATSC4se', 'ATSC4pe', 'ATSC4are'] 
     [1;46mCluster[0m 23 ['nBridgehead', 'C3SP3', 'Xch-7dv', 'Xc-6d', 'Xc-5dv', 'Xc-6dv', 'Xpc-6dv', 'NsssCH', 'NssssC', 'SMR_VSA4', 'nARing', 'n5ARing'] 
     [1;46mCluster[0m 24 ['AATSC0dv', 'BCUTdv-1h', 'BCUTs-1h', 'BCUTse-1h', 'BCUTpe-1h', 'BCUTare-1h'] 
     [1;46mCluster[0m 25 ['AATSC0se'] 
     [1;46mCluster[0m 26 ['ATSC1se', 'ATSC1pe', 'ATSC1are'] 
     [1;46mCluster[0m 27 ['EState_VSA7'] 
     [1;46mCluster[0m 28 ['nAtom', 'ATS0se', 'ATS1se', 'ATS2se', 'ATS0pe', 'ATS1pe', 'ATS2pe', 'ATS0are', 'ATS1are', 'ATS2are', 'ATS0i', 'ATS1i', 'ATS2i', 'ATSC0v', 'nBonds', 'nBondsS', 'nBondsKS', 'Sse', 'Spe', 'Sare', 'Si', 'bpol'] 
     [1;46mCluster[0m 29 ['AATSC1Z', 'AATSC1m'] 
     [1;46mCluster[0m 30 ['VR3_A'] 
     [1;46mCluster[0m 31 ['PEOE_VSA5', 'AMID_X'] 
     [1;46mCluster[0m 32 ['SsssCH'] 
     [1;46mCluster[0m 33 ['AATS0s', 'AATS1s', 'AATS2s', 'AATSC0s'] 
     [1;46mCluster[0m 34 ['nHBAcc', 'TopoPSA'] 
     [1;46mCluster[0m 35 ['PEOE_VSA6', 'SLogP'] 
     [1;46mCluster[0m 36 ['SIC1', 'SIC2', 'BIC1', 'BIC2'] 
     [1;46mCluster[0m 37 ['SpAbs_Dzi', 'SpMax_Dzi', 'SpAD_Dzi', 'LogEE_Dzi'] 
     [1;46mCluster[0m 38 ['ETA_epsilon_3', 'AMID', 'piPC3', 'piPC4', 'piPC5', 'piPC6', 'TpiPC10'] 
     [1;46mCluster[0m 39 ['AATSC0Z', 'AATSC0m'] 
     [1;46mCluster[0m 40 ['ATS7s', 'ATS7se', 'ATS7pe', 'ATS7are', 'ATS7i'] 
     [1;46mCluster[0m 41 ['VE1_Dzse', 'VE1_Dzpe', 'VE1_Dzare'] 
     [1;46mCluster[0m 42 ['VR2_Dzm'] 
     [1;46mCluster[0m 43 ['Xp-5dv', 'Xp-6dv', 'Xp-7dv'] 
     [1;46mCluster[0m 44 ['VE1_Dzi'] 
     [1;46mCluster[0m 45 ['PEOE_VSA3'] 
     [1;46mCluster[0m 46 ['SaaNH'] 
     [1;46mCluster[0m 47 ['ATS1s', 'ATSC0dv'] 
     [1;46mCluster[0m 48 ['VE1_D'] 
     [1;46mCluster[0m 49 ['AATS0i', 'AATS1i', 'Mi'] 
     [1;46mCluster[0m 50 ['SpDiam_Dzp'] 
     [1;46mCluster[0m 51 ['ATSC7d', 'ATSC8d'] 
     [1;46mCluster[0m 52 ['AATS0v', 'AATS1v', 'AATS2v', 'Mv'] 
     [1;46mCluster[0m 53 ['VE1_A'] 
     [1;46mCluster[0m 54 ['MATS1pe', 'MATS1are'] 
     [1;46mCluster[0m 55 ['ATSC8Z', 'ATSC8m'] 
     [1;46mCluster[0m 56 ['LogEE_A'] 
     [1;46mCluster[0m 57 ['AATSC1pe', 'AATSC1are'] 
     [1;46mCluster[0m 58 ['ATSC6dv', 'ATSC5d', 'ATSC6d'] 
     [1;46mCluster[0m 59 ['NsOH', 'SsOH'] 
     [1;46mCluster[0m 60 ['nRot', 'RotRatio'] 
     [1;46mCluster[0m 61 ['AATS2se', 'AATS2pe'] 
     [1;46mCluster[0m 62 ['ATSC0se', 'ATSC0pe', 'ATSC0are', 'VSA_EState1'] 
     [1;46mCluster[0m 63 ['nHeavyAtom', 'Xp-1d', 'ETA_eta_RL', 'mZagreb2'] 
     [1;46mCluster[0m 64 ['GATS2c'] 
     [1;46mCluster[0m 65 ['Xc-4d', 'Xc-4dv'] 
     [1;46mCluster[0m 66 ['EState_VSA3'] 
     [1;46mCluster[0m 67 ['ATS0Z', 'ATS0m', 'ATSC0Z', 'ATSC0m', 'BCUTZ-1h', 'BCUTm-1h'] 
     [1;46mCluster[0m 68 ['ETA_epsilon_4'] 
     [1;46mCluster[0m 69 ['ATSC4Z', 'ATSC4m', 'ATSC4v', 'ATSC4p', 'ATSC4i'] 
     [1;46mCluster[0m 70 ['NaaO', 'SaaO'] 
     [1;46mCluster[0m 71 ['SpMAD_A'] 
     [1;46mCluster[0m 72 ['nAromAtom', 'nAromBond', 'nBondsA', 'n6Ring', 'naRing', 'n6aRing'] 
     [1;46mCluster[0m 73 ['MATS2i'] 
     [1;46mCluster[0m 74 ['SpDiam_DzZ', 'SpDiam_Dzm'] 
     [1;46mCluster[0m 75 ['ATSC6se', 'ATSC6pe', 'ATSC6are'] 
     [1;46mCluster[0m 76 ['NdssC', 'ETA_eta', 'ETA_eta_L'] 
     [1;46mCluster[0m 77 ['nBase', 'C1SP2'] 
     [1;46mCluster[0m 78 ['n6AHRing'] 
     [1;46mCluster[0m 79 ['VR1_Dzv', 'VR1_Dzp'] 
     [1;46mCluster[0m 80 ['ATS4d', 'ATS4v'] 
     [1;46mCluster[0m 81 ['TSRW10'] 
     [1;46mCluster[0m 82 ['Xp-0d', 'ETA_beta_s', 'ETA_eta_R', 'Kier1'] 
     [1;46mCluster[0m 83 ['AATS0Z', 'AATS2Z', 'AATS0m', 'AATS2m', 'AETA_alpha', 'ETA_dAlpha_A'] 
     [1;46mCluster[0m 84 ['MZ', 'Mm', 'AMW'] 
     [1;46mCluster[0m 85 ['AATS1Z', 'AATS1m'] 
     [1;46mCluster[0m 86 ['CIC2', 'CIC3', 'CIC4', 'CIC5'] 
     [1;46mCluster[0m 87 ['n9FaRing', 'n9FaHRing'] 
     [1;46mCluster[0m 88 ['ATSC2i', 'AATSC2i'] 
     [1;46mCluster[0m 89 ['ATS2Z', 'ATS2m', 'Xp-2dv'] 
     [1;46mCluster[0m 90 ['JGI2'] 
     [1;46mCluster[0m 91 ['ATS3i'] 
     [1;46mCluster[0m 92 ['MIC4'] 
     [1;46mCluster[0m 93 ['VE3_Dzv', 'VE3_Dzi'] 
     [1;46mCluster[0m 94 ['AATSC0v', 'AETA_eta_RL', 'ETA_dEpsilon_C', 'AMID_C'] 
     [1;46mCluster[0m 95 ['SpMAD_DzZ', 'SpMAD_Dzm'] 
     [1;46mCluster[0m 96 ['NdCH2', 'SdCH2'] 
     [1;46mCluster[0m 97 ['AATS0p', 'AATS1p', 'AATS2p', 'Mp'] 
     [1;46mCluster[0m 98 ['fragCpx'] 
     [1;46mCluster[0m 99 ['ATS2s', 'ATS3s', 'ATS4s', 'ATS5s'] 
     [1;46mCluster[0m 100 ['AATS2i', 'GATS1i', 'BCUTZ-1l', 'BCUTm-1l', 'BCUTse-1l', 'FCSP3', 'AXp-0d'] 
     [1;46mCluster[0m 101 ['AATSC1d', 'MATS1Z', 'MATS1m'] 
     [1;46mCluster[0m 102 ['piPC1'] 
     [1;46mCluster[0m 103 ['VR2_DzZ'] 
     [1;46mCluster[0m 104 ['MATS2s'] 
     [1;46mCluster[0m 105 ['ATS5d', 'ATS6d', 'ATS5v', 'ATS6v', 'ATS5p', 'ATS6p', 'GGI6'] 
     [1;46mCluster[0m 106 ['ATSC1s', 'ATSC5s', 'AATSC1s', 'SddsN', 'SddssS'] 
     [1;46mCluster[0m 107 ['nBondsT', 'NtN', 'StN', 'SMR_VSA2'] 
     [1;46mCluster[0m 108 ['EState_VSA8'] 
     [1;46mCluster[0m 109 ['Xc-5d', 'Xpc-4d', 'Xpc-5d', 'Xpc-6d', 'GGI2', 'GGI3'] 
     [1;46mCluster[0m 110 ['LogEE_Dzm'] 
     [1;46mCluster[0m 111 ['VR2_Dzse', 'VR2_Dzpe', 'VR2_Dzare', 'VR2_Dzi', 'VR2_D', 'TMWC10'] 
     [1;46mCluster[0m 112 ['MPC3', 'GGI4', 'WPol'] 
     [1;46mCluster[0m 113 ['GATS1c', 'PEOE_VSA1'] 
     [1;46mCluster[0m 114 ['nG12FRing', 'nG12FHRing'] 
     [1;46mCluster[0m 115 ['SIC3', 'SIC4', 'SIC5'] 
     [1;46mCluster[0m 116 ['TIC0', 'TIC1'] 
     [1;46mCluster[0m 117 ['VAdjMat', 'SRW02'] 
     [1;46mCluster[0m 118 ['WPath'] 
     [1;46mCluster[0m 119 ['IC1', 'IC2'] 
     [1;46mCluster[0m 120 ['NaaS', 'SaaS'] 
     [1;46mCluster[0m 121 ['Xp-5d', 'Xp-6d', 'MPC4', 'MPC5'] 
     [1;46mCluster[0m 122 ['NdS', 'SdS'] 
     [1;46mCluster[0m 123 ['nN', 'BCUTi-1h', 'SMR_VSA3', 'MID_N', 'AMID_N'] 
     [1;46mCluster[0m 124 ['ATSC7v'] 
     [1;46mCluster[0m 125 ['ATSC1c', 'AATSC1c', 'MATS1c'] 
     [1;46mCluster[0m 126 ['AETA_eta_BR', 'JGT10'] 
     [1;46mCluster[0m 127 ['SpAbs_Dzv', 'SpAD_Dzv'] 
     [1;46mCluster[0m 128 ['EState_VSA1'] 
     [1;46mCluster[0m 129 ['nHRing', 'n6HRing', 'naHRing', 'n6aHRing'] 
     [1;46mCluster[0m 130 ['SIC0', 'BIC0'] 
     [1;46mCluster[0m 131 ['SpMax_Dzv', 'LogEE_Dzv'] 
     [1;46mCluster[0m 132 ['ATSC3d', 'ATSC4d', 'ATSC5Z', 'ATSC5m', 'ATSC5v', 'ATSC5p', 'ATSC5i'] 
     [1;46mCluster[0m 133 ['Xp-7d'] 
     [1;46mCluster[0m 134 ['VR1_A', 'VR2_A'] 
     [1;46mCluster[0m 135 ['IC0'] 
     [1;46mCluster[0m 136 ['SpMax_A', 'SpDiam_A', 'MWC10', 'SRW10'] 
     [1;46mCluster[0m 137 ['AATS0d', 'AATS1d', 'AATS2d', 'ETA_shape_y', 'ETA_dEpsilon_B'] 
     [1;46mCluster[0m 138 ['SpDiam_Dzse', 'SpDiam_Dzpe', 'SpDiam_Dzare', 'SpDiam_Dzi'] 
     [1;46mCluster[0m 139 ['ATS0dv', 'ATS1dv', 'ATS2dv', 'ETA_eta_F'] 
     [1;46mCluster[0m 140 ['AMID_h'] 
     [1;46mCluster[0m 141 ['SpAbs_Dzse', 'SpMax_Dzse', 'SpAD_Dzse', 'LogEE_Dzse'] 
     [1;46mCluster[0m 142 ['NsNH2', 'SsNH2', 'VSA_EState4'] 
     [1;46mCluster[0m 143 ['AATS2are'] 
     [1;46mCluster[0m 144 ['EState_VSA2'] 
     [1;46mCluster[0m 145 ['ATSC8i', 'SssssC'] 
     [1;46mCluster[0m 146 ['ATS1Z', 'ATS1m', 'ATS0p', 'ATSC0p', 'Xp-1dv'] 
     [1;46mCluster[0m 147 ['AXp-0dv', 'AXp-1dv', 'AETA_eta_L'] 
     [1;46mCluster[0m 148 ['SpAbs_DzZ', 'SpAbs_Dzm'] 
     [1;46mCluster[0m 149 ['AATSC0p', 'BCUTv-1h'] 
     [1;46mCluster[0m 150 ['ATSC7i'] 
     [1;46mCluster[0m 151 ['piPC8', 'piPC9', 'piPC10'] 
     [1;46mCluster[0m 152 ['ATS4i'] 
     [1;46mCluster[0m 153 ['JGI6'] 
     [1;46mCluster[0m 154 ['VR2_Dzv', 'VR2_Dzp'] 
     [1;46mCluster[0m 155 ['ETA_dAlpha_B'] 
     [1;46mCluster[0m 156 ['NaaN', 'SaaN', 'PEOE_VSA12'] 
     [1;46mCluster[0m 157 ['ATS7Z', 'ATS7m', 'JGI7', 'JGI8'] 
     [1;46mCluster[0m 158 ['MIC3'] 
     [1;46mCluster[0m 159 ['C1SP3', 'SMR_VSA5', 'SlogP_VSA5', 'EState_VSA4', 'VSA_EState8'] 
     [1;46mCluster[0m 160 ['AETA_eta_FL'] 
     [1;46mCluster[0m 161 ['ATSC6c', 'ATSC6s'] 
     [1;46mCluster[0m 162 ['VR3_DzZ', 'VR3_Dzm', 'VR3_Dzse', 'VR3_Dzpe', 'VR3_Dzare', 'VR3_Dzi', 'VR3_D'] 
     [1;46mCluster[0m 163 ['BCUTd-1l', 'NaaaC', 'SaaaC', 'nFRing', 'nFHRing', 'nFaRing', 'nFaHRing'] 
     [1;46mCluster[0m 164 ['TIC3', 'TIC4', 'TIC5'] 
     [1;46mCluster[0m 165 ['SpMAD_Dzi'] 
     [1;46mCluster[0m 166 ['AATS0dv', 'AATS1dv', 'AATS2dv', 'BCUTc-1h'] 
     [1;46mCluster[0m 167 ['SdssS'] 
     [1;46mCluster[0m 168 ['BCUTpe-1l', 'BCUTare-1l', 'SdsssP'] 
     [1;46mCluster[0m 169 ['SpAD_DzZ', 'SpAD_Dzm'] 
     [1;46mCluster[0m 170 ['NaasN', 'SaasN'] 
     [1;46mCluster[0m 171 ['AATSC0d'] 
     [1;46mCluster[0m 172 ['SpMAD_Dzse', 'SpMAD_Dzpe', 'SpMAD_Dzare'] 
     [1;46mCluster[0m 173 ['NddC', 'NdsN', 'SddC', 'SdsN'] 
     [1;46mCluster[0m 174 ['VE3_DzZ', 'VE3_Dzm'] 
     [1;46mCluster[0m 175 ['nCl', 'nX', 'NsCl', 'SsCl', 'EState_VSA9', 'MID_X'] 
     [1;46mCluster[0m 176 ['VE3_Dzse', 'VE3_Dzpe', 'VE3_Dzare'] 
     [1;46mCluster[0m 177 ['VE1_Dzv', 'VE1_Dzp'] 
     [1;46mCluster[0m 178 ['GATS2se', 'GATS2pe', 'GATS2are'] 
     [1;46mCluster[0m 179 ['SaasC'] 
     [1;46mCluster[0m 180 ['GATS1dv', 'BCUTc-1l', 'BCUTv-1l', 'BCUTp-1l', 'SM1_Dzv', 'SM1_Dzp'] 
     [1;46mCluster[0m 181 ['AATS0are', 'AATS1are', 'Mare'] 
     [1;46mCluster[0m 182 ['ATSC2se', 'ATSC2pe', 'ATSC2are'] 
     [1;46mCluster[0m 183 ['AATSC0c', 'AATSC0pe', 'AATSC0are', 'ETA_dPsi_A'] 
     [1;46mCluster[0m 184 ['MATS2se', 'MATS2pe', 'MATS2are'] 
     [1;46mCluster[0m 185 ['MATS1v', 'MATS1i', 'GATS2i'] 
     [1;46mCluster[0m 186 ['ATSC8c', 'ATSC8dv', 'ATSC8s', 'ATSC8se', 'ATSC8pe', 'ATSC8are'] 
     [1;46mCluster[0m 187 ['NtsC', 'StsC'] 
     [1;46mCluster[0m 188 ['VE3_A'] 
     [1;46mCluster[0m 189 ['ATSC8v', 'ATSC8p'] 
     [1;46mCluster[0m 190 ['ATSC0d', 'Xc-3d', 'ETA_eta_B', 'ETA_eta_BR', 'GGI1', 'mZagreb1'] 
     [1;46mCluster[0m 191 ['ATSC2Z', 'ATSC2m', 'AATSC2Z', 'AATSC2m', 'NsCH3', 'SsCH3', 'PEOE_VSA8'] 
     [1;46mCluster[0m 192 ['ATSC6Z', 'ATSC6m', 'ATSC6v', 'ATSC6p', 'ATSC6i'] 
     [1;46mCluster[0m 193 ['NdsCH', 'SdsCH'] 
     [1;46mCluster[0m 194 ['VE2_A', 'VE2_DzZ', 'VE2_Dzm', 'VE2_Dzv', 'VE2_Dzse', 'VE2_Dzpe', 'VE2_Dzare', 'VE2_Dzp', 'VE2_Dzi', 'RNCG', 'AXp-1d', 'VE2_D', 'FilterItLogS'] 
     [1;46mCluster[0m 195 ['AATS1se', 'AATS1pe'] 
     [1;46mCluster[0m 196 ['Xch-6d', 'Xch-7d', 'Xch-6dv'] 
     [1;46mCluster[0m 197 ['BCUTd-1h'] 
     [1;46mCluster[0m 198 ['nG12FaRing', 'nG12FaHRing'] 
     [1;46mCluster[0m 199 ['IC3', 'IC4', 'IC5'] 
     [1;46mCluster[0m 200 ['ATS5se', 'ATS6se', 'ATS5pe', 'ATS6pe', 'ATS5are', 'ATS6are', 'ATS5i', 'ATS6i'] 
     [1;46mCluster[0m 201 ['NssO', 'SssO', 'SlogP_VSA3'] 
     [1;46mCluster[0m 202 ['SlogP_VSA1'] 
     [1;46mCluster[0m 203 ['ATS3p', 'ATS4p'] 
     [1;46mCluster[0m 204 ['MWC03', 'MWC04', 'SRW04', 'SRW06'] 
     [1;46mCluster[0m 205 ['SpDiam_Dzv'] 
     [1;46mCluster[0m 206 ['ATS1d', 'ATS2d', 'ATS3d'] 
     [1;46mCluster[0m 207 ['ETA_epsilon_2', 'ETA_epsilon_5'] 
     [1;46mCluster[0m 208 ['Xch-3d', 'Xch-3dv', 'n3Ring', 'n3HRing', 'n3ARing', 'n3AHRing', 'SRW03'] 
     [1;46mCluster[0m 209 ['BertzCT', 'ETA_beta', 'piPC7'] 
     [1;46mCluster[0m 210 ['SpMax_DzZ', 'LogEE_DzZ', 'SpMax_Dzm'] 
     [1;46mCluster[0m 211 ['SpAbs_Dzpe'] 
     [1;46mCluster[0m 212 ['SpMax_Dzpe', 'LogEE_Dzpe'] 
     [1;46mCluster[0m 213 ['nC', 'MID_C'] 
     [1;46mCluster[0m 214 ['ATSC1dv', 'AATSC1dv', 'MATS1dv'] 
     [1;46mCluster[0m 215 ['VR3_Dzv', 'VR3_Dzp'] 
     [1;46mCluster[0m 216 ['n5Ring', 'n5HRing', 'SRW05', 'SRW07', 'SRW09'] 
     [1;46mCluster[0m 217 ['SpAbs_A', 'SpAD_A'] 
     [1;46mCluster[0m 218 ['BCUTdv-1l'] 
     [1;46mCluster[0m 219 ['n6ARing'] 
     [1;46mCluster[0m 220 ['nBr', 'NsBr', 'SsBr', 'ETA_psi_1', 'ETA_dPsi_B'] 
     [1;46mCluster[0m 221 ['ATS1p', 'ATS2p'] 
     [1;46mCluster[0m 222 ['n10FRing', 'n10FHRing', 'n10FaRing', 'n10FaHRing'] 
     [1;46mCluster[0m 223 ['piPC2', 'MWC05', 'MWC06', 'MWC07', 'MWC08', 'MWC09', 'SRW08'] 
     [1;46mCluster[0m 224 ['MATS1se'] 
     [1;46mCluster[0m 225 ['NddssS', 'n5AHRing', 'n9FARing', 'n9FAHRing'] 
     [1;46mCluster[0m 226 ['EState_VSA5'] 
     [1;46mCluster[0m 227 ['SsF'] 
     [1;46mCluster[0m 228 ['nF', 'NsF'] 
     [1;46mCluster[0m 229 ['ATS5Z', 'ATS6Z', 'ATS5m', 'ATS6m'] 
     [1;46mCluster[0m 230 ['SpAD_Dzpe'] 
     [1;46mCluster[0m 231 ['ATS7dv', 'ATS8dv', 'GGI8', 'GGI9', 'GGI10', 'JGI10'] 
     [1;46mCluster[0m 232 ['GATS1Z', 'GATS1m', 'GATS1v', 'GATS1p'] 
     [1;46mCluster[0m 233 ['nBondsD', 'NdO', 'SdO', 'VSA_EState2'] 
     [1;46mCluster[0m 234 ['ATSC7Z', 'ATSC7m'] 
     [1;46mCluster[0m 235 ['ATSC1v', 'ATSC1i', 'AATSC1v', 'AATSC1i'] 
     [1;46mCluster[0m 236 ['ATSC7p'] 
     [1;46mCluster[0m 237 ['ATS3dv', 'ATS4dv'] 
     [1;46mCluster[0m 238 ['C4SP3', 'Xch-4d', 'Xch-5d', 'Xch-4dv', 'Xch-5dv'] 
     [1;46mCluster[0m 239 ['AETA_eta_R'] 
     [1;46mCluster[0m 240 ['ATSC2c', 'AATSC2c', 'MATS2c'] 
     [1;46mCluster[0m 241 ['ATS5dv', 'ATS6dv', 'ATS6s', 'GGI5'] 
     [1;46mCluster[0m 242 ['fMF', 'nRing'] 
     [1;46mCluster[0m 243 ['C1SP1'] 
     [1;46mCluster[0m 244 ['ATSC7c'] 
     [1;46mCluster[0m 245 ['SMR_VSA9'] 
     [1;46mCluster[0m 246 ['nI', 'C2SP1', 'NsI', 'SsI'] 
     [1;46mCluster[0m 247 ['TopoShapeIndex', 'PetitjeanIndex'] 
     [1;46mCluster[0m 248 ['nBondsKD', 'ETA_beta_ns', 'AETA_eta_F', 'ETA_eta_FL'] 
     [1;46mCluster[0m 249 ['GATS2dv'] 
     [1;46mCluster[0m 250 ['AATSC0i'] 
     [1;46mCluster[0m 251 ['MATS1s'] 
     [1;46mCluster[0m 252 ['ATSC0i'] 
     [1;46mCluster[0m 253 ['AETA_eta'] 
     [1;46mCluster[0m 254 ['nHBDon'] 
     [1;46mCluster[0m 255 ['SM1_Dzare', 'SM1_Dzi', 'SlogP_VSA10', 'TopoPSA(NO)'] 
     [1;46mCluster[0m 256 ['SpAbs_Dzare', 'SpMax_Dzare', 'SpAD_Dzare', 'LogEE_Dzare'] 
     [1;46mCluster[0m 257 ['NsssN', 'SsssN', 'VSA_EState9'] 
     [1;46mCluster[0m 258 ['PEOE_VSA7'] 
     [1;46mCluster[0m 259 ['MATS2Z', 'MATS2m'] 
     [1;46mCluster[0m 260 ['AETA_beta_s'] 
     [1;46mCluster[0m 261 ['JGI9'] 
     [1;46mCluster[0m 262 ['ATSC3c'] 
     [1;46mCluster[0m 263 ['ATSC3s', 'ATSC3se', 'ATSC3pe', 'ATSC3are'] 
     [1;46mCluster[0m 264 ['nH', 'CIC0', 'CIC1'] 
     [1;46mCluster[0m 265 ['AATSC1se'] 
     [1;46mCluster[0m 266 ['ATSC2d', 'AATSC2d'] 
     [1;46mCluster[0m 267 ['JGI5'] 
     [1;46mCluster[0m 268 ['n5aRing', 'n5aHRing'] 
     [1;46mCluster[0m 269 ['nBondsM', 'C2SP2', 'NaaCH', 'SaaCH', 'SMR_VSA7', 'SlogP_VSA6', 'VSA_EState6'] 
     [1;46mCluster[0m 270 ['NaasC', 'ETA_beta_ns_d', 'AETA_beta_ns_d', 'SlogP_VSA7'] 
     [1;46mCluster[0m 271 ['ATS8Z', 'ATS8m', 'EState_VSA6', 'GGI7'] 
     [1;46mCluster[0m 272 ['ATS3Z', 'ATS4Z', 'ATS3m', 'ATS4m', 'Sm', 'ZMIC0', 'ZMIC1', 'MW'] 
     [1;46mCluster[0m 273 ['Xpc-4dv', 'Xpc-5dv', 'Xp-3dv', 'Xp-4dv'] 
     [1;46mCluster[0m 274 ['ATS0s', 'ATSC0s', 'ATSC2s', 'AATSC2s', 'EState_VSA10'] 
     [1;46mCluster[0m 275 ['RPCG', 'ETA_shape_p'] 
     [1;46mCluster[0m 276 ['NddsN', 'PEOE_VSA2', 'PEOE_VSA13', 'SlogP_VSA4', 'VSA_EState3'] 
     [1;46mCluster[0m 277 ['MIC0'] 
     [1;46mCluster[0m 278 ['SpMAD_Dzv', 'SpMAD_Dzp'] 
     [1;46mCluster[0m 279 ['ATS7d', 'ATS8d', 'ATS7v', 'ATS8v', 'ATS8se', 'ATS8pe', 'ATS8are', 'ATS7p', 'ATS8p', 'ATS8i'] 
     [1;46mCluster[0m 280 ['nS'] 
     [1;46mCluster[0m 281 ['PEOE_VSA10', 'SlogP_VSA11'] 
     [1;46mCluster[0m 282 ['BalabanJ', 'SMR_VSA6', 'SlogP_VSA2'] 
     [1;46mCluster[0m 283 ['ATS0v', 'Xp-0dv', 'SZ', 'Sv', 'Sp', 'ETA_alpha', 'VMcGowan', 'LabuteASA', 'apol', 'SMR'] 
     [1;46mCluster[0m 284 ['SpMAD_D', 'Diameter', 'Radius'] 
     [1;46mCluster[0m 285 ['ATSC5c', 'ATSC5dv', 'ATSC5se', 'ATSC5pe', 'ATSC5are'] 
     [1;46mCluster[0m 286 ['PEOE_VSA11', 'nAHRing', 'nFARing', 'nG12FARing', 'nFAHRing', 'nG12FAHRing'] 
     [1;46mCluster[0m 287 ['BCUTp-1h'] 
     [1;46mCluster[0m 288 ['ECIndex'] 
     [1;46mCluster[0m 289 ['ATSC2v', 'ATSC2p', 'AATSC2v', 'AATSC2p', 'MATS2v', 'MATS2p'] 
     [1;46mCluster[0m 290 ['VE1_DzZ', 'VE1_Dzm'] 
     [1;46mCluster[0m 291 ['ATSC7dv', 'ATSC7s', 'ATSC7se', 'ATSC7pe', 'ATSC7are'] 
     [1;46mCluster[0m 292 ['MPC6', 'MPC7', 'MPC8', 'MPC9', 'MPC10', 'TMPC10'] 
     [1;46mCluster[0m 293 ['ATSC1p', 'AATSC1p', 'MATS1p'] 
     [1;46mCluster[0m 294 ['ATSC2dv', 'ATSC3dv'] 
     [1;46mCluster[0m 295 ['nHetero', 'SM1_Dzse', 'SM1_Dzpe', 'PEOE_VSA9', 'MID_h'] 
     [1;46mCluster[0m 296 ['nP', 'ATSC1Z', 'ATSC1m', 'NdsssP'] 
     [1;46mCluster[0m 297 ['NdssS', 'n7Ring', 'n7HRing', 'n7ARing', 'n7AHRing', 'n12FRing', 'n12FHRing', 'n12FARing', 'n12FAHRing'] 
     [1;46mCluster[0m 298 ['BIC3', 'BIC4', 'BIC5'] 
     [1;46mCluster[0m 299 ['JGI4'] 
     [1;46mCluster[0m 300 ['TIC2'] 
     [1;46mCluster[0m 301 ['MWC02'] 
     [1;46mCluster[0m 302 ['ATSC3Z', 'ATSC3m', 'ATSC3v', 'ATSC3p', 'ATSC3i'] 
     [1;46mCluster[0m 303 ['SM1_DzZ', 'SM1_Dzm', 'ZMIC2', 'ZMIC3', 'ZMIC4', 'ZMIC5'] 
     [1;46mCluster[0m 304 ['AATSC2se', 'AATSC2pe', 'AATSC2are'] 
     [1;46mCluster[0m 305 ['ATS0d', 'ATS1v', 'nBondsO', 'MID', 'MWC01'] 
     [1;46mCluster[0m 306 ['NsSH', 'SsSH', 'PEOE_VSA4'] 
     [1;46mCluster[0m 307 ['AETA_beta', 'AETA_beta_ns', 'ETA_dBeta', 'AETA_dBeta'] 
     [1;46mCluster[0m 308 ['ETA_shape_x'] 
     [1;46mCluster[0m 309 ['VR1_DzZ', 'VR1_Dzm', 'VR1_Dzse', 'VR1_Dzpe', 'VR1_Dzare', 'VR1_Dzi', 'VR1_D'] 
     [1;46mCluster[0m 310 ['AATS0se', 'AATS0pe', 'Mse', 'Mpe', 'ETA_epsilon_1', 'ETA_dEpsilon_A'] 
     [1;46mCluster[0m 311 ['GATS2s'] 
     [1;46mCluster[0m 312 ['C3SP2', 'SlogP_VSA8'] 
     [1;46mCluster[0m 313 ['GATS2Z', 'GATS2m', 'GATS2v', 'GATS2p'] 
     [1;46mCluster[0m 314 ['Xc-3dv'] 
     [1;46mCluster[0m 315 ['ETA_dEpsilon_D'] 
     [1;46mCluster[0m 316 ['NssS', 'SssS'] 
    Cluster info file [1mqdb1_F_s_som.cluster[0m file saved



    
![png](output_44_2.png)
    



```python
# numbef of clusters within the distance range 
clust.cluster_dist()
```


    
![png](output_45_0.png)
    



    interactive(children=(ToggleButtons(description='Bins Index', options=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), value=0)…


<h2> 2.2-3. Hierarchical </h2>
<h3> mt.FeatureCluster(X_data, method, depth) : Apply S.O.M Class </h3>

    *method = linkage method
        Possible values = 'average','complete','single'
    *depth = Standard cut branch
    
    cluster_dist = clust evaluation by histogram
    set_cluster() : get cluster value as array
    


```python
dt.man('hierarchical.png')
```


    
![png](output_47_0.png)
    



```python
# Calculating Cophenet values, chose the largest one
# Not working with BLOCK split 
clust = mt.FeatureCluster(X_train)
clust.cophenetic()
```

    average linkage cophenet: 0.8213907973584261
    complete linkage cophenet: 0.8070560318312705
    single linkage cophenet: 0.35883754064130347



```python
# Hierarchical Clustering 
clust_info = clust.set_cluster('average',1)
```

    
     [1;46mCluster[0m 0 ['ATSC0pe'] 
     [1;46mCluster[0m 1 ['ATSC0are'] 
     [1;46mCluster[0m 2 ['ATSC0se'] 
     [1;46mCluster[0m 3 ['ATS1s'] 
     [1;46mCluster[0m 4 ['ATS3s'] 
     [1;46mCluster[0m 5 ['ATS5s'] 
     [1;46mCluster[0m 6 ['ATSC0dv'] 
     [1;46mCluster[0m 7 ['ATS4dv'] 
     [1;46mCluster[0m 8 ['ATS5dv'] 
     [1;46mCluster[0m 9 ['ATS2s', 'ATS4s'] 
     [1;46mCluster[0m 10 ['ATS0dv'] 
     [1;46mCluster[0m 11 ['ATS6dv'] 
     [1;46mCluster[0m 12 ['GGI5'] 
     [1;46mCluster[0m 13 ['nC', 'MID_C'] 
     [1;46mCluster[0m 14 ['Xp-3d', 'MPC2', 'WPol', 'Zagreb2'] 
     [1;46mCluster[0m 15 ['Xp-4d'] 
     [1;46mCluster[0m 16 ['MPC3'] 
     [1;46mCluster[0m 17 ['TSRW10'] 
     [1;46mCluster[0m 18 ['LogEE_A', 'VE3_A', 'VE3_DzZ', 'VR2_DzZ', 'VR3_DzZ', 'VE3_Dzm', 'VR2_Dzm', 'VR3_Dzm', 'VE3_Dzv', 'VR2_Dzv', 'VR3_Dzv', 'VE3_Dzse', 'VR2_Dzse', 'VR3_Dzse', 'VE3_Dzpe', 'VR2_Dzpe', 'VR3_Dzpe', 'VE3_Dzare', 'VR2_Dzare', 'VR3_Dzare', 'VE3_Dzp', 'VR2_Dzp', 'VR3_Dzp', 'VE3_Dzi', 'VR2_Dzi', 'VR3_Dzi', 'VE3_D', 'VR2_D', 'VR3_D', 'VAdjMat', 'TMWC10', 'SRW02'] 
     [1;46mCluster[0m 19 ['MWC02', 'MWC03', 'SRW04'] 
     [1;46mCluster[0m 20 ['VE1_A', 'VE2_A', 'VE2_DzZ', 'VE2_Dzm', 'VE2_Dzv', 'VE2_Dzse', 'VE2_Dzpe', 'VE2_Dzare', 'VE2_Dzp', 'VE2_Dzi', 'VE2_D'] 
     [1;46mCluster[0m 21 ['piPC1'] 
     [1;46mCluster[0m 22 ['ATSC0d'] 
     [1;46mCluster[0m 23 ['AETA_eta_R'] 
     [1;46mCluster[0m 24 ['ETA_alpha', 'LabuteASA'] 
     [1;46mCluster[0m 25 ['ATS0v', 'SMR'] 
     [1;46mCluster[0m 26 ['ATS1p', 'Sp', 'apol'] 
     [1;46mCluster[0m 27 ['VMcGowan'] 
     [1;46mCluster[0m 28 ['ATS1v'] 
     [1;46mCluster[0m 29 ['Sv'] 
     [1;46mCluster[0m 30 ['ATS4v'] 
     [1;46mCluster[0m 31 ['ATS2d', 'ATS3d'] 
     [1;46mCluster[0m 32 ['ATS2v', 'ATS3v'] 
     [1;46mCluster[0m 33 ['ATS2p', 'ATS3p'] 
     [1;46mCluster[0m 34 ['ATS4p'] 
     [1;46mCluster[0m 35 ['SpMAD_D', 'mZagreb2'] 
     [1;46mCluster[0m 36 ['Kier1'] 
     [1;46mCluster[0m 37 ['VR1_DzZ', 'VR1_Dzm', 'VR1_Dzv', 'VR1_Dzse', 'VR1_Dzpe', 'VR1_Dzare', 'VR1_Dzp', 'VR1_Dzi', 'VR1_D'] 
     [1;46mCluster[0m 38 ['Xp-0d', 'ETA_eta_R'] 
     [1;46mCluster[0m 39 ['ABCGG', 'ETA_beta_s'] 
     [1;46mCluster[0m 40 ['Xp-2d', 'Zagreb1'] 
     [1;46mCluster[0m 41 ['VR3_A', 'VE1_DzZ', 'VE1_Dzm', 'VE1_Dzv', 'VE1_Dzse', 'VE1_Dzpe', 'VE1_Dzare', 'VE1_Dzp', 'VE1_Dzi', 'VE1_D'] 
     [1;46mCluster[0m 42 ['ABC', 'SpAbs_A', 'SpAD_A', 'nHeavyAtom', 'ATS0d', 'ATS1d', 'nBondsO', 'Xp-1d', 'ETA_eta_RL', 'MID', 'MWC01'] 
     [1;46mCluster[0m 43 ['ATS2dv', 'ATS3dv'] 
     [1;46mCluster[0m 44 ['ATS1dv'] 
     [1;46mCluster[0m 45 ['ETA_eta_F'] 
     [1;46mCluster[0m 46 ['BertzCT'] 
     [1;46mCluster[0m 47 ['ETA_beta'] 
     [1;46mCluster[0m 48 ['ETA_eta_FL'] 
     [1;46mCluster[0m 49 ['GGI4'] 
     [1;46mCluster[0m 50 ['ATS8Z', 'ATS8m', 'GGI7'] 
     [1;46mCluster[0m 51 ['ATS7Z', 'ATS7m'] 
     [1;46mCluster[0m 52 ['GGI6'] 
     [1;46mCluster[0m 53 ['ATS7se', 'ATS7pe', 'ATS7are', 'ATS7i'] 
     [1;46mCluster[0m 54 ['ATS7d', 'ATS7v', 'ATS7p'] 
     [1;46mCluster[0m 55 ['ATS8se', 'ATS8pe', 'ATS8are', 'ATS8i'] 
     [1;46mCluster[0m 56 ['ATS8d', 'ATS8v', 'ATS8p'] 
     [1;46mCluster[0m 57 ['ATS8s'] 
     [1;46mCluster[0m 58 ['ATS7s'] 
     [1;46mCluster[0m 59 ['ATS2i', 'nBondsKS'] 
     [1;46mCluster[0m 60 ['ATS3i'] 
     [1;46mCluster[0m 61 ['ATSC0v'] 
     [1;46mCluster[0m 62 ['ATS6se', 'ATS6pe', 'ATS6are'] 
     [1;46mCluster[0m 63 ['ATS6i'] 
     [1;46mCluster[0m 64 ['ATSC0i'] 
     [1;46mCluster[0m 65 ['SpMAD_Dzv'] 
     [1;46mCluster[0m 66 ['SpMAD_Dzp'] 
     [1;46mCluster[0m 67 ['SpMAD_Dzse', 'SpMAD_Dzpe', 'SpMAD_Dzare', 'SpMAD_Dzi'] 
     [1;46mCluster[0m 68 ['ATS6d', 'ATS6v', 'ATS6p'] 
     [1;46mCluster[0m 69 ['ATS5d', 'ATS5v'] 
     [1;46mCluster[0m 70 ['ATS5p'] 
     [1;46mCluster[0m 71 ['TIC1', 'TIC2'] 
     [1;46mCluster[0m 72 ['TIC0'] 
     [1;46mCluster[0m 73 ['ATS3se', 'ATS3pe', 'ATS3are'] 
     [1;46mCluster[0m 74 ['ATS4d'] 
     [1;46mCluster[0m 75 ['ATS4se', 'ATS4pe', 'ATS4are'] 
     [1;46mCluster[0m 76 ['ATS4i'] 
     [1;46mCluster[0m 77 ['ATS0se', 'ATS0pe', 'ATS0are'] 
     [1;46mCluster[0m 78 ['ATS1se', 'ATS2se', 'ATS1pe', 'ATS2pe', 'ATS1are', 'ATS2are', 'ATS1i', 'nBonds', 'Sse', 'Spe', 'Sare'] 
     [1;46mCluster[0m 79 ['nAtom', 'ATS0i', 'Si'] 
     [1;46mCluster[0m 80 ['ATS5se', 'ATS5pe', 'ATS5are', 'ATS5i'] 
     [1;46mCluster[0m 81 ['SpAbs_Dzv', 'SpMax_Dzv', 'SpDiam_Dzv', 'SpAD_Dzv', 'LogEE_Dzv', 'SpAbs_Dzp', 'SpMax_Dzp', 'SpDiam_Dzp', 'SpAD_Dzp', 'LogEE_Dzp'] 
     [1;46mCluster[0m 82 ['SpAbs_Dzse', 'SpMax_Dzse', 'SpDiam_Dzse', 'SpAD_Dzse', 'LogEE_Dzse', 'SpAbs_Dzpe', 'SpMax_Dzpe', 'SpDiam_Dzpe', 'SpAD_Dzpe', 'LogEE_Dzpe', 'SpAbs_Dzare', 'SpMax_Dzare', 'SpDiam_Dzare', 'SpAD_Dzare', 'LogEE_Dzare', 'SpAbs_Dzi', 'SpMax_Dzi', 'SpDiam_Dzi', 'SpAD_Dzi', 'LogEE_Dzi', 'SpAbs_D', 'SpMax_D', 'SpDiam_D', 'SpAD_D', 'LogEE_D'] 
     [1;46mCluster[0m 83 ['SpAbs_DzZ', 'SpMax_DzZ', 'SpDiam_DzZ', 'SpAD_DzZ', 'LogEE_DzZ', 'SpAbs_Dzm', 'SpMax_Dzm', 'SpDiam_Dzm', 'SpAD_Dzm', 'LogEE_Dzm', 'ECIndex'] 
     [1;46mCluster[0m 84 ['WPath'] 
     [1;46mCluster[0m 85 ['Diameter'] 
     [1;46mCluster[0m 86 ['Radius'] 
     [1;46mCluster[0m 87 ['SpMAD_DzZ', 'SpMAD_Dzm'] 
     [1;46mCluster[0m 88 ['TIC3', 'TIC4', 'TIC5'] 
     [1;46mCluster[0m 89 ['GGI1'] 
     [1;46mCluster[0m 90 ['mZagreb1'] 
     [1;46mCluster[0m 91 ['ATS0p', 'Xp-1dv'] 
     [1;46mCluster[0m 92 ['Xp-0dv'] 
     [1;46mCluster[0m 93 ['ATS1Z', 'ATS1m'] 
     [1;46mCluster[0m 94 ['SZ', 'Sm', 'MW'] 
     [1;46mCluster[0m 95 ['ZMIC0'] 
     [1;46mCluster[0m 96 ['ETA_eta'] 
     [1;46mCluster[0m 97 ['ETA_eta_L'] 
     [1;46mCluster[0m 98 ['ATS6Z', 'ATS6m'] 
     [1;46mCluster[0m 99 ['ATS5Z', 'ATS5m'] 
     [1;46mCluster[0m 100 ['Xp-6dv'] 
     [1;46mCluster[0m 101 ['Xp-7dv'] 
     [1;46mCluster[0m 102 ['Xp-3dv'] 
     [1;46mCluster[0m 103 ['Xp-4dv'] 
     [1;46mCluster[0m 104 ['Xp-5dv'] 
     [1;46mCluster[0m 105 ['Xp-2dv'] 
     [1;46mCluster[0m 106 ['ATS2Z'] 
     [1;46mCluster[0m 107 ['ATS2m'] 
     [1;46mCluster[0m 108 ['Xc-3d'] 
     [1;46mCluster[0m 109 ['ETA_eta_BR'] 
     [1;46mCluster[0m 110 ['ATS4Z', 'ATS4m'] 
     [1;46mCluster[0m 111 ['Xpc-5d', 'GGI3'] 
     [1;46mCluster[0m 112 ['Xpc-4d'] 
     [1;46mCluster[0m 113 ['GGI2'] 
     [1;46mCluster[0m 114 ['ATS3Z', 'ATS3m'] 
     [1;46mCluster[0m 115 ['Xpc-6d'] 
     [1;46mCluster[0m 116 ['MWC04', 'MWC05', 'MWC06', 'MWC07', 'MWC08', 'MWC09', 'MWC10', 'SRW06', 'SRW08', 'SRW10'] 
     [1;46mCluster[0m 117 ['SpMax_A', 'SpDiam_A'] 
     [1;46mCluster[0m 118 ['piPC2'] 
     [1;46mCluster[0m 119 ['IC3', 'IC4', 'IC5'] 
     [1;46mCluster[0m 120 ['piPC3', 'piPC4'] 
     [1;46mCluster[0m 121 ['piPC5', 'piPC6'] 
     [1;46mCluster[0m 122 ['piPC7'] 
     [1;46mCluster[0m 123 ['TpiPC10'] 
     [1;46mCluster[0m 124 ['nBondsKD'] 
     [1;46mCluster[0m 125 ['ETA_beta_ns'] 
     [1;46mCluster[0m 126 ['nBondsM'] 
     [1;46mCluster[0m 127 ['AETA_eta_F'] 
     [1;46mCluster[0m 128 ['FilterItLogS'] 
     [1;46mCluster[0m 129 ['SLogP'] 
     [1;46mCluster[0m 130 ['ZMIC1'] 
     [1;46mCluster[0m 131 ['piPC9'] 
     [1;46mCluster[0m 132 ['piPC10'] 
     [1;46mCluster[0m 133 ['piPC8'] 
     [1;46mCluster[0m 134 ['Xp-6d'] 
     [1;46mCluster[0m 135 ['MPC4'] 
     [1;46mCluster[0m 136 ['Xp-5d'] 
     [1;46mCluster[0m 137 ['MPC5'] 
     [1;46mCluster[0m 138 ['Xp-7d'] 
     [1;46mCluster[0m 139 ['MPC6', 'TMPC10'] 
     [1;46mCluster[0m 140 ['MPC7'] 
     [1;46mCluster[0m 141 ['fragCpx'] 
     [1;46mCluster[0m 142 ['AATS0i', 'Mi'] 
     [1;46mCluster[0m 143 ['AATS1i'] 
     [1;46mCluster[0m 144 ['AETA_eta_L'] 
     [1;46mCluster[0m 145 ['AETA_eta_FL'] 
     [1;46mCluster[0m 146 ['FCSP3'] 
     [1;46mCluster[0m 147 ['AETA_beta_ns'] 
     [1;46mCluster[0m 148 ['AETA_beta'] 
     [1;46mCluster[0m 149 ['AETA_dBeta'] 
     [1;46mCluster[0m 150 ['AATSC1v'] 
     [1;46mCluster[0m 151 ['MATS1v'] 
     [1;46mCluster[0m 152 ['AATS2i'] 
     [1;46mCluster[0m 153 ['GATS1i'] 
     [1;46mCluster[0m 154 ['ETA_dEpsilon_B'] 
     [1;46mCluster[0m 155 ['fMF'] 
     [1;46mCluster[0m 156 ['GATS1v'] 
     [1;46mCluster[0m 157 ['GATS1p'] 
     [1;46mCluster[0m 158 ['AETA_beta_ns_d'] 
     [1;46mCluster[0m 159 ['SlogP_VSA7'] 
     [1;46mCluster[0m 160 ['CIC4', 'CIC5'] 
     [1;46mCluster[0m 161 ['CIC3'] 
     [1;46mCluster[0m 162 ['SIC4', 'SIC5'] 
     [1;46mCluster[0m 163 ['BIC4', 'BIC5'] 
     [1;46mCluster[0m 164 ['SIC3'] 
     [1;46mCluster[0m 165 ['BIC3'] 
     [1;46mCluster[0m 166 ['AATSC1i'] 
     [1;46mCluster[0m 167 ['MATS1i'] 
     [1;46mCluster[0m 168 ['AATSC1p'] 
     [1;46mCluster[0m 169 ['MATS1p'] 
     [1;46mCluster[0m 170 ['JGI4'] 
     [1;46mCluster[0m 171 ['AATSC0se'] 
     [1;46mCluster[0m 172 ['AATSC0pe'] 
     [1;46mCluster[0m 173 ['AATS0are', 'Mare'] 
     [1;46mCluster[0m 174 ['AATS1are'] 
     [1;46mCluster[0m 175 ['AATS2are'] 
     [1;46mCluster[0m 176 ['IC0'] 
     [1;46mCluster[0m 177 ['AATS0se', 'Mse'] 
     [1;46mCluster[0m 178 ['AATS0pe', 'Mpe'] 
     [1;46mCluster[0m 179 ['ETA_epsilon_1', 'ETA_dEpsilon_A'] 
     [1;46mCluster[0m 180 ['AATS1se'] 
     [1;46mCluster[0m 181 ['AATS1pe'] 
     [1;46mCluster[0m 182 ['AATS2pe'] 
     [1;46mCluster[0m 183 ['ETA_epsilon_4', 'ETA_dEpsilon_C'] 
     [1;46mCluster[0m 184 ['AATS2se'] 
     [1;46mCluster[0m 185 ['ETA_epsilon_2'] 
     [1;46mCluster[0m 186 ['ETA_epsilon_5'] 
     [1;46mCluster[0m 187 ['AMID_h'] 
     [1;46mCluster[0m 188 ['MIC0'] 
     [1;46mCluster[0m 189 ['SIC2'] 
     [1;46mCluster[0m 190 ['BIC2'] 
     [1;46mCluster[0m 191 ['SIC1'] 
     [1;46mCluster[0m 192 ['BIC1'] 
     [1;46mCluster[0m 193 ['AETA_eta_RL'] 
     [1;46mCluster[0m 194 ['JGI1'] 
     [1;46mCluster[0m 195 ['BCUTd-1l'] 
     [1;46mCluster[0m 196 ['EState_VSA8'] 
     [1;46mCluster[0m 197 ['AATSC0v'] 
     [1;46mCluster[0m 198 ['AATS0v'] 
     [1;46mCluster[0m 199 ['Mv'] 
     [1;46mCluster[0m 200 ['AATS1v'] 
     [1;46mCluster[0m 201 ['AATS2v'] 
     [1;46mCluster[0m 202 ['nX', 'MID_X'] 
     [1;46mCluster[0m 203 ['nCl', 'NsCl', 'SsCl'] 
     [1;46mCluster[0m 204 ['EState_VSA9'] 
     [1;46mCluster[0m 205 ['BCUTZ-1h', 'BCUTm-1h'] 
     [1;46mCluster[0m 206 ['AATSC0p'] 
     [1;46mCluster[0m 207 ['BCUTv-1h'] 
     [1;46mCluster[0m 208 ['BCUTdv-1l'] 
     [1;46mCluster[0m 209 ['MZ', 'Mm', 'AMW'] 
     [1;46mCluster[0m 210 ['AATS1Z', 'AATS1m'] 
     [1;46mCluster[0m 211 ['AATS1p'] 
     [1;46mCluster[0m 212 ['AATS2p'] 
     [1;46mCluster[0m 213 ['AATS0p'] 
     [1;46mCluster[0m 214 ['Mp'] 
     [1;46mCluster[0m 215 ['AMID_X'] 
     [1;46mCluster[0m 216 ['nBr', 'NsBr', 'SsBr', 'ETA_dPsi_B'] 
     [1;46mCluster[0m 217 ['AETA_alpha'] 
     [1;46mCluster[0m 218 ['ETA_dAlpha_A'] 
     [1;46mCluster[0m 219 ['AATSC0Z', 'AATSC0m'] 
     [1;46mCluster[0m 220 ['AATS2Z', 'AATS2m'] 
     [1;46mCluster[0m 221 ['AATS0Z', 'AATS0m'] 
     [1;46mCluster[0m 222 ['MATS1s'] 
     [1;46mCluster[0m 223 ['ETA_dEpsilon_D'] 
     [1;46mCluster[0m 224 ['ATSC1se', 'ATSC1are'] 
     [1;46mCluster[0m 225 ['ATSC1pe'] 
     [1;46mCluster[0m 226 ['MATS1pe', 'MATS1are'] 
     [1;46mCluster[0m 227 ['MATS1se'] 
     [1;46mCluster[0m 228 ['AATSC1se'] 
     [1;46mCluster[0m 229 ['AATSC1pe'] 
     [1;46mCluster[0m 230 ['AATSC1are'] 
     [1;46mCluster[0m 231 ['GATS1se', 'GATS1pe'] 
     [1;46mCluster[0m 232 ['GATS1are'] 
     [1;46mCluster[0m 233 ['GATS1s'] 
     [1;46mCluster[0m 234 ['GATS2pe'] 
     [1;46mCluster[0m 235 ['GATS2are'] 
     [1;46mCluster[0m 236 ['GATS2se'] 
     [1;46mCluster[0m 237 ['GATS2s'] 
     [1;46mCluster[0m 238 ['GATS2dv'] 
     [1;46mCluster[0m 239 ['ATSC3i'] 
     [1;46mCluster[0m 240 ['SdssC'] 
     [1;46mCluster[0m 241 ['AATSC0are'] 
     [1;46mCluster[0m 242 ['BCUTs-1h'] 
     [1;46mCluster[0m 243 ['AMID_O'] 
     [1;46mCluster[0m 244 ['AATS0s', 'AATSC0s'] 
     [1;46mCluster[0m 245 ['NddsN', 'SddsN'] 
     [1;46mCluster[0m 246 ['PEOE_VSA13'] 
     [1;46mCluster[0m 247 ['VSA_EState3'] 
     [1;46mCluster[0m 248 ['ATSC5s'] 
     [1;46mCluster[0m 249 ['AATSC2s'] 
     [1;46mCluster[0m 250 ['PEOE_VSA2'] 
     [1;46mCluster[0m 251 ['AATS2s'] 
     [1;46mCluster[0m 252 ['SlogP_VSA4'] 
     [1;46mCluster[0m 253 ['AATSC0c'] 
     [1;46mCluster[0m 254 ['ETA_dPsi_A'] 
     [1;46mCluster[0m 255 ['ETA_psi_1'] 
     [1;46mCluster[0m 256 ['AATSC1s'] 
     [1;46mCluster[0m 257 ['ETA_dAlpha_B'] 
     [1;46mCluster[0m 258 ['AATS1s'] 
     [1;46mCluster[0m 259 ['BalabanJ'] 
     [1;46mCluster[0m 260 ['ETA_dBeta'] 
     [1;46mCluster[0m 261 ['NsCH3', 'SsCH3'] 
     [1;46mCluster[0m 262 ['SlogP_VSA2'] 
     [1;46mCluster[0m 263 ['NssCH2'] 
     [1;46mCluster[0m 264 ['SssCH2'] 
     [1;46mCluster[0m 265 ['C2SP3'] 
     [1;46mCluster[0m 266 ['C1SP3'] 
     [1;46mCluster[0m 267 ['VSA_EState8'] 
     [1;46mCluster[0m 268 ['RotRatio'] 
     [1;46mCluster[0m 269 ['ATSC1i'] 
     [1;46mCluster[0m 270 ['SMR_VSA6'] 
     [1;46mCluster[0m 271 ['ATSC1v'] 
     [1;46mCluster[0m 272 ['n3HRing', 'n3AHRing'] 
     [1;46mCluster[0m 273 ['nAHRing'] 
     [1;46mCluster[0m 274 ['SRW05', 'SRW07', 'SRW09'] 
     [1;46mCluster[0m 275 ['Xch-3dv', 'n3Ring', 'n3ARing', 'SRW03'] 
     [1;46mCluster[0m 276 ['Xch-3d'] 
     [1;46mCluster[0m 277 ['NdssS', 'SdssS', 'n7Ring', 'n7HRing', 'n7ARing', 'n7AHRing', 'n12FRing', 'n12FHRing', 'n12FARing', 'n12FAHRing'] 
     [1;46mCluster[0m 278 ['nG12FARing', 'nG12FAHRing'] 
     [1;46mCluster[0m 279 ['nFARing', 'nFAHRing'] 
     [1;46mCluster[0m 280 ['Xc-6d'] 
     [1;46mCluster[0m 281 ['Xc-6dv'] 
     [1;46mCluster[0m 282 ['nBridgehead'] 
     [1;46mCluster[0m 283 ['n5ARing'] 
     [1;46mCluster[0m 284 ['nARing'] 
     [1;46mCluster[0m 285 ['n5Ring'] 
     [1;46mCluster[0m 286 ['Xc-5dv'] 
     [1;46mCluster[0m 287 ['nG12FRing', 'nG12FHRing'] 
     [1;46mCluster[0m 288 ['nFRing'] 
     [1;46mCluster[0m 289 ['nFHRing'] 
     [1;46mCluster[0m 290 ['SlogP_VSA8'] 
     [1;46mCluster[0m 291 ['ATSC2c'] 
     [1;46mCluster[0m 292 ['MATS2c'] 
     [1;46mCluster[0m 293 ['PEOE_VSA1'] 
     [1;46mCluster[0m 294 ['nN', 'MID_N'] 
     [1;46mCluster[0m 295 ['BCUTv-1l'] 
     [1;46mCluster[0m 296 ['BCUTi-1h'] 
     [1;46mCluster[0m 297 ['AATSC0i'] 
     [1;46mCluster[0m 298 ['ATSC3dv'] 
     [1;46mCluster[0m 299 ['ATSC6dv'] 
     [1;46mCluster[0m 300 ['TopoShapeIndex', 'PetitjeanIndex'] 
     [1;46mCluster[0m 301 ['PEOE_VSA3'] 
     [1;46mCluster[0m 302 ['PEOE_VSA12'] 
     [1;46mCluster[0m 303 ['AATSC2d'] 
     [1;46mCluster[0m 304 ['MATS1c'] 
     [1;46mCluster[0m 305 ['GATS1c'] 
     [1;46mCluster[0m 306 ['AATSC1c'] 
     [1;46mCluster[0m 307 ['NssNH', 'SssNH'] 
     [1;46mCluster[0m 308 ['ATSC3c'] 
     [1;46mCluster[0m 309 ['AATSC2c'] 
     [1;46mCluster[0m 310 ['SlogP_VSA1'] 
     [1;46mCluster[0m 311 ['SMR_VSA3'] 
     [1;46mCluster[0m 312 ['ATSC7pe', 'ATSC7are'] 
     [1;46mCluster[0m 313 ['ATSC7se'] 
     [1;46mCluster[0m 314 ['ATSC7s'] 
     [1;46mCluster[0m 315 ['ATSC6pe'] 
     [1;46mCluster[0m 316 ['ATSC6are'] 
     [1;46mCluster[0m 317 ['ATSC6se'] 
     [1;46mCluster[0m 318 ['ATSC6c'] 
     [1;46mCluster[0m 319 ['ATSC6s'] 
     [1;46mCluster[0m 320 ['ATSC4s'] 
     [1;46mCluster[0m 321 ['ATSC4c'] 
     [1;46mCluster[0m 322 ['NaaN', 'SaaN'] 
     [1;46mCluster[0m 323 ['n6HRing', 'n6aHRing'] 
     [1;46mCluster[0m 324 ['naHRing'] 
     [1;46mCluster[0m 325 ['nHRing'] 
     [1;46mCluster[0m 326 ['PEOE_VSA11'] 
     [1;46mCluster[0m 327 ['VSA_EState5'] 
     [1;46mCluster[0m 328 ['VSA_EState7'] 
     [1;46mCluster[0m 329 ['NaasN', 'SaasN'] 
     [1;46mCluster[0m 330 ['n6AHRing'] 
     [1;46mCluster[0m 331 ['nBase'] 
     [1;46mCluster[0m 332 ['n6ARing'] 
     [1;46mCluster[0m 333 ['NddC', 'SddC'] 
     [1;46mCluster[0m 334 ['NdsN', 'SdsN'] 
     [1;46mCluster[0m 335 ['NsSH', 'SsSH'] 
     [1;46mCluster[0m 336 ['nAcid'] 
     [1;46mCluster[0m 337 ['NdCH2', 'SdCH2'] 
     [1;46mCluster[0m 338 ['NdsCH', 'SdsCH'] 
     [1;46mCluster[0m 339 ['NddssS', 'SddssS'] 
     [1;46mCluster[0m 340 ['n5AHRing'] 
     [1;46mCluster[0m 341 ['n9FARing', 'n9FAHRing'] 
     [1;46mCluster[0m 342 ['n9FRing', 'n9FHRing'] 
     [1;46mCluster[0m 343 ['n9FaRing', 'n9FaHRing'] 
     [1;46mCluster[0m 344 ['NaaNH', 'SaaNH'] 
     [1;46mCluster[0m 345 ['NaaS', 'SaaS'] 
     [1;46mCluster[0m 346 ['n5aRing', 'n5aHRing'] 
     [1;46mCluster[0m 347 ['n5HRing'] 
     [1;46mCluster[0m 348 ['NaaO', 'SaaO'] 
     [1;46mCluster[0m 349 ['n10FRing', 'n10FaRing'] 
     [1;46mCluster[0m 350 ['n10FHRing', 'n10FaHRing'] 
     [1;46mCluster[0m 351 ['nG12FaRing', 'nG12FaHRing'] 
     [1;46mCluster[0m 352 ['NaaaC', 'SaaaC', 'nFaRing'] 
     [1;46mCluster[0m 353 ['nFaHRing'] 
     [1;46mCluster[0m 354 ['ATSC3se'] 
     [1;46mCluster[0m 355 ['ATSC3are'] 
     [1;46mCluster[0m 356 ['ATSC3pe'] 
     [1;46mCluster[0m 357 ['ATSC3s'] 
     [1;46mCluster[0m 358 ['ATSC5c'] 
     [1;46mCluster[0m 359 ['nHBDon'] 
     [1;46mCluster[0m 360 ['AMID_N'] 
     [1;46mCluster[0m 361 ['NsOH', 'SsOH'] 
     [1;46mCluster[0m 362 ['NsNH2', 'SsNH2'] 
     [1;46mCluster[0m 363 ['VSA_EState4'] 
     [1;46mCluster[0m 364 ['ATSC1p'] 
     [1;46mCluster[0m 365 ['ATSC4i'] 
     [1;46mCluster[0m 366 ['PEOE_VSA5'] 
     [1;46mCluster[0m 367 ['EState_VSA5'] 
     [1;46mCluster[0m 368 ['nI', 'C2SP1', 'NsI', 'SsI'] 
     [1;46mCluster[0m 369 ['NtsC', 'StsC'] 
     [1;46mCluster[0m 370 ['NtN', 'StN', 'SMR_VSA2'] 
     [1;46mCluster[0m 371 ['nBondsT'] 
     [1;46mCluster[0m 372 ['C1SP1'] 
     [1;46mCluster[0m 373 ['C1SP2'] 
     [1;46mCluster[0m 374 ['PEOE_VSA4'] 
     [1;46mCluster[0m 375 ['SsssCH'] 
     [1;46mCluster[0m 376 ['ATSC7i'] 
     [1;46mCluster[0m 377 ['ATSC7c'] 
     [1;46mCluster[0m 378 ['ATSC8se', 'ATSC8pe', 'ATSC8are'] 
     [1;46mCluster[0m 379 ['ATSC8s'] 
     [1;46mCluster[0m 380 ['ATSC8dv'] 
     [1;46mCluster[0m 381 ['ATSC8c'] 
     [1;46mCluster[0m 382 ['ATSC8Z', 'ATSC8m'] 
     [1;46mCluster[0m 383 ['ATSC8v'] 
     [1;46mCluster[0m 384 ['ATSC8p'] 
     [1;46mCluster[0m 385 ['ATSC7d'] 
     [1;46mCluster[0m 386 ['ATSC8d'] 
     [1;46mCluster[0m 387 ['NssS', 'SssS'] 
     [1;46mCluster[0m 388 ['NdS', 'SdS'] 
     [1;46mCluster[0m 389 ['nS'] 
     [1;46mCluster[0m 390 ['BCUTp-1h'] 
     [1;46mCluster[0m 391 ['BCUTi-1l'] 
     [1;46mCluster[0m 392 ['ATSC4Z', 'ATSC4m'] 
     [1;46mCluster[0m 393 ['ATSC4p'] 
     [1;46mCluster[0m 394 ['ATSC4v'] 
     [1;46mCluster[0m 395 ['NsssN', 'SsssN'] 
     [1;46mCluster[0m 396 ['VSA_EState9'] 
     [1;46mCluster[0m 397 ['ATSC5p'] 
     [1;46mCluster[0m 398 ['ATSC5i'] 
     [1;46mCluster[0m 399 ['ATSC5v'] 
     [1;46mCluster[0m 400 ['GATS2i'] 
     [1;46mCluster[0m 401 ['GATS2c'] 
     [1;46mCluster[0m 402 ['ATSC7Z', 'ATSC7m'] 
     [1;46mCluster[0m 403 ['ATSC7v'] 
     [1;46mCluster[0m 404 ['ATSC7p'] 
     [1;46mCluster[0m 405 ['ATSC3Z', 'ATSC3m'] 
     [1;46mCluster[0m 406 ['ATSC3v'] 
     [1;46mCluster[0m 407 ['ATSC3p'] 
     [1;46mCluster[0m 408 ['ATSC2d'] 
     [1;46mCluster[0m 409 ['SpMAD_A'] 
     [1;46mCluster[0m 410 ['AXp-0d'] 
     [1;46mCluster[0m 411 ['ETA_shape_p'] 
     [1;46mCluster[0m 412 ['ETA_epsilon_3'] 
     [1;46mCluster[0m 413 ['AMID'] 
     [1;46mCluster[0m 414 ['AXp-1d'] 
     [1;46mCluster[0m 415 ['nRing'] 
     [1;46mCluster[0m 416 ['BCUTZ-1l', 'BCUTm-1l'] 
     [1;46mCluster[0m 417 ['BCUTse-1l'] 
     [1;46mCluster[0m 418 ['SaaCH', 'VSA_EState6'] 
     [1;46mCluster[0m 419 ['NaaCH', 'SlogP_VSA6'] 
     [1;46mCluster[0m 420 ['nAromAtom', 'nAromBond', 'nBondsA', 'n6Ring', 'naRing', 'n6aRing'] 
     [1;46mCluster[0m 421 ['C2SP2'] 
     [1;46mCluster[0m 422 ['SMR_VSA7'] 
     [1;46mCluster[0m 423 ['AATS1dv'] 
     [1;46mCluster[0m 424 ['AATS2dv'] 
     [1;46mCluster[0m 425 ['AATS0dv'] 
     [1;46mCluster[0m 426 ['MIC3', 'MIC4', 'MIC5'] 
     [1;46mCluster[0m 427 ['MIC2'] 
     [1;46mCluster[0m 428 ['MIC1'] 
     [1;46mCluster[0m 429 ['AATSC0d'] 
     [1;46mCluster[0m 430 ['AATS0d', 'AATS1d'] 
     [1;46mCluster[0m 431 ['AATS2d'] 
     [1;46mCluster[0m 432 ['NaasC'] 
     [1;46mCluster[0m 433 ['ETA_shape_y'] 
     [1;46mCluster[0m 434 ['nBondsS'] 
     [1;46mCluster[0m 435 ['bpol'] 
     [1;46mCluster[0m 436 ['ATSC0p'] 
     [1;46mCluster[0m 437 ['nRot'] 
     [1;46mCluster[0m 438 ['nH'] 
     [1;46mCluster[0m 439 ['ATSC1d'] 
     [1;46mCluster[0m 440 ['CIC0'] 
     [1;46mCluster[0m 441 ['ATS7dv'] 
     [1;46mCluster[0m 442 ['ATS8dv'] 
     [1;46mCluster[0m 443 ['GGI8'] 
     [1;46mCluster[0m 444 ['GGI9'] 
     [1;46mCluster[0m 445 ['GGI10'] 
     [1;46mCluster[0m 446 ['JGI10'] 
     [1;46mCluster[0m 447 ['NssO', 'SssO'] 
     [1;46mCluster[0m 448 ['SlogP_VSA3'] 
     [1;46mCluster[0m 449 ['VSA_EState1'] 
     [1;46mCluster[0m 450 ['ATSC1c'] 
     [1;46mCluster[0m 451 ['JGI7'] 
     [1;46mCluster[0m 452 ['JGI8'] 
     [1;46mCluster[0m 453 ['JGI6'] 
     [1;46mCluster[0m 454 ['PEOE_VSA7'] 
     [1;46mCluster[0m 455 ['RNCG'] 
     [1;46mCluster[0m 456 ['ATSC1dv'] 
     [1;46mCluster[0m 457 ['BCUTc-1l'] 
     [1;46mCluster[0m 458 ['ATSC2pe', 'ATSC2are'] 
     [1;46mCluster[0m 459 ['ATSC2se'] 
     [1;46mCluster[0m 460 ['AATSC2dv'] 
     [1;46mCluster[0m 461 ['MATS2dv'] 
     [1;46mCluster[0m 462 ['ATSC2dv'] 
     [1;46mCluster[0m 463 ['MATS2pe', 'MATS2are'] 
     [1;46mCluster[0m 464 ['AATSC2se'] 
     [1;46mCluster[0m 465 ['AATSC2pe'] 
     [1;46mCluster[0m 466 ['AATSC2are'] 
     [1;46mCluster[0m 467 ['MATS2se'] 
     [1;46mCluster[0m 468 ['MATS2s'] 
     [1;46mCluster[0m 469 ['ATS0s'] 
     [1;46mCluster[0m 470 ['ATS6s'] 
     [1;46mCluster[0m 471 ['nO', 'MID_O'] 
     [1;46mCluster[0m 472 ['SMR_VSA1'] 
     [1;46mCluster[0m 473 ['ATSC0c'] 
     [1;46mCluster[0m 474 ['nHBAcc'] 
     [1;46mCluster[0m 475 ['EState_VSA1'] 
     [1;46mCluster[0m 476 ['nBondsD'] 
     [1;46mCluster[0m 477 ['TopoPSA'] 
     [1;46mCluster[0m 478 ['SlogP_VSA10'] 
     [1;46mCluster[0m 479 ['JGI5'] 
     [1;46mCluster[0m 480 ['SM1_Dzv', 'SM1_Dzp'] 
     [1;46mCluster[0m 481 ['ATSC0s'] 
     [1;46mCluster[0m 482 ['EState_VSA10'] 
     [1;46mCluster[0m 483 ['ATSC2s'] 
     [1;46mCluster[0m 484 ['NdO', 'SdO'] 
     [1;46mCluster[0m 485 ['VSA_EState2'] 
     [1;46mCluster[0m 486 ['TopoPSA(NO)'] 
     [1;46mCluster[0m 487 ['AATSC0dv'] 
     [1;46mCluster[0m 488 ['BCUTdv-1h'] 
     [1;46mCluster[0m 489 ['BCUTse-1h'] 
     [1;46mCluster[0m 490 ['BCUTpe-1h'] 
     [1;46mCluster[0m 491 ['BCUTare-1h'] 
     [1;46mCluster[0m 492 ['BCUTp-1l'] 
     [1;46mCluster[0m 493 ['BCUTpe-1l'] 
     [1;46mCluster[0m 494 ['BCUTare-1l'] 
     [1;46mCluster[0m 495 ['BCUTd-1h'] 
     [1;46mCluster[0m 496 ['AETA_beta_s'] 
     [1;46mCluster[0m 497 ['SM1_Dzse', 'SM1_Dzpe'] 
     [1;46mCluster[0m 498 ['SM1_Dzare'] 
     [1;46mCluster[0m 499 ['SM1_Dzi'] 
     [1;46mCluster[0m 500 ['nHetero', 'MID_h'] 
     [1;46mCluster[0m 501 ['ETA_eta_B'] 
     [1;46mCluster[0m 502 ['IC1'] 
     [1;46mCluster[0m 503 ['IC2'] 
     [1;46mCluster[0m 504 ['BCUTc-1h'] 
     [1;46mCluster[0m 505 ['SM1_DzZ', 'SM1_Dzm'] 
     [1;46mCluster[0m 506 ['ATS0Z'] 
     [1;46mCluster[0m 507 ['ATS0m'] 
     [1;46mCluster[0m 508 ['ZMIC4', 'ZMIC5'] 
     [1;46mCluster[0m 509 ['ZMIC3'] 
     [1;46mCluster[0m 510 ['ZMIC2'] 
     [1;46mCluster[0m 511 ['Xc-3dv'] 
     [1;46mCluster[0m 512 ['ATSC0Z'] 
     [1;46mCluster[0m 513 ['ATSC0m'] 
     [1;46mCluster[0m 514 ['AETA_eta'] 
     [1;46mCluster[0m 515 ['Xch-7d'] 
     [1;46mCluster[0m 516 ['Xch-6dv'] 
     [1;46mCluster[0m 517 ['MPC9', 'MPC10'] 
     [1;46mCluster[0m 518 ['MPC8'] 
     [1;46mCluster[0m 519 ['Xch-6d'] 
     [1;46mCluster[0m 520 ['NssssC'] 
     [1;46mCluster[0m 521 ['SssssC'] 
     [1;46mCluster[0m 522 ['C3SP3'] 
     [1;46mCluster[0m 523 ['Xch-7dv'] 
     [1;46mCluster[0m 524 ['NsssCH'] 
     [1;46mCluster[0m 525 ['Xpc-4dv'] 
     [1;46mCluster[0m 526 ['Xpc-5dv'] 
     [1;46mCluster[0m 527 ['Xpc-6dv'] 
     [1;46mCluster[0m 528 ['Xc-5d'] 
     [1;46mCluster[0m 529 ['ATSC3d'] 
     [1;46mCluster[0m 530 ['AXp-0dv'] 
     [1;46mCluster[0m 531 ['AXp-1dv'] 
     [1;46mCluster[0m 532 ['GATS1dv'] 
     [1;46mCluster[0m 533 ['AATSC2Z', 'AATSC2m'] 
     [1;46mCluster[0m 534 ['AATSC1Z', 'AATSC1m'] 
     [1;46mCluster[0m 535 ['AATSC1dv'] 
     [1;46mCluster[0m 536 ['MATS1dv'] 
     [1;46mCluster[0m 537 ['AATSC1d'] 
     [1;46mCluster[0m 538 ['CIC1'] 
     [1;46mCluster[0m 539 ['CIC2'] 
     [1;46mCluster[0m 540 ['SIC0'] 
     [1;46mCluster[0m 541 ['BIC0'] 
     [1;46mCluster[0m 542 ['RPCG'] 
     [1;46mCluster[0m 543 ['C3SP2'] 
     [1;46mCluster[0m 544 ['SMR_VSA9'] 
     [1;46mCluster[0m 545 ['SaasC'] 
     [1;46mCluster[0m 546 ['EState_VSA7'] 
     [1;46mCluster[0m 547 ['ATSC4d'] 
     [1;46mCluster[0m 548 ['PEOE_VSA6'] 
     [1;46mCluster[0m 549 ['AMID_C'] 
     [1;46mCluster[0m 550 ['ATSC5pe'] 
     [1;46mCluster[0m 551 ['ATSC5are'] 
     [1;46mCluster[0m 552 ['ATSC5se'] 
     [1;46mCluster[0m 553 ['ATSC5dv'] 
     [1;46mCluster[0m 554 ['EState_VSA6'] 
     [1;46mCluster[0m 555 ['JGI9'] 
     [1;46mCluster[0m 556 ['PEOE_VSA10'] 
     [1;46mCluster[0m 557 ['SlogP_VSA11'] 
     [1;46mCluster[0m 558 ['EState_VSA3'] 
     [1;46mCluster[0m 559 ['nF', 'NsF', 'SsF'] 
     [1;46mCluster[0m 560 ['ATSC8i'] 
     [1;46mCluster[0m 561 ['VR1_A', 'VR2_A'] 
     [1;46mCluster[0m 562 ['Xch-5d', 'Xch-5dv'] 
     [1;46mCluster[0m 563 ['Xch-4d', 'Xch-4dv'] 
     [1;46mCluster[0m 564 ['SMR_VSA4'] 
     [1;46mCluster[0m 565 ['C4SP3'] 
     [1;46mCluster[0m 566 ['ATSC6Z', 'ATSC6m'] 
     [1;46mCluster[0m 567 ['ATSC6v'] 
     [1;46mCluster[0m 568 ['ATSC6p'] 
     [1;46mCluster[0m 569 ['ATSC6i'] 
     [1;46mCluster[0m 570 ['ATSC6d'] 
     [1;46mCluster[0m 571 ['NdssC'] 
     [1;46mCluster[0m 572 ['SMR_VSA5'] 
     [1;46mCluster[0m 573 ['SlogP_VSA5'] 
     [1;46mCluster[0m 574 ['EState_VSA4'] 
     [1;46mCluster[0m 575 ['ATSC7dv'] 
     [1;46mCluster[0m 576 ['ATSC5d'] 
     [1;46mCluster[0m 577 ['AETA_eta_BR'] 
     [1;46mCluster[0m 578 ['JGT10'] 
     [1;46mCluster[0m 579 ['AETA_eta_B'] 
     [1;46mCluster[0m 580 ['GATS2v'] 
     [1;46mCluster[0m 581 ['GATS2p'] 
     [1;46mCluster[0m 582 ['JGI2'] 
     [1;46mCluster[0m 583 ['PEOE_VSA9'] 
     [1;46mCluster[0m 584 ['EState_VSA2'] 
     [1;46mCluster[0m 585 ['GATS2Z', 'GATS2m'] 
     [1;46mCluster[0m 586 ['JGI3'] 
     [1;46mCluster[0m 587 ['GATS1Z', 'GATS1m'] 
     [1;46mCluster[0m 588 ['ETA_beta_ns_d'] 
     [1;46mCluster[0m 589 ['ATSC5Z', 'ATSC5m'] 
     [1;46mCluster[0m 590 ['Xc-4d'] 
     [1;46mCluster[0m 591 ['ETA_shape_x'] 
     [1;46mCluster[0m 592 ['Xc-4dv'] 
     [1;46mCluster[0m 593 ['AATSC2v', 'MATS2v', 'MATS2p'] 
     [1;46mCluster[0m 594 ['AATSC2p'] 
     [1;46mCluster[0m 595 ['ATSC2v'] 
     [1;46mCluster[0m 596 ['ATSC2p'] 
     [1;46mCluster[0m 597 ['AATSC2i', 'MATS2i'] 
     [1;46mCluster[0m 598 ['ATSC2i'] 
     [1;46mCluster[0m 599 ['PEOE_VSA8'] 
     [1;46mCluster[0m 600 ['MATS1Z', 'MATS1m'] 
     [1;46mCluster[0m 601 ['ATSC1Z', 'ATSC1m'] 
     [1;46mCluster[0m 602 ['MATS2Z', 'MATS2m'] 
     [1;46mCluster[0m 603 ['ATSC2Z', 'ATSC2m'] 
     [1;46mCluster[0m 604 ['nP', 'NdsssP', 'SdsssP'] 
     [1;46mCluster[0m 605 ['BCUTs-1l'] 
     [1;46mCluster[0m 606 ['ATSC4pe', 'ATSC4are'] 
     [1;46mCluster[0m 607 ['ATSC4se'] 
     [1;46mCluster[0m 608 ['ATSC4dv'] 
     [1;46mCluster[0m 609 ['ATSC1s'] 
    
    Cluster info file [1mqdb1_F_s_hierarchical.cluster[0m file saved


     [1;46mCluster[0m 524 ['ATSC6v'] 
     [1;46mCluster[0m 525 ['ATSC6p'] 
     [1;46mCluster[0m 526 ['VR1_A', 'VR2_A'] 
     [1;46mCluster[0m 527 ['Xch-5d', 'Xch-5dv'] 
     [1;46mCluster[0m 528 ['Xch-4d', 'Xch-4dv'] 
     [1;46mCluster[0m 529 ['SMR_VSA4'] 
     [1;46mCluster[0m 530 ['PEOE_VSA10'] 
     [1;46mCluster[0m 531 ['NdssC'] 
     [1;46mCluster[0m 532 ['ATSC8Z', 'ATSC8m'] 
     [1;46mCluster[0m 533 ['ATSC8v'] 
     [1;46mCluster[0m 534 ['ATSC8p'] 
     [1;46mCluster[0m 535 ['ATSC8i'] 
     [1;46mCluster[0m 536 ['ATSC7d'] 
     [1;46mCluster[0m 537 ['EState_VSA3'] 
     [1;46mCluster[0m 538 ['EState_VSA6'] 
     [1;46mCluster[0m 539 ['SMR_VSA9'] 
     [1;46mCluster[0m 540 ['SlogP_VSA11'] 
     [1;46mCluster[0m 541 ['C3SP2'] 
     [1;46mCluster[0m 542 ['SaasC'] 
     [1;46mCluster[0m 543 ['EState_VSA7'] 
     [1;46mCluster[0m 544 ['PEOE_VSA6'] 
     [1;46mCluster[0m 545 ['AATS1dv'] 
     [1;46mCluster[0m 546 ['AATS2dv'] 
     [1;46mCluster[0m 547 ['AATS0dv'] 
     [1;46mCluster[0m 548 ['SIC0'] 
     [1;46mCluster[0m 549 ['BIC0'] 
     [1;46mCluster[0m 550 ['CIC1'] 
     [1;46mCluster[0m 551 ['AXp-0dv'] 
     [1;46mCluster[0m 552 ['AXp-1dv'] 
     [1;46mCluster[0m 553 ['MATS1s'] 
     [1;46mCluster[0m 554 ['GATS1s'] 
     [1;46mCluster[0m 555 ['GATS1dv'] 
     [1;46mCluster[0m 556 ['MATS1Z', 'MATS1m'] 
     [1;46mCluster[0m 557 ['ATSC1Z', 'ATSC1m'] 
     [1;46mCluster[0m 558 ['AATSC1dv'] 
     [1;46mCluster[0m 559 ['MATS1dv'] 
     [1;46mCluster[0m 560 ['AATSC1d'] 
     [1;46mCluster[0m 561 ['RPCG'] 
     [1;46mCluster[0m 562 ['BCUTd-1h'] 
     [1;46mCluster[0m 563 ['ETA_eta_B'] 
     [1;46mCluster[0m 564 ['nHetero', 'MID_h'] 
     [1;46mCluster[0m 565 ['SM1_Dzse', 'SM1_Dzpe'] 
     [1;46mCluster[0m 566 ['SM1_Dzare'] 
     [1;46mCluster[0m 567 ['SM1_Dzi'] 
     [1;46mCluster[0m 568 ['BCUTc-1h'] 
     [1;46mCluster[0m 569 ['nBondsD'] 
     [1;46mCluster[0m 570 ['TopoPSA'] 
     [1;46mCluster[0m 571 ['ATSC5pe', 'ATSC5are'] 
     [1;46mCluster[0m 572 ['ATSC5se'] 
     [1;46mCluster[0m 573 ['ATSC5dv'] 
     [1;46mCluster[0m 574 ['ATSC2se', 'ATSC2pe', 'ATSC2are'] 
     [1;46mCluster[0m 575 ['AATSC2dv', 'MATS2dv'] 
     [1;46mCluster[0m 576 ['ATSC2dv'] 
     [1;46mCluster[0m 577 ['AATSC2pe', 'AATSC2are'] 
     [1;46mCluster[0m 578 ['AATSC2se'] 
     [1;46mCluster[0m 579 ['MATS2pe', 'MATS2are'] 
     [1;46mCluster[0m 580 ['MATS2se'] 
     [1;46mCluster[0m 581 ['MATS2s'] 
     [1;46mCluster[0m 582 ['SlogP_VSA4'] 
     [1;46mCluster[0m 583 ['EState_VSA1'] 
     [1;46mCluster[0m 584 ['JGI5'] 
     [1;46mCluster[0m 585 ['SM1_Dzv', 'SM1_Dzp'] 
     [1;46mCluster[0m 586 ['TopoPSA(NO)'] 
     [1;46mCluster[0m 587 ['ATSC0s'] 
     [1;46mCluster[0m 588 ['EState_VSA10'] 
     [1;46mCluster[0m 589 ['NdO', 'SdO'] 
     [1;46mCluster[0m 590 ['VSA_EState2'] 
     [1;46mCluster[0m 591 ['nO', 'MID_O'] 
     [1;46mCluster[0m 592 ['SMR_VSA1'] 
     [1;46mCluster[0m 593 ['nHBAcc'] 
     [1;46mCluster[0m 594 ['ATS0s'] 
     [1;46mCluster[0m 595 ['BCUTv-1l'] 
     [1;46mCluster[0m 596 ['BCUTi-1h'] 
     [1;46mCluster[0m 597 ['AATSC0dv'] 
     [1;46mCluster[0m 598 ['BCUTdv-1h'] 
     [1;46mCluster[0m 599 ['BCUTse-1h'] 
     [1;46mCluster[0m 600 ['BCUTpe-1h'] 
     [1;46mCluster[0m 601 ['BCUTare-1h'] 
     [1;46mCluster[0m 602 ['BCUTp-1l'] 
     [1;46mCluster[0m 603 ['BCUTs-1h'] 
    
    Cluster info file [1mqdb1_F_s_hierarchical.cluster[0m file saved



```python
# numbef of clusters within the distance range
clust.cluster_dist()
```


    
![png](output_50_0.png)
    



    interactive(children=(ToggleButtons(description='Bins Index', options=(0, 3, 4, 5, 6, 7, 8, 9), value=0), Outp…



```python

```

<h1> 3. Feature Selection </h1>
<h2>  3.1. Monte Carlo </h2>

    selection(X_train, y_train, clust_info, model, pop_info, learning, bank, component, pntinterval)
    X_train : Descriptor Data
    y_train : end point
    clust_info : clustering informatino after clustering
    model : model algorithm
        Possible : 'PLS', 'MLR'
    pop_info
        None = start with empty
        population = used previous data continuously
    learning : # learning
    bank : # bank
    component : # descriptor
    pntiniterval : print value interval
    
    * return
        - select : selected feature
        - population : save point to restart


```python
from pyqsar import data_tools as dt
from pyqsar import model_tools as mt
from pyqsar import draw_mol
import pandas as pd
import numpy as np
```

### 3.1-1 Monte Carlo Feature Selection


```python
# split train data into EP and Descriptors 
X_train, y_train = mt.split_xy('.train')
```

    [1mPath :/home/jun/pyqsar/pyqsar3r25[0m
    Index   File Name
       0    qdb1_F_s.train
       1    qdb1_F_s_B.train


    Enter .train file index to load :  1


    
    



```python
# Feature Selection with MC + MLR. 
descriptor = 5
model_algo='MLR'   # MLR or PLS 
N_learning = 500
N_bank = 200
N_component = descriptor
select, population = mt.selection_mc(X_train,y_train,model=model_algo,
                                     #pop_info=None,               
                                     pop_info=population,            
                                     learning=N_learning,
                                     bank=N_bank,
                                     component=N_component,
                                     pntinterval=20 )
```

    Start time :  10:54:22
    Index   File Name
    0       qdb1_F_s_B_kmeans.cluster
    1       qdb1_F_s_B_som.cluster
    2       qdb1_F_s_kmeans.cluster
    3       qdb1_F_s_som.cluster
    4       qdb1_F_s_hierarchical.cluster
    
    


    Enter Clust Info File Index :  1


    [48;5;226mqdb1_F_s_B_som.cluster[0m    file selected
    
    [1;42mMLR[0m
                           R^2                  RMSE



      0%|          | 0/500 [00:00<?, ?it/s]


         20 => 10:54:32 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
         40 => 10:54:35 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
         60 => 10:54:39 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
         80 => 10:54:43 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        100 => 10:54:47 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        120 => 10:54:50 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        140 => 10:54:54 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        160 => 10:54:58 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        180 => 10:55:01 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        200 => 10:55:05 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        220 => 10:55:09 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        240 => 10:55:13 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        260 => 10:55:16 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        280 => 10:55:20 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        300 => 10:55:24 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        320 => 10:55:27 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        340 => 10:55:31 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        360 => 10:55:35 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        380 => 10:55:39 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        400 => 10:55:42 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        420 => 10:55:46 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        440 => 10:55:50 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        460 => 10:55:54 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        480 => 10:55:57 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
        500 => 10:56:01 [ 0.7613, 0.8539] ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
    R2       :  0.7613
    RMSE     :  0.8539
    Cluster  : ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']
    Model's cluster info [131, 70, 250, 216, 283]
    qdb1_F_s_B_mc_MLR_som.log  is updated!


##### Epoch Graph


```python
# Draw epoch vs. metrics
mt.Draw_epoch()
```

    Index   File Name
    0       qdb1_F_s_B_mc_MLR_som.log
    
    


    Enter Clust Info File Index :  0


    [48;5;226mqdb1_F_s_B_mc_MLR_som.log[0m    file selected
    



    
![png](output_58_3.png)
    


#### 3.1.2. Training with whole train set


```python
model = mt.GetModel(select,model_algo,n_component=descriptor)
```

    [1mPath :/home/jun/pyqsar/pyqsar3r25[0m
    Index   File Name
       0    qdb1_F_s.train
       1    qdb1_F_s_B.train


    Enter .train file index to load :  1


    
    



```python
# training with whole train set 
model.train_plot()
```


<div class="bk-root">
        <a href="https://bokeh.org" target="_blank" class="bk-logo bk-logo-small bk-logo-notebook"></a>
        <span id="1108">Loading BokehJS ...</span>
    </div>






<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>GATS1i</th>
      <th>NddC</th>
      <th>SMR_VSA9</th>
      <th>SpAbs_A</th>
      <th>Xp-0dv</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Coef Value</th>
      <td>0.650119</td>
      <td>-0.307476</td>
      <td>-0.37946</td>
      <td>0.800424</td>
      <td>-1.744295</td>
    </tr>
  </tbody>
</table>
</div>




<div class="bk-root" id="3fb6bcc7-c966-409b-bc32-80baeb2abb1f" data-root-id="1111"></div>





<h3> 3.1.3. Train with cross validation </h3>


```python
# Cross Validation 
model.k_fold() # default : k=5, run=1000
```

    sklearn R^2CV mean: 0.762429
    sklearn Q^2CV mean: 0.752753
    RMSE CV : 0.873862
    Features set = ['GATS1i', 'NddC', 'SMR_VSA9', 'SpAbs_A', 'Xp-0dv']




<div class="bk-root" id="d36be0cc-fb57-4345-8e1c-efe52e458665" data-root-id="1229"></div>





#### 3.1.4. Analysis


```python
# Values for selected descriptor, EP, predicted  values, Error 
model.features_table()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>GATS1i</th>
      <th>NddC</th>
      <th>SMR_VSA9</th>
      <th>SpAbs_A</th>
      <th>Xp-0dv</th>
      <th>EP</th>
      <th>Predict</th>
      <th>Error</th>
    </tr>
    <tr>
      <th>ID</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>11859</th>
      <td>-0.995487</td>
      <td>-0.171398</td>
      <td>0.602304</td>
      <td>-0.257506</td>
      <td>-0.091129</td>
      <td>-5.46</td>
      <td>-5.472863</td>
      <td>0.012863</td>
    </tr>
    <tr>
      <th>11899</th>
      <td>-0.972881</td>
      <td>-0.171398</td>
      <td>-0.575593</td>
      <td>-0.112802</td>
      <td>-0.154874</td>
      <td>-4.66</td>
      <td>-4.784189</td>
      <td>0.124189</td>
    </tr>
    <tr>
      <th>2723704</th>
      <td>0.197154</td>
      <td>-0.171398</td>
      <td>-0.575593</td>
      <td>-1.157113</td>
      <td>-1.029226</td>
      <td>-3.98</td>
      <td>-3.334291</td>
      <td>-0.645709</td>
    </tr>
    <tr>
      <th>3032338</th>
      <td>0.249884</td>
      <td>-0.171398</td>
      <td>-0.575593</td>
      <td>-1.026855</td>
      <td>-0.822031</td>
      <td>-4.00</td>
      <td>-3.557157</td>
      <td>-0.442843</td>
    </tr>
    <tr>
      <th>2346</th>
      <td>-0.581662</td>
      <td>4.833418</td>
      <td>-0.575593</td>
      <td>-0.076781</td>
      <td>-0.307311</td>
      <td>-6.54</td>
      <td>-5.773984</td>
      <td>-0.766016</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>6348</th>
      <td>0.643948</td>
      <td>4.833418</td>
      <td>-0.575593</td>
      <td>-1.493278</td>
      <td>-1.279054</td>
      <td>-4.56</td>
      <td>-4.415984</td>
      <td>-0.144016</td>
    </tr>
    <tr>
      <th>2723790</th>
      <td>-0.985418</td>
      <td>-0.171398</td>
      <td>-0.575593</td>
      <td>-1.404159</td>
      <td>-1.299579</td>
      <td>-3.84</td>
      <td>-3.829269</td>
      <td>-0.010731</td>
    </tr>
    <tr>
      <th>2723949</th>
      <td>-0.812529</td>
      <td>-0.171398</td>
      <td>-0.575593</td>
      <td>-1.404159</td>
      <td>-1.175735</td>
      <td>-3.64</td>
      <td>-3.932891</td>
      <td>0.292891</td>
    </tr>
    <tr>
      <th>969491</th>
      <td>-0.539767</td>
      <td>-0.171398</td>
      <td>-0.575593</td>
      <td>1.643416</td>
      <td>1.924606</td>
      <td>-6.28</td>
      <td>-6.724121</td>
      <td>0.444121</td>
    </tr>
    <tr>
      <th>6228</th>
      <td>3.392620</td>
      <td>-0.171398</td>
      <td>-0.575593</td>
      <td>-1.157113</td>
      <td>-1.137432</td>
      <td>-0.70</td>
      <td>-1.068115</td>
      <td>0.368115</td>
    </tr>
  </tbody>
</table>
<p>194 rows × 8 columns</p>
</div>




```python
# Correlation between EP and Selected Descriptors 
model.feature_corr()
```


    
![png](output_66_0.png)
    





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>GATS1i</th>
      <th>NddC</th>
      <th>SMR_VSA9</th>
      <th>SpAbs_A</th>
      <th>Xp-0dv</th>
      <th>EP</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>GATS1i</th>
      <td>1.000000</td>
      <td>-0.049496</td>
      <td>-0.323854</td>
      <td>-0.333575</td>
      <td>-0.176288</td>
      <td>0.479926</td>
    </tr>
    <tr>
      <th>NddC</th>
      <td>-0.049496</td>
      <td>1.000000</td>
      <td>-0.112251</td>
      <td>-0.102673</td>
      <td>-0.113369</td>
      <td>-0.130381</td>
    </tr>
    <tr>
      <th>SMR_VSA9</th>
      <td>-0.323854</td>
      <td>-0.112251</td>
      <td>1.000000</td>
      <td>0.561949</td>
      <td>0.468914</td>
      <td>-0.536796</td>
    </tr>
    <tr>
      <th>SpAbs_A</th>
      <td>-0.333575</td>
      <td>-0.102673</td>
      <td>0.561949</td>
      <td>1.000000</td>
      <td>0.926396</td>
      <td>-0.716174</td>
    </tr>
    <tr>
      <th>Xp-0dv</th>
      <td>-0.176288</td>
      <td>-0.113369</td>
      <td>0.468914</td>
      <td>0.926396</td>
      <td>1.000000</td>
      <td>-0.746864</td>
    </tr>
    <tr>
      <th>EP</th>
      <td>0.479926</td>
      <td>-0.130381</td>
      <td>-0.536796</td>
      <td>-0.716174</td>
      <td>-0.746864</td>
      <td>1.000000</td>
    </tr>
  </tbody>
</table>
</div>




```python
# Draw 2D images of selected molecules 
# to get a molecular ID, place mouse point on the dot  
mol = draw_mol.DrawMols(ID=['16115','6129','10107','40585'])
mol.show()
```

    [1mPath :/home/jun/pyqsar/pyqsar3r25[0m
    Index   File Name
       0    qdb1.sdf


    Enter .sdf file index to load :  0


    
    





    
![png](output_67_3.png)
    




```python
# maximum common substructure of selected molecules
substr = mol.common_substr()
mol.show_substr(substr)
```




    
![png](output_68_0.png)
    




```python
# 3D Image
mol.show_3D()
```


    interactive(children=(Dropdown(description='ID', options=('16115', '6129', '10107', '40585'), value='16115'), …


### 3.1.5. Model Save


```python
model.save()
```

    [1mPath :/home/jun/pyqsar/pyqsar3r25[0m
    Index   File Name
       0    qdb1_F_standard.info


    Enter .info file index to load :  0


    
    


    Selected Algorithm in Feature Selection (MC/GA) :  MC
    Selected Model (MLR/PLS) :  MLR
    Enter model name to  save as
    (default) qdb1_F_s_B_MC_MLR.model
    - :  


    qdb1_F_s_B_MC_MLR.model file saved


## 4. Model Test 

#### ModelTest(test=True, scaled=True)
    * Model Test Part
    * If test is True, then show test R2, RMSE, and train/test plot
    * if test if False, Just show abut Train
    
    * Scaled False means test data is external, (not from split module) so if scaled=False, scale the test data by scale_obj above stetp


```python
from pyqsar import data_tools as dt
from pyqsar import model_tools as mt
from pyqsar import draw_mol
import pandas as pd
import numpy as np
```

## 4.1. Test on test set


```python
ModelTest = mt.ModelTest()
```

    [1mPath :/home/jun/pyqsar/pyqsar3r25[0m
    Index   File Name
       0    qdb1_F_s_B_MC_MLR.model


    Enter .model file index to load :  0


    
    



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>GATS1i</th>
      <th>NddC</th>
      <th>SMR_VSA9</th>
      <th>SpAbs_A</th>
      <th>Xp-0dv</th>
      <th>Description</th>
    </tr>
    <tr>
      <th>Row</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Scale</th>
      <td>0.277416</td>
      <td>0.878603</td>
      <td>21.954656</td>
      <td>0.294521</td>
      <td>0.066528</td>
      <td>Per feature relative scaling of the data to ac...</td>
    </tr>
    <tr>
      <th>Mean</th>
      <td>0.918409</td>
      <td>0.284247</td>
      <td>26.083008</td>
      <td>0.047945</td>
      <td>0.496183</td>
      <td>The mean value for each feature in the trainin...</td>
    </tr>
    <tr>
      <th>Var</th>
      <td>0.076960</td>
      <td>0.771944</td>
      <td>482.006927</td>
      <td>0.086742</td>
      <td>0.004426</td>
      <td>The variance for each feature in the training set</td>
    </tr>
    <tr>
      <th>Coef</th>
      <td>0.650119</td>
      <td>-0.307476</td>
      <td>-0.379460</td>
      <td>0.800424</td>
      <td>-1.744295</td>
      <td>Coef value of each feature</td>
    </tr>
    <tr>
      <th>Padel Index</th>
      <td>350.000000</td>
      <td>581.000000</td>
      <td>763.000000</td>
      <td>3.000000</td>
      <td>537.000000</td>
      <td>PaDEL descriptor index</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>GATS1i</th>
      <th>NddC</th>
      <th>SMR_VSA9</th>
      <th>SpAbs_A</th>
      <th>Xp-0dv</th>
      <th>EP</th>
      <th>Predict</th>
      <th>Error</th>
    </tr>
    <tr>
      <th>ID</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>11859</th>
      <td>-0.995487</td>
      <td>-0.171398</td>
      <td>0.602304</td>
      <td>-0.257506</td>
      <td>-0.091129</td>
      <td>-5.46</td>
      <td>-5.472863</td>
      <td>0.012863</td>
    </tr>
    <tr>
      <th>11899</th>
      <td>-0.972881</td>
      <td>-0.171398</td>
      <td>-0.575593</td>
      <td>-0.112802</td>
      <td>-0.154874</td>
      <td>-4.66</td>
      <td>-4.784189</td>
      <td>0.124189</td>
    </tr>
    <tr>
      <th>2723704</th>
      <td>0.197154</td>
      <td>-0.171398</td>
      <td>-0.575593</td>
      <td>-1.157113</td>
      <td>-1.029226</td>
      <td>-3.98</td>
      <td>-3.334291</td>
      <td>-0.645709</td>
    </tr>
    <tr>
      <th>3032338</th>
      <td>0.249884</td>
      <td>-0.171398</td>
      <td>-0.575593</td>
      <td>-1.026855</td>
      <td>-0.822031</td>
      <td>-4.00</td>
      <td>-3.557157</td>
      <td>-0.442843</td>
    </tr>
    <tr>
      <th>2346</th>
      <td>-0.581662</td>
      <td>4.833418</td>
      <td>-0.575593</td>
      <td>-0.076781</td>
      <td>-0.307311</td>
      <td>-6.54</td>
      <td>-5.773984</td>
      <td>-0.766016</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>6348</th>
      <td>0.643948</td>
      <td>4.833418</td>
      <td>-0.575593</td>
      <td>-1.493278</td>
      <td>-1.279054</td>
      <td>-4.56</td>
      <td>-4.415984</td>
      <td>-0.144016</td>
    </tr>
    <tr>
      <th>2723790</th>
      <td>-0.985418</td>
      <td>-0.171398</td>
      <td>-0.575593</td>
      <td>-1.404159</td>
      <td>-1.299579</td>
      <td>-3.84</td>
      <td>-3.829269</td>
      <td>-0.010731</td>
    </tr>
    <tr>
      <th>2723949</th>
      <td>-0.812529</td>
      <td>-0.171398</td>
      <td>-0.575593</td>
      <td>-1.404159</td>
      <td>-1.175735</td>
      <td>-3.64</td>
      <td>-3.932891</td>
      <td>0.292891</td>
    </tr>
    <tr>
      <th>969491</th>
      <td>-0.539767</td>
      <td>-0.171398</td>
      <td>-0.575593</td>
      <td>1.643416</td>
      <td>1.924606</td>
      <td>-6.28</td>
      <td>-6.724121</td>
      <td>0.444121</td>
    </tr>
    <tr>
      <th>6228</th>
      <td>3.392620</td>
      <td>-0.171398</td>
      <td>-0.575593</td>
      <td>-1.157113</td>
      <td>-1.137432</td>
      <td>-0.70</td>
      <td>-1.068115</td>
      <td>0.368115</td>
    </tr>
  </tbody>
</table>
<p>194 rows × 8 columns</p>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Value</th>
    </tr>
    <tr>
      <th>Row</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>R2_CV</th>
      <td>0.7624294186002666</td>
    </tr>
    <tr>
      <th>Q2_CV</th>
      <td>0.7527528997291147</td>
    </tr>
    <tr>
      <th>RMSE_CV</th>
      <td>0.8738619202147496</td>
    </tr>
    <tr>
      <th>Model_Algo</th>
      <td>MLR</td>
    </tr>
    <tr>
      <th>Preprocessing</th>
      <td>standard</td>
    </tr>
  </tbody>
</table>
</div>



```python
ModelTest.model_test(scaled=True)
```

    [1mPath :/home/jun/pyqsar/pyqsar3r25[0m
    Index   File Name
       0    qdb1_F_s.test
       1    qdb1_F_s_B.test


    Enter .test file index to load :  1


    
    




<div class="bk-root" id="8fb5248e-f1ba-4f54-bbe6-f2bec30e21e7" data-root-id="1385"></div>





    Test R2 : 0.562727
    Test RMSE : 1.177887



```python
# Draw 2D images of selected molecules 
# to get a molecular ID, place mouse point on the dot  
mol = draw_mol.DrawMols(ID=['61198','7847','9395'])
mol.show()
```

    [1mPath :/home/jun/pyqsar/pyqsar3r25[0m
    Index   File Name
       0    qdb1.sdf


    Enter .sdf file index to load :  0


    
    





    
![png](output_78_3.png)
    




```python
# maximum common substructure of selected molecules
substr = mol.common_substr()
mol.show_substr(substr)
```




    
![png](output_79_0.png)
    



## 4.2. Test on new molecules


```python
ModelTest.editor()
```






<style>.bk-root, .bk-root .bk:before, .bk-root .bk:after {
  font-family: var(--jp-ui-font-size1);
  font-size: var(--jp-ui-font-size1);
  color: var(--jp-ui-font-color1);
}
</style>





<div id='1551'>
  <div class="bk-root" id="3bbcfb37-7c17-42ba-8c59-321c59d01d78" data-root-id="1551"></div>
</div>
<script type="application/javascript">(function(root) {
  function embed_document(root) {
    var docs_json = {"8c58dbed-5aeb-4108-8f9f-7243b423d783":{"defs":[{"extends":null,"module":null,"name":"ReactiveHTML1","overrides":[],"properties":[]},{"extends":null,"module":null,"name":"FlexBox1","overrides":[],"properties":[{"default":"flex-start","kind":null,"name":"align_content"},{"default":"flex-start","kind":null,"name":"align_items"},{"default":"row","kind":null,"name":"flex_direction"},{"default":"wrap","kind":null,"name":"flex_wrap"},{"default":"flex-start","kind":null,"name":"justify_content"}]},{"extends":null,"module":null,"name":"GridStack1","overrides":[],"properties":[{"default":"warn","kind":null,"name":"mode"},{"default":null,"kind":null,"name":"ncols"},{"default":null,"kind":null,"name":"nrows"},{"default":true,"kind":null,"name":"allow_resize"},{"default":true,"kind":null,"name":"allow_drag"},{"default":[],"kind":null,"name":"state"}]},{"extends":null,"module":null,"name":"click1","overrides":[],"properties":[{"default":"","kind":null,"name":"terminal_output"},{"default":"","kind":null,"name":"debug_name"},{"default":0,"kind":null,"name":"clears"}]},{"extends":null,"module":null,"name":"NotificationAreaBase1","overrides":[],"properties":[{"default":"bottom-right","kind":null,"name":"position"},{"default":0,"kind":null,"name":"_clear"}]},{"extends":null,"module":null,"name":"NotificationArea1","overrides":[],"properties":[{"default":[],"kind":null,"name":"notifications"},{"default":"bottom-right","kind":null,"name":"position"},{"default":0,"kind":null,"name":"_clear"},{"default":[{"background":"#ffc107","icon":{"className":"fas fa-exclamation-triangle","color":"white","tagName":"i"},"type":"warning"},{"background":"#007bff","icon":{"className":"fas fa-info-circle","color":"white","tagName":"i"},"type":"info"}],"kind":null,"name":"types"}]},{"extends":null,"module":null,"name":"Notification","overrides":[],"properties":[{"default":null,"kind":null,"name":"background"},{"default":3000,"kind":null,"name":"duration"},{"default":null,"kind":null,"name":"icon"},{"default":"","kind":null,"name":"message"},{"default":null,"kind":null,"name":"notification_type"},{"default":false,"kind":null,"name":"_destroyed"}]},{"extends":null,"module":null,"name":"TemplateActions1","overrides":[],"properties":[{"default":0,"kind":null,"name":"open_modal"},{"default":0,"kind":null,"name":"close_modal"}]},{"extends":null,"module":null,"name":"MaterialTemplateActions1","overrides":[],"properties":[{"default":0,"kind":null,"name":"open_modal"},{"default":0,"kind":null,"name":"close_modal"}]}],"roots":{"references":[{"attributes":{"children":[{"id":"1555"},{"id":"1558"},{"id":"1559"},{"id":"1560"}],"margin":[0,0,0,0],"name":"Column00112","sizing_mode":"stretch_width"},"id":"1554","type":"Column"},{"attributes":{"height":30,"margin":[5,10,5,10],"max_length":5000,"placeholder":"Enter File Name; File Extension is .csv ","sizing_mode":"fixed","width":400},"id":"1559","type":"TextAreaInput"},{"attributes":{"children":[{"id":"1556"},{"id":"1557"}],"margin":[5,5,5,5],"name":"","sizing_mode":"stretch_width","width":300},"id":"1555","type":"Column"},{"attributes":{"disabled":true,"height":50,"margin":[5,10,5,10],"max_length":5000,"sizing_mode":"fixed","title":"SMILES","width":400},"id":"1557","type":"TextAreaInput"},{"attributes":{"button_type":"primary","icon":null,"label":"Evaluation","margin":[5,10,5,10],"sizing_mode":"stretch_width","subscribed_events":["button_click"]},"id":"1561","type":"Button"},{"attributes":{"margin":[5,10,5,10],"name":"","sizing_mode":"stretch_width","text":"<b></b>"},"id":"1556","type":"Div"},{"attributes":{"active":[0],"labels":["New File?"],"margin":[5,10,5,10],"sizing_mode":"stretch_width"},"id":"1558","type":"CheckboxGroup"},{"attributes":{"height":320,"margin":[5,10,5,10],"max_length":5000,"sizing_mode":"fixed","value":"Ready","width":400},"id":"1560","type":"TextAreaInput"},{"attributes":{"children":[{"id":"1552"},{"id":"1561"}],"margin":[0,0,0,0],"name":"Column00114","sizing_mode":"stretch_width"},"id":"1551","type":"Column"},{"attributes":{"format":"smiles","guicolor":"#c0c0c0","height":500,"jme":"Not Subscribed","margin":[5,10,5,10],"mol":"Not Subscribed","mol3000":"Not Subscribed","name":"","sdf":"Not Subscribed","sizing_mode":"stretch_width","smiles":"Not Subscribed"},"id":"1553","type":"panel_chemistry.bokeh_extensions.jsme_editor.JSMEEditor"},{"attributes":{"children":[{"id":"1553"},{"id":"1554"}],"margin":[0,0,0,0],"name":"Row00113","sizing_mode":"stretch_width"},"id":"1552","type":"Row"},{"attributes":{"client_comm_id":"98666e793bfe471a9114460cdb3ff05a","comm_id":"700eabd67a974f0e9e3084056a735a52","plot_id":"1551"},"id":"1562","type":"panel.models.comm_manager.CommManager"},{"attributes":{"reload":false},"id":"1563","type":"panel.models.location.Location"}],"root_ids":["1551","1562","1563"]},"title":"Bokeh Application","version":"2.4.3"}};
    var render_items = [{"docid":"8c58dbed-5aeb-4108-8f9f-7243b423d783","root_ids":["1551"],"roots":{"1551":"3bbcfb37-7c17-42ba-8c59-321c59d01d78"}}];
    root.Bokeh.embed.embed_items_notebook(docs_json, render_items);
    for (const render_item of render_items) {
      for (const root_id of render_item.root_ids) {
	const id_el = document.getElementById(root_id)
	if (id_el.children.length && (id_el.children[0].className === 'bk-root')) {
	  const root_el = id_el.children[0]
	  root_el.id = root_el.id + '-rendered'
	}
      }
    }
  }
  if (root.Bokeh !== undefined && root.Bokeh.Panel !== undefined) {
    embed_document(root);
  } else {
    var attempts = 0;
    var timer = setInterval(function(root) {
      if (root.Bokeh !== undefined && root.Bokeh.Panel !== undefined) {
        clearInterval(timer);
        embed_document(root);
      } else if (document.readyState == "complete") {
        attempts++;
        if (attempts > 200) {
          clearInterval(timer);
          console.log("Bokeh: ERROR: Unable to run BokehJS code because BokehJS library is missing");
        }
      }
    }, 25, root)
  }
})(window);</script>



```python
ModelTest.new_molecule('new.csv')
```

# Below is the exmaple for Genetic Algorithm 

# Sample for GA + MLR 

### Genetic Algorithm
    selection(X_train, y_train, clust_info, clustering, model, pop_info, n_pop_size, N_generation, component)
    X_train : Descriptor Data
    y_train : end point
    clust_info : clustering informatino after clustering
    clustering : method of clustering
        Possible : 'hierarchical', 'kmeans', 'som'
    model : model algorithm
        Possible : 'PLS', 'MLR'
    pop_info
        None = start with empty
        population = used previous data continuously
    N_generation : # learning
    n_pop_size : # bank
    component : # descriptor


```python
from pyqsar import data_tools as dt
from pyqsar import model_tools as mt
from pyqsar import draw_mol
import pandas as pd
import numpy as np
# Train 데이터를 Descriptor에 해당하는 부분과 EP에 해당하는 값으로 나누는 단계
X_train, y_train = mt.split_xy('.train')
```

    [1mPath :/home/jun/pyqsar/pyqsar3r25[0m
    Index   File Name
       0    qdb1_F_s.train
       1    qdb1_F_s_B.train


    Enter .train file index to load :  1


    
    



```python
# GA + PLS을 이용한 Feature Selection 예
descriptor = 5
model_algo='MLR'
N_learning = 200
N_bank = 200
N_component = descriptor
select, population = mt.selection_ga(X_train,y_train,model=model_algo,
                                     pop_info=None,
                                     #pop_info=population,
                                     n_pop_size=N_bank,
                                     N_generation=N_learning,
                                     component=N_component,
                                    )
```

    [1;42mMLR[0m
    Index   File Name
    0       qdb1_F_s_B_kmeans.cluster
    1       qdb1_F_s_B_som.cluster
    2       qdb1_F_s_kmeans.cluster
    3       qdb1_F_s_som.cluster
    4       qdb1_F_s_hierarchical.cluster
    
    


    Enter Clust Info File Index :  0


    [48;5;226mqdb1_F_s_B_kmeans.cluster[0m    file selected
    
    ---------------  Initialization  ---------------



      0%|          | 0/200 [00:00<?, ?it/s]


      1.023459   0.657182   0.896157% 
    Generation        RMSD        R2      Portion
             1 :   0.932336   0.715509   1.075042% 
             2 :   0.932336   0.715509   0.916662% 
             3 :   0.920281   0.722819   0.858753% 
             4 :   0.920281   0.722819   0.772628% 
             5 :   0.920281   0.722819   0.723299% 
             6 :   0.918280   0.724022   0.667369% 
             7 :   0.909156   0.729480   0.656658% 
             8 :   0.909156   0.729480   0.615554% 
             9 :   0.904750   0.732095   0.610864% 
            10 :   0.904750   0.732095   0.601702% 
            11 :   0.904750   0.732095   0.575212% 
            12 :   0.904750   0.732095   0.574346% 
            13 :   0.904750   0.732095   0.559688% 
            14 :   0.904750   0.732095   0.536423% 
            15 :   0.890763   0.740315   0.567041% 
            16 :   0.890763   0.740315   0.566254% 
            17 :   0.890763   0.740315   0.569545% 
            18 :   0.890763   0.740315   0.557126% 
            19 :   0.890763   0.740315   0.567638% 
            20 :   0.886450   0.742823   0.578561% 
            21 :   0.884638   0.743873   0.555912% 
            22 :   0.875943   0.748884   0.572390% 
            23 :   0.860984   0.757387   0.616812% 
            24 :   0.860984   0.757387   0.594476% 
            25 :   0.860984   0.757387   0.601345% 
            26 :   0.858075   0.759024   0.604617% 
            27 :   0.858075   0.759024   0.602086% 
            28 :   0.858075   0.759024   0.608674% 
            29 :   0.858075   0.759024   0.602948% 
            30 :   0.858075   0.759024   0.611481% 
            31 :   0.858075   0.759024   0.612383% 
            32 :   0.858075   0.759024   0.593617% 
            33 :   0.858075   0.759024   0.595281% 
            34 :   0.858075   0.759024   0.566684% 
            35 :   0.858075   0.759024   0.588956% 
            36 :   0.858075   0.759024   0.598461% 
            37 :   0.858075   0.759024   0.590095% 
            38 :   0.858075   0.759024   0.594964% 
            39 :   0.858075   0.759024   0.577448% 
            40 :   0.858075   0.759024   0.594622% 
            41 :   0.858075   0.759024   0.602755% 
            42 :   0.850381   0.763326   0.589734% 
            43 :   0.850381   0.763326   0.601404% 
            44 :   0.850381   0.763326   0.638427% 
            45 :   0.850381   0.763326   0.631942% 
            46 :   0.850381   0.763326   0.610559% 
            47 :   0.850381   0.763326   0.595082% 
            48 :   0.850215   0.763419   0.567399% 
            49 :   0.850215   0.763419   0.594726% 
            50 :   0.850215   0.763419   0.599074% 
            51 :   0.850215   0.763419   0.601234% 
            52 :   0.850215   0.763419   0.588719% 
            53 :   0.850215   0.763419   0.606841% 
            54 :   0.850215   0.763419   0.588262% 
            55 :   0.850215   0.763419   0.610245% 
            56 :   0.850215   0.763419   0.596253% 
            57 :   0.850215   0.763419   0.605206% 
            58 :   0.850215   0.763419   0.567937% 
            59 :   0.850215   0.763419   0.614141% 
            60 :   0.850215   0.763419   0.593419% 
            61 :   0.850215   0.763419   0.610180% 
            62 :   0.850215   0.763419   0.595413% 
            63 :   0.850215   0.763419   0.606065% 
            64 :   0.850215   0.763419   0.595287% 
            65 :   0.844810   0.766417   0.599449% 
            66 :   0.844810   0.766417   0.588064% 
            67 :   0.844810   0.766417   0.610841% 
            68 :   0.844810   0.766417   0.578432% 
            69 :   0.842796   0.767529   0.614661% 
            70 :   0.842796   0.767529   0.616812% 
            71 :   0.842796   0.767529   0.611121% 
            72 :   0.842796   0.767529   0.586034% 
            73 :   0.842796   0.767529   0.608772% 
            74 :   0.842796   0.767529   0.637771% 
            75 :   0.842796   0.767529   0.653855% 
            76 :   0.842796   0.767529   0.630697% 
            77 :   0.842796   0.767529   0.656531% 
            78 :   0.842796   0.767529   0.596226% 
            79 :   0.842796   0.767529   0.645798% 
            80 :   0.842796   0.767529   0.624628% 
            81 :   0.842796   0.767529   0.628077% 
            82 :   0.841839   0.768057   0.643363% 
            83 :   0.841839   0.768057   0.640988% 
            84 :   0.841839   0.768057   0.647340% 
            85 :   0.841839   0.768057   0.646072% 
            86 :   0.841839   0.768057   0.638042% 
            87 :   0.841839   0.768057   0.622135% 
            88 :   0.841839   0.768057   0.643509% 
            89 :   0.841839   0.768057   0.649068% 
            90 :   0.841839   0.768057   0.650728% 
            91 :   0.841839   0.768057   0.626636% 
            92 :   0.841839   0.768057   0.661544% 
            93 :   0.841839   0.768057   0.681305% 
            94 :   0.841839   0.768057   0.646065% 
            95 :   0.841839   0.768057   0.621469% 
            96 :   0.834891   0.771870   0.650749% 
            97 :   0.834891   0.771870   0.692281% 
            98 :   0.834891   0.771870   0.652003% 
            99 :   0.834891   0.771870   0.656738% 
           100 :   0.834891   0.771870   0.621883% 
           101 :   0.834891   0.771870   0.683390% 
           102 :   0.834891   0.771870   0.688958% 
           103 :   0.834891   0.771870   0.707486% 
           104 :   0.834891   0.771870   0.712590% 
           105 :   0.834891   0.771870   0.699436% 
           106 :   0.834891   0.771870   0.696882% 
           107 :   0.834891   0.771870   0.688858% 
           108 :   0.834891   0.771870   0.688719% 
           109 :   0.834891   0.771870   0.689320% 
           110 :   0.834891   0.771870   0.673938% 
           111 :   0.834891   0.771870   0.688842% 
           112 :   0.834891   0.771870   0.731186% 
           113 :   0.834891   0.771870   0.739667% 
           114 :   0.834891   0.771870   0.685925% 
           115 :   0.834891   0.771870   0.677661% 
           116 :   0.834891   0.771870   0.701339% 
           117 :   0.834891   0.771870   0.700543% 
           118 :   0.834891   0.771870   0.745563% 
           119 :   0.834891   0.771870   0.756046% 
           120 :   0.834891   0.771870   0.706712% 
           121 :   0.834891   0.771870   0.734583% 
           122 :   0.834891   0.771870   0.731948% 
           123 :   0.834891   0.771870   0.708713% 
           124 :   0.834891   0.771870   0.716997% 
           125 :   0.834891   0.771870   0.730688% 
           126 :   0.834891   0.771870   0.719982% 
           127 :   0.834891   0.771870   0.694916% 
           128 :   0.834891   0.771870   0.698458% 
           129 :   0.834891   0.771870   0.742783% 
           130 :   0.834891   0.771870   0.710483% 
           131 :   0.834891   0.771870   0.710101% 
           132 :   0.834891   0.771870   0.709370% 
           133 :   0.834891   0.771870   0.724714% 
           134 :   0.834891   0.771870   0.696862% 
           135 :   0.834891   0.771870   0.729289% 
           136 :   0.834891   0.771870   0.740347% 
           137 :   0.834891   0.771870   0.763258% 
           138 :   0.834891   0.771870   0.730993% 
           139 :   0.834891   0.771870   0.748140% 
           140 :   0.834891   0.771870   0.774732% 
           141 :   0.834891   0.771870   0.759414% 
           142 :   0.834891   0.771870   0.735120% 
           143 :   0.834891   0.771870   0.734554% 
           144 :   0.834891   0.771870   0.742146% 
           145 :   0.834891   0.771870   0.710556% 
           146 :   0.834891   0.771870   0.716994% 
           147 :   0.834891   0.771870   0.744108% 
           148 :   0.834891   0.771870   0.741567% 
           149 :   0.834891   0.771870   0.744557% 
           150 :   0.834891   0.771870   0.723289% 
           151 :   0.834891   0.771870   0.797403% 
           152 :   0.834891   0.771870   0.777516% 
           153 :   0.834891   0.771870   0.727969% 
           154 :   0.834891   0.771870   0.748080% 
           155 :   0.834891   0.771870   0.734471% 
           156 :   0.834891   0.771870   0.725796% 
           157 :   0.834891   0.771870   0.689813% 
           158 :   0.834891   0.771870   0.719328% 
           159 :   0.834891   0.771870   0.739740% 
           160 :   0.834891   0.771870   0.748382% 
           161 :   0.834891   0.771870   0.762900% 
           162 :   0.834891   0.771870   0.769227% 
           163 :   0.834891   0.771870   0.748044% 
           164 :   0.834891   0.771870   0.770813% 
           165 :   0.834891   0.771870   2.003274% 
           166 :   0.834891   0.771870   1.135646% 
           167 :   0.834891   0.771870   0.848887% 
           168 :   0.834891   0.771870   0.732654% 
           169 :   0.834891   0.771870   0.698804% 
           170 :   0.834891   0.771870   0.641889% 
           171 :   0.834891   0.771870   0.630723% 
           172 :   0.834891   0.771870   0.607518% 
           173 :   0.834891   0.771870   0.621363% 
           174 :   0.834891   0.771870   0.658822% 
           175 :   0.834891   0.771870   0.622790% 
           176 :   0.834891   0.771870   0.649694% 
           177 :   0.834891   0.771870   0.613736% 
           178 :   0.834891   0.771870   0.597820% 
           179 :   0.834891   0.771870   0.623882% 
           180 :   0.834891   0.771870   0.652260% 
           181 :   0.834891   0.771870   0.661051% 
           182 :   0.834891   0.771870   0.640270% 
           183 :   0.834891   0.771870   0.656483% 
           184 :   0.834891   0.771870   0.651284% 
           185 :   0.834891   0.771870   0.646368% 
           186 :   0.834891   0.771870   0.631881% 
           187 :   0.834891   0.771870   0.678613% 
           188 :   0.834891   0.771870   0.675008% 
           189 :   0.834891   0.771870   0.654332% 
           190 :   0.834891   0.771870   0.696828% 
           191 :   0.834891   0.771870   0.642341% 
           192 :   0.834891   0.771870   0.635671% 
           193 :   0.834891   0.771870   0.649271% 
           194 :   0.834891   0.771870   0.640641% 
           195 :   0.834891   0.771870   0.677181% 
           196 :   0.834891   0.771870   0.661407% 
           197 :   0.834891   0.771870   0.662011% 
           198 :   0.834891   0.771870   0.654082% 
           199 :   0.834891   0.771870   0.646612% 
           200 : ---------------  End of Generation  ---------------
    Fit Val   RMSE    R2      Portion 
    1.5886    0.8349  0.7719  0.647%  ['ATSC4p', 'MIC1', 'NddC', 'SlogP_VSA2', 'Xp-0dv']
    1.5354    0.8411  0.7685  0.625%  ['JGI9', 'MIC1', 'NddC', 'SlogP_VSA2', 'Xp-0dv']
    1.5344    0.8412  0.7684  0.625%  ['MIC1', 'NddC', 'SlogP_VSA2', 'VSA_EState9', 'Xp-0dv']
    qdb1_F_s_B_ga_MLR_kmeans.log  is saved!


#### Epoch Graph


```python
mt.Draw_epoch()
```

    Index   File Name
    0       qdb1_F_s_B_mc_MLR_som.log
    1       qdb1_F_s_B_ga_MLR_kmeans.log
    
    


    Enter Clust Info File Index :  1


    [48;5;226mqdb1_F_s_B_ga_MLR_kmeans.log[0m    file selected
    



    
![png](output_89_3.png)
    


#### Model Information and Save


```python
model = mt.GetModel(select,model_algo,descriptor)
```

    [1mPath :/home/jun/pyqsar/pyqsar3r25[0m
    Index   File Name
       0    qdb1_F_s.train
       1    qdb1_F_s_B.train


    Enter .train file index to load :  1


    
    



```python
model.train_plot()
```


<div class="bk-root">
        <a href="https://bokeh.org" target="_blank" class="bk-logo bk-logo-small bk-logo-notebook"></a>
        <span id="1582">Loading BokehJS ...</span>
    </div>






<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ATSC4p</th>
      <th>MIC1</th>
      <th>NddC</th>
      <th>SlogP_VSA2</th>
      <th>Xp-0dv</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Coef Value</th>
      <td>0.206766</td>
      <td>-0.349078</td>
      <td>-0.313306</td>
      <td>0.633987</td>
      <td>-1.280102</td>
    </tr>
  </tbody>
</table>
</div>




<div class="bk-root" id="1d840afc-be03-4065-ab77-c70d4b2864cb" data-root-id="1585"></div>






```python
model.k_fold()
```

    sklearn R^2CV mean: 0.773204
    sklearn Q^2CV mean: 0.759639
    RMSE CV : 0.85897
    Features set = ['ATSC4p', 'MIC1', 'NddC', 'SlogP_VSA2', 'Xp-0dv']




<div class="bk-root" id="4ed496e7-0a26-45b4-be6e-0e7aa3e1c28a" data-root-id="1741"></div>






```python
model.features_table()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ATSC4p</th>
      <th>MIC1</th>
      <th>NddC</th>
      <th>SlogP_VSA2</th>
      <th>Xp-0dv</th>
      <th>EP</th>
      <th>Predict</th>
      <th>Error</th>
    </tr>
    <tr>
      <th>ID</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>11859</th>
      <td>0.400891</td>
      <td>1.010308</td>
      <td>-0.171398</td>
      <td>-0.521548</td>
      <td>-0.091129</td>
      <td>-5.46</td>
      <td>-5.041991</td>
      <td>-0.418009</td>
    </tr>
    <tr>
      <th>11899</th>
      <td>0.160976</td>
      <td>1.472546</td>
      <td>-0.171398</td>
      <td>-0.537675</td>
      <td>-0.154874</td>
      <td>-4.66</td>
      <td>-5.181579</td>
      <td>0.521579</td>
    </tr>
    <tr>
      <th>2723704</th>
      <td>-0.430709</td>
      <td>-0.026408</td>
      <td>-0.171398</td>
      <td>0.099304</td>
      <td>-1.029226</td>
      <td>-3.98</td>
      <td>-3.257571</td>
      <td>-0.722429</td>
    </tr>
    <tr>
      <th>3032338</th>
      <td>0.105913</td>
      <td>-0.062641</td>
      <td>-0.171398</td>
      <td>0.055038</td>
      <td>-0.822031</td>
      <td>-4.00</td>
      <td>-3.427263</td>
      <td>-0.572737</td>
    </tr>
    <tr>
      <th>2346</th>
      <td>0.271452</td>
      <td>-0.094217</td>
      <td>4.833418</td>
      <td>-0.516730</td>
      <td>-0.307311</td>
      <td>-6.54</td>
      <td>-5.971437</td>
      <td>-0.568563</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>6348</th>
      <td>0.375887</td>
      <td>-0.982260</td>
      <td>4.833418</td>
      <td>-0.591293</td>
      <td>-1.279054</td>
      <td>-4.56</td>
      <td>-4.443189</td>
      <td>-0.116811</td>
    </tr>
    <tr>
      <th>2723790</th>
      <td>0.716717</td>
      <td>-0.407856</td>
      <td>-0.171398</td>
      <td>-0.521028</td>
      <td>-1.299579</td>
      <td>-3.84</td>
      <td>-2.934370</td>
      <td>-0.905630</td>
    </tr>
    <tr>
      <th>2723949</th>
      <td>0.899870</td>
      <td>-0.222806</td>
      <td>-0.171398</td>
      <td>-0.531936</td>
      <td>-1.175735</td>
      <td>-3.64</td>
      <td>-3.126547</td>
      <td>-0.513453</td>
    </tr>
    <tr>
      <th>969491</th>
      <td>3.190576</td>
      <td>2.116495</td>
      <td>-0.171398</td>
      <td>1.342976</td>
      <td>1.924606</td>
      <td>-6.28</td>
      <td>-6.249588</td>
      <td>-0.030412</td>
    </tr>
    <tr>
      <th>6228</th>
      <td>0.902360</td>
      <td>-1.017566</td>
      <td>-0.171398</td>
      <td>1.265142</td>
      <td>-1.137432</td>
      <td>-0.70</td>
      <td>-1.758305</td>
      <td>1.058305</td>
    </tr>
  </tbody>
</table>
<p>194 rows × 8 columns</p>
</div>




```python
model.feature_corr()
```


    
![png](output_95_0.png)
    





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ATSC4p</th>
      <th>MIC1</th>
      <th>NddC</th>
      <th>SlogP_VSA2</th>
      <th>Xp-0dv</th>
      <th>EP</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ATSC4p</th>
      <td>1.000000</td>
      <td>-0.164075</td>
      <td>0.064593</td>
      <td>-0.207991</td>
      <td>-0.077051</td>
      <td>0.128154</td>
    </tr>
    <tr>
      <th>MIC1</th>
      <td>-0.164075</td>
      <td>1.000000</td>
      <td>-0.007424</td>
      <td>0.162557</td>
      <td>0.540639</td>
      <td>-0.563949</td>
    </tr>
    <tr>
      <th>NddC</th>
      <td>0.064593</td>
      <td>-0.007424</td>
      <td>1.000000</td>
      <td>-0.045657</td>
      <td>-0.113369</td>
      <td>-0.130381</td>
    </tr>
    <tr>
      <th>SlogP_VSA2</th>
      <td>-0.207991</td>
      <td>0.162557</td>
      <td>-0.045657</td>
      <td>1.000000</td>
      <td>0.284031</td>
      <td>0.120575</td>
    </tr>
    <tr>
      <th>Xp-0dv</th>
      <td>-0.077051</td>
      <td>0.540639</td>
      <td>-0.113369</td>
      <td>0.284031</td>
      <td>1.000000</td>
      <td>-0.746864</td>
    </tr>
    <tr>
      <th>EP</th>
      <td>0.128154</td>
      <td>-0.563949</td>
      <td>-0.130381</td>
      <td>0.120575</td>
      <td>-0.746864</td>
      <td>1.000000</td>
    </tr>
  </tbody>
</table>
</div>




```python
model.save()
```

    [1mPath :/home/jun/pyqsar/pyqsar3r25[0m
    Index   File Name
       0    qdb1_F_standard.info


    Enter .info file index to load :  0


    
    


    Selected Algorithm in Feature Selection (MC/GA) :  MC
    Selected Model (MLR/PLS) :  MLR
    Enter model name to  save as
    (default) qdb1_F_s_B_MC_MLR.model
    - :  


    qdb1_F_s_B_MC_MLR.model file saved



```python
mol = draw_mol.DrawMols(ID=['16115','6129'])
mol.show()
```

    [1mPath :/home/jun/pyqsar/pyqsar3r25[0m
    Index   File Name
       0    qdb1.sdf


    Enter .sdf file index to load :  0


    
    





    
![png](output_97_3.png)
    




```python
substr = mol.common_substr()
mol.show_substr(substr)
```




    
![png](output_98_0.png)
    




```python
mol.show_3D()
```


    interactive(children=(Dropdown(description='ID', options=('16115', '6129'), value='16115'), Dropdown(descripti…


## Model Information
-----------------------------------------------


```python
ModelTest = mt.ModelTest()
```

    [1mPath :/home/jun/pyqsar/pyqsar3r25[0m
    Index   File Name
       0    qdb1_F_s_B_MC_MLR.model


    Enter .model file index to load :  0


    
    



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ATSC4p</th>
      <th>MIC1</th>
      <th>NddC</th>
      <th>SlogP_VSA2</th>
      <th>Xp-0dv</th>
      <th>Description</th>
    </tr>
    <tr>
      <th>Row</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Scale</th>
      <td>3.065640</td>
      <td>7.410414</td>
      <td>0.878603</td>
      <td>3.948729</td>
      <td>0.066528</td>
      <td>Per feature relative scaling of the data to ac...</td>
    </tr>
    <tr>
      <th>Mean</th>
      <td>-0.977846</td>
      <td>19.200235</td>
      <td>0.284247</td>
      <td>2.481973</td>
      <td>0.496183</td>
      <td>The mean value for each feature in the trainin...</td>
    </tr>
    <tr>
      <th>Var</th>
      <td>9.398148</td>
      <td>54.914237</td>
      <td>0.771944</td>
      <td>15.592459</td>
      <td>0.004426</td>
      <td>The variance for each feature in the training set</td>
    </tr>
    <tr>
      <th>Coef</th>
      <td>0.206766</td>
      <td>-0.349078</td>
      <td>-0.313306</td>
      <td>0.633987</td>
      <td>-1.280102</td>
      <td>Coef value of each feature</td>
    </tr>
    <tr>
      <th>Padel Index</th>
      <td>258.000000</td>
      <td>728.000000</td>
      <td>581.000000</td>
      <td>765.000000</td>
      <td>537.000000</td>
      <td>PaDEL descriptor index</td>
    </tr>
  </tbody>
</table>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ATSC4p</th>
      <th>MIC1</th>
      <th>NddC</th>
      <th>SlogP_VSA2</th>
      <th>Xp-0dv</th>
      <th>EP</th>
      <th>Predict</th>
      <th>Error</th>
    </tr>
    <tr>
      <th>ID</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>11859</th>
      <td>0.400891</td>
      <td>1.010308</td>
      <td>-0.171398</td>
      <td>-0.521548</td>
      <td>-0.091129</td>
      <td>-5.46</td>
      <td>-5.041991</td>
      <td>-0.418009</td>
    </tr>
    <tr>
      <th>11899</th>
      <td>0.160976</td>
      <td>1.472546</td>
      <td>-0.171398</td>
      <td>-0.537675</td>
      <td>-0.154874</td>
      <td>-4.66</td>
      <td>-5.181579</td>
      <td>0.521579</td>
    </tr>
    <tr>
      <th>2723704</th>
      <td>-0.430709</td>
      <td>-0.026408</td>
      <td>-0.171398</td>
      <td>0.099304</td>
      <td>-1.029226</td>
      <td>-3.98</td>
      <td>-3.257571</td>
      <td>-0.722429</td>
    </tr>
    <tr>
      <th>3032338</th>
      <td>0.105913</td>
      <td>-0.062641</td>
      <td>-0.171398</td>
      <td>0.055038</td>
      <td>-0.822031</td>
      <td>-4.00</td>
      <td>-3.427263</td>
      <td>-0.572737</td>
    </tr>
    <tr>
      <th>2346</th>
      <td>0.271452</td>
      <td>-0.094217</td>
      <td>4.833418</td>
      <td>-0.516730</td>
      <td>-0.307311</td>
      <td>-6.54</td>
      <td>-5.971437</td>
      <td>-0.568563</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>6348</th>
      <td>0.375887</td>
      <td>-0.982260</td>
      <td>4.833418</td>
      <td>-0.591293</td>
      <td>-1.279054</td>
      <td>-4.56</td>
      <td>-4.443189</td>
      <td>-0.116811</td>
    </tr>
    <tr>
      <th>2723790</th>
      <td>0.716717</td>
      <td>-0.407856</td>
      <td>-0.171398</td>
      <td>-0.521028</td>
      <td>-1.299579</td>
      <td>-3.84</td>
      <td>-2.934370</td>
      <td>-0.905630</td>
    </tr>
    <tr>
      <th>2723949</th>
      <td>0.899870</td>
      <td>-0.222806</td>
      <td>-0.171398</td>
      <td>-0.531936</td>
      <td>-1.175735</td>
      <td>-3.64</td>
      <td>-3.126547</td>
      <td>-0.513453</td>
    </tr>
    <tr>
      <th>969491</th>
      <td>3.190576</td>
      <td>2.116495</td>
      <td>-0.171398</td>
      <td>1.342976</td>
      <td>1.924606</td>
      <td>-6.28</td>
      <td>-6.249588</td>
      <td>-0.030412</td>
    </tr>
    <tr>
      <th>6228</th>
      <td>0.902360</td>
      <td>-1.017566</td>
      <td>-0.171398</td>
      <td>1.265142</td>
      <td>-1.137432</td>
      <td>-0.70</td>
      <td>-1.758305</td>
      <td>1.058305</td>
    </tr>
  </tbody>
</table>
<p>194 rows × 8 columns</p>
</div>



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Value</th>
    </tr>
    <tr>
      <th>Row</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>R2_CV</th>
      <td>0.7732035864945177</td>
    </tr>
    <tr>
      <th>Q2_CV</th>
      <td>0.759638991396595</td>
    </tr>
    <tr>
      <th>RMSE_CV</th>
      <td>0.858969974417255</td>
    </tr>
    <tr>
      <th>Model_Algo</th>
      <td>MLR</td>
    </tr>
    <tr>
      <th>Preprocessing</th>
      <td>standard</td>
    </tr>
  </tbody>
</table>
</div>



```python
ModelTest.model_test(scaled=True)
```

    [1mPath :/home/jun/pyqsar/pyqsar3r25[0m
    Index   File Name
       0    qdb1_F_s.test
       1    qdb1_F_s_B.test


    Enter .test file index to load :  1


    
    




<div class="bk-root" id="3cfe5354-895c-4661-bb31-4bac7f3ba729" data-root-id="2142"></div>





    Test R2 : 0.598894
    Test RMSE : 1.123991



```python
ModelTest.editor()
```






<style>.bk-root, .bk-root .bk:before, .bk-root .bk:after {
  font-family: var(--jp-ui-font-size1);
  font-size: var(--jp-ui-font-size1);
  color: var(--jp-ui-font-color1);
}
</style>





<div id='2359'>
  <div class="bk-root" id="dc42ebcf-cac2-4e18-94fe-5670c6014d14" data-root-id="2359"></div>
</div>
<script type="application/javascript">(function(root) {
  function embed_document(root) {
    var docs_json = {"4b9f0a9f-af56-4a1e-bce3-a85be840d880":{"defs":[{"extends":null,"module":null,"name":"ReactiveHTML1","overrides":[],"properties":[]},{"extends":null,"module":null,"name":"FlexBox1","overrides":[],"properties":[{"default":"flex-start","kind":null,"name":"align_content"},{"default":"flex-start","kind":null,"name":"align_items"},{"default":"row","kind":null,"name":"flex_direction"},{"default":"wrap","kind":null,"name":"flex_wrap"},{"default":"flex-start","kind":null,"name":"justify_content"}]},{"extends":null,"module":null,"name":"GridStack1","overrides":[],"properties":[{"default":"warn","kind":null,"name":"mode"},{"default":null,"kind":null,"name":"ncols"},{"default":null,"kind":null,"name":"nrows"},{"default":true,"kind":null,"name":"allow_resize"},{"default":true,"kind":null,"name":"allow_drag"},{"default":[],"kind":null,"name":"state"}]},{"extends":null,"module":null,"name":"click1","overrides":[],"properties":[{"default":"","kind":null,"name":"terminal_output"},{"default":"","kind":null,"name":"debug_name"},{"default":0,"kind":null,"name":"clears"}]},{"extends":null,"module":null,"name":"NotificationAreaBase1","overrides":[],"properties":[{"default":"bottom-right","kind":null,"name":"position"},{"default":0,"kind":null,"name":"_clear"}]},{"extends":null,"module":null,"name":"NotificationArea1","overrides":[],"properties":[{"default":[],"kind":null,"name":"notifications"},{"default":"bottom-right","kind":null,"name":"position"},{"default":0,"kind":null,"name":"_clear"},{"default":[{"background":"#ffc107","icon":{"className":"fas fa-exclamation-triangle","color":"white","tagName":"i"},"type":"warning"},{"background":"#007bff","icon":{"className":"fas fa-info-circle","color":"white","tagName":"i"},"type":"info"}],"kind":null,"name":"types"}]},{"extends":null,"module":null,"name":"Notification","overrides":[],"properties":[{"default":null,"kind":null,"name":"background"},{"default":3000,"kind":null,"name":"duration"},{"default":null,"kind":null,"name":"icon"},{"default":"","kind":null,"name":"message"},{"default":null,"kind":null,"name":"notification_type"},{"default":false,"kind":null,"name":"_destroyed"}]},{"extends":null,"module":null,"name":"TemplateActions1","overrides":[],"properties":[{"default":0,"kind":null,"name":"open_modal"},{"default":0,"kind":null,"name":"close_modal"}]},{"extends":null,"module":null,"name":"MaterialTemplateActions1","overrides":[],"properties":[{"default":0,"kind":null,"name":"open_modal"},{"default":0,"kind":null,"name":"close_modal"}]}],"roots":{"references":[{"attributes":{"active":[0],"labels":["New File?"],"margin":[5,10,5,10],"sizing_mode":"stretch_width"},"id":"2366","type":"CheckboxGroup"},{"attributes":{"children":[{"id":"2363"},{"id":"2366"},{"id":"2367"},{"id":"2368"}],"margin":[0,0,0,0],"name":"Column00130","sizing_mode":"stretch_width"},"id":"2362","type":"Column"},{"attributes":{"format":"smiles","guicolor":"#c0c0c0","height":500,"jme":"Not Subscribed","margin":[5,10,5,10],"mol":"Not Subscribed","mol3000":"Not Subscribed","name":"","sdf":"Not Subscribed","sizing_mode":"stretch_width","smiles":"Not Subscribed"},"id":"2361","type":"panel_chemistry.bokeh_extensions.jsme_editor.JSMEEditor"},{"attributes":{"children":[{"id":"2360"},{"id":"2369"}],"margin":[0,0,0,0],"name":"Column00132","sizing_mode":"stretch_width"},"id":"2359","type":"Column"},{"attributes":{"height":30,"margin":[5,10,5,10],"max_length":5000,"placeholder":"Enter File Name; File Extension is .csv ","sizing_mode":"fixed","width":400},"id":"2367","type":"TextAreaInput"},{"attributes":{"client_comm_id":"17bfff7910b143f7813e40086c48c32f","comm_id":"c82d20ff19d94eadb379ca5c491bbe7c","plot_id":"2359"},"id":"2370","type":"panel.models.comm_manager.CommManager"},{"attributes":{"height":320,"margin":[5,10,5,10],"max_length":5000,"sizing_mode":"fixed","value":"Ready","width":400},"id":"2368","type":"TextAreaInput"},{"attributes":{"children":[{"id":"2361"},{"id":"2362"}],"margin":[0,0,0,0],"name":"Row00131","sizing_mode":"stretch_width"},"id":"2360","type":"Row"},{"attributes":{"margin":[5,10,5,10],"name":"","sizing_mode":"stretch_width","text":"<b></b>"},"id":"2364","type":"Div"},{"attributes":{"button_type":"primary","icon":null,"label":"Evaluation","margin":[5,10,5,10],"sizing_mode":"stretch_width","subscribed_events":["button_click"]},"id":"2369","type":"Button"},{"attributes":{"children":[{"id":"2364"},{"id":"2365"}],"margin":[5,5,5,5],"name":"","sizing_mode":"stretch_width","width":300},"id":"2363","type":"Column"},{"attributes":{"disabled":true,"height":50,"margin":[5,10,5,10],"max_length":5000,"sizing_mode":"fixed","title":"SMILES","width":400},"id":"2365","type":"TextAreaInput"}],"root_ids":["2359","2370"]},"title":"Bokeh Application","version":"2.4.3"}};
    var render_items = [{"docid":"4b9f0a9f-af56-4a1e-bce3-a85be840d880","root_ids":["2359"],"roots":{"2359":"dc42ebcf-cac2-4e18-94fe-5670c6014d14"}}];
    root.Bokeh.embed.embed_items_notebook(docs_json, render_items);
    for (const render_item of render_items) {
      for (const root_id of render_item.root_ids) {
	const id_el = document.getElementById(root_id)
	if (id_el.children.length && (id_el.children[0].className === 'bk-root')) {
	  const root_el = id_el.children[0]
	  root_el.id = root_el.id + '-rendered'
	}
      }
    }
  }
  if (root.Bokeh !== undefined && root.Bokeh.Panel !== undefined) {
    embed_document(root);
  } else {
    var attempts = 0;
    var timer = setInterval(function(root) {
      if (root.Bokeh !== undefined && root.Bokeh.Panel !== undefined) {
        clearInterval(timer);
        embed_document(root);
      } else if (document.readyState == "complete") {
        attempts++;
        if (attempts > 200) {
          clearInterval(timer);
          console.log("Bokeh: ERROR: Unable to run BokehJS code because BokehJS library is missing");
        }
      }
    }, 25, root)
  }
})(window);</script>



```python
ModelTest.new_molecule('asdf.csv')
```


    ---------------------------------------------------------------------------

    FileNotFoundError                         Traceback (most recent call last)

    /tmp/ipykernel_17625/2019858943.py in <module>
    ----> 1 ModelTest.new_molecule('asdf.csv')
    

    ~/pyqsar/pyqsar3r25/pyqsar/model_tools.py in new_molecule(self, file_name)
       2069 
       2070     def new_molecule(self,file_name):
    -> 2071         data = pd.read_csv(file_name)
       2072         min_val = np.min(self.Feature['EP'].values)
       2073         work_path = os.getcwd().replace(sys.path[0]+'/','')


    ~/miniforge3/envs/pq3/lib/python3.9/site-packages/pandas/util/_decorators.py in wrapper(*args, **kwargs)
        209                 else:
        210                     kwargs[new_arg_name] = new_arg_value
    --> 211             return func(*args, **kwargs)
        212 
        213         return cast(F, wrapper)


    ~/miniforge3/envs/pq3/lib/python3.9/site-packages/pandas/util/_decorators.py in wrapper(*args, **kwargs)
        329                     stacklevel=find_stack_level(),
        330                 )
    --> 331             return func(*args, **kwargs)
        332 
        333         # error: "Callable[[VarArg(Any), KwArg(Any)], Any]" has no


    ~/miniforge3/envs/pq3/lib/python3.9/site-packages/pandas/io/parsers/readers.py in read_csv(filepath_or_buffer, sep, delimiter, header, names, index_col, usecols, squeeze, prefix, mangle_dupe_cols, dtype, engine, converters, true_values, false_values, skipinitialspace, skiprows, skipfooter, nrows, na_values, keep_default_na, na_filter, verbose, skip_blank_lines, parse_dates, infer_datetime_format, keep_date_col, date_parser, dayfirst, cache_dates, iterator, chunksize, compression, thousands, decimal, lineterminator, quotechar, quoting, doublequote, escapechar, comment, encoding, encoding_errors, dialect, error_bad_lines, warn_bad_lines, on_bad_lines, delim_whitespace, low_memory, memory_map, float_precision, storage_options)
        948     kwds.update(kwds_defaults)
        949 
    --> 950     return _read(filepath_or_buffer, kwds)
        951 
        952 


    ~/miniforge3/envs/pq3/lib/python3.9/site-packages/pandas/io/parsers/readers.py in _read(filepath_or_buffer, kwds)
        603 
        604     # Create the parser.
    --> 605     parser = TextFileReader(filepath_or_buffer, **kwds)
        606 
        607     if chunksize or iterator:


    ~/miniforge3/envs/pq3/lib/python3.9/site-packages/pandas/io/parsers/readers.py in __init__(self, f, engine, **kwds)
       1440 
       1441         self.handles: IOHandles | None = None
    -> 1442         self._engine = self._make_engine(f, self.engine)
       1443 
       1444     def close(self) -> None:


    ~/miniforge3/envs/pq3/lib/python3.9/site-packages/pandas/io/parsers/readers.py in _make_engine(self, f, engine)
       1733                 if "b" not in mode:
       1734                     mode += "b"
    -> 1735             self.handles = get_handle(
       1736                 f,
       1737                 mode,


    ~/miniforge3/envs/pq3/lib/python3.9/site-packages/pandas/io/common.py in get_handle(path_or_buf, mode, encoding, compression, memory_map, is_text, errors, storage_options)
        854         if ioargs.encoding and "b" not in ioargs.mode:
        855             # Encoding
    --> 856             handle = open(
        857                 handle,
        858                 ioargs.mode,


    FileNotFoundError: [Errno 2] No such file or directory: 'asdf.csv'



```python

```
