# plane-fitting-gwflow

##  Overview
Computation of groundwater flow gradient and direction using N-point plane fitting solution. Computes true fit for the 3-point problem and a least-squares approximation for 4+ point solutions. 


Points are used to solve or estimate plane equation of the form:

<img src="https://latex.codecogs.com/svg.latex?Ax+By+Cz-D=0" title="formula" size=100px/>

the gradient and direction are then calculated as below.

<b> Gradient </b>


<img src="https://latex.codecogs.com/svg.latex?gradient=\sqrt{\frac{A^2+B^2}{C^2}}" title="formula" size=100px/>

<b> Direction </b>

<img src="https://latex.codecogs.com/svg.latex?\theta=arctan\left(\frac{B}{A}\right)" title="formula2" size=100px/>

where theta is degrees from the x-axis, then converted to compass degrees depending on the quadrant

## Uses
<b> Python </b>
- Pandas
- Numpy
- Scipy

## Usage
```
from PlaneFit import PlaneFit

clusters = [
    ["MW-4","MW-13","MW-15"],
    ["MW-12","MW-13","MW-15"],
    ["MW-4","MW-6","MW-15"],
    ["MW-6","MW-12","MW-15"],
]


cd = PlaneFit(
        filepath = 'WellElevations.xlsx',
        locnamecol = 'Location_ID'
        xcol = 'X',
        ycol = 'Y',
        zcol = "Water Elevation",
        tcol = "DateTime",
        clusters = clusters
)

resultsdf = cd.analyzeAllClusters()

```

## Graphed Histogram of GW flow Direction Output Looks like

![Alt Text](hists.png)
