import pandas as pd
import numpy as np
from scipy.optimize import leastsq,least_squares


class PlaneFit():

    def __init__(self,filepath,locnamecol, xcol,ycol,zcol,tcol,clusters):

        self.filepath = filepath
        self.locnamecol = locnamecol
        self.xcol = xcol
        self.ycol = ycol
        self.zcol = zcol
        self.tcol = tcol
        self.clusters = clusters
        self.data = pd.read_excel(self.filepath)


    def planeFit3Point(self,xs,ys,zs):  
    
        x1,x2,x3 = xs
        y1,y2,y3 = ys
        z1,z2,z3 = zs
        
        # get twp different points in the plane. 
        v1 = np.array([x1-x2,y1-y2,z1-z2])
        v2 = np.array([x3-x2,y3-y2,z3-z2])
        
        # Compute the cross product of the two obtained vectors This is the normal vector of the plane
        solve = np.cross(v1,v2)
        dp = solve[0] *x1 + solve[1]*y1 + solve[2]*z1
        d = dp * -1
        a,b,c = solve
        return a,b,c,d

    def planeFitLeastSquares(self,xs,ys,zs):
    
        XYZ = np.vstack([xs,ys,zs])
        # Inital guess of the plane
        p0 = [0.005, 0.005, 0.005, 0.05]

        def f_min(X,p):
            plane_xyz = p[0:3]
            distance = (plane_xyz*X.T).sum(axis=1) + p[3]
            return distance / np.linalg.norm(plane_xyz)

        def residuals(params, signal, X):
            return f_min(X, params)

        sol = least_squares(residuals, p0, args=(None, XYZ))
        fit = sol.x
        residual = sol.fun
        se = residual.var()
        sez = np.array(zs).var()
        r2 = 1 - se/sez
        #r = (f_min(XYZ, sol)**2).sum()
        # calculate the R2 value
        
        return fit, r2


    def gradient(self,fit):
        a,b,c,d = fit
        return float(np.sqrt(((a**2 + b**2)/c**2)))

    def direction(self,fit):
        
        a,b,c,d = fit    
        theta = np.arctan(b/a) * 180/np.pi

        if a > 0:
            angle = 90 - theta
        elif a < 0:
            angle = 270 - theta
        elif a == 0:
            if b > 0:
                angle = 0
            elif b < 0:
                angle = 180
            else:
                angle = None
        if c < 0:
            angle = angle + 180
            
        return float(np.mod(angle,360))
    
    def analyzeCluster(self,data,wells,wellCol,timeCol,xCol,yCol,zCol):
    
        dfsub = data.loc[data[wellCol].isin(wells)]
        times = pd.unique(dfsub[timeCol])
        fits = []
        grads = []
        directions = []
        r2s = []
        ftimes = []
        
        for time in times:
            dataTime = dfsub.loc[dfsub[timeCol] == time]
            x = dataTime[xCol].values.tolist()
            y = dataTime[yCol].values.tolist()
            z = dataTime[zCol].values.tolist()
            
            if len(x) < 3:
                pass
            
            elif len(x) == 3:
                fit = self.planeFit3Point(x,y,z)
                r2 = 1
                
            else:
                fit,r2 = self.planeFitLeastSquares(x,y,z)
            
            r2s.append(r2)
            fits.append(fit)
            grads.append(self.gradient(fit))
            directions.append(self.direction(fit))
            ftimes.append(time)
        
        clusterdf = pd.DataFrame(data = [ftimes,fits,grads,directions,r2s])
        clusterdf = clusterdf.T
        clusterdf.columns = ['DateTime','Plane Fit Coefficients','Gradient','Direction','R2']
        return clusterdf

    def analyzeAllClusters(self):
        
        results = []
        for c in self.clusters:
            r = self.analyzeCluster(self.data,c,self.locnamecol,self.tcol,self.xcol,self.ycol,self.zcol)
            results.append(r)

        resultsdf = pd.DataFrame(data=[self.clusters,results])
        resultsdf = resultsdf.T
        resultsdf.columns = ['Well Cluster','Results DataFrame']

        return resultsdf