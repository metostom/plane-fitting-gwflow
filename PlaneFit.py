import pandas as pd
import numpy as np
from scipy.optimize import leastsq,least_squares


class PlaneFit():

    """
    PlaneFit Data Class

    Attributes
    ----------
    locnamecol : str
        location name column header
    xcol : str
        x-coordinate column header
    ycol : str
        y-coordinate column header
    zcol : str
        z-coordinate column header
    tcol : str
        time stamp column header   
    clusters : list
        list of cluster lists to run the plane fitting on
    """

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
        '''
        Solves the equation of a plane (Ax + By + Cz - D = 0) that intersects the 3 points

            Parameters:
                    xs (array) : array of x-coordinates for the points:
                    ys (array) : array of y-coordinates for the points:
                    zs (array) : array of z-coordinates for the points:
                                
            Returns:
                    fit (array) : array of the plane coefficients (A,B,C,D)
        '''
    
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

        '''
        Optimizes the equation of a plane (Ax + By + Cz - D = 0) that intersects the given points

            Parameters:
                    xs (array) : array of x-coordinates for the points:
                    ys (array) : array of y-coordinates for the points:
                    zs (array) : array of z-coordinates for the points:
                                
            Returns:
                    fit (array) : array of the plane coefficients (A,B,C,D)
                    r2 (float) : sum of residuals squared - goodness of fit metric
        '''
    
        XYZ = np.vstack([xs,ys,zs])
        # Inital guess of the plane
        p0 = [0.005, 0.005, 0.005, 0.05]

        # minimizing function
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
        
        return fit, r2


    def gradient(self,fit):
        '''
        Calculates the gradient of a plane given the coefficients of the plane equation (Ax + By + Cz - D = 0)

            Parameters:
                    fit (array) : array of the plane coefficients (A,B,C,D)
                                
            Returns:
                    gradient (float) : gradient of the plane
        '''
        a,b,c,d = fit
        gradient = float(np.sqrt(((a**2 + b**2)/c**2)))

        return gradient

    def direction(self,fit):
        '''
        Calculates the direction of maximum gradient of a plane given the coefficients of the plane equation (Ax + By + Cz - D = 0) 
        and converts the angle to compass degrees. 

            Parameters:
                    fit (array) : array of the plane coefficients (A,B,C,D)
                                
            Returns:
                    direction (float) : compass direction of maximum gradient, in degrees from north
        '''
        
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
        
        # modulus to 360 degrees
        direction = float(np.mod(angle,360))

        return direction
    
    def analyzeCluster(self,data,wells,wellCol,timeCol,xCol,yCol,zCol):
        '''
        Computes the direction, gradient, R2, and fit for each time the well cluster has data 

            Parameters:
                    data (dataframe) : df of all x,y,z,t data
                    wells (list) : cluster to analyze
                        locnamecol : str
                    wellCol : str
                        location name column header
                    xCol : str
                        x-coordinate column header
                    yCol : str
                        y-coordinate column header
                    zCol : str
                        z-coordinate column header
                    timeCol : str
                        time stamp column header  

                                
            Returns:
                    clusterdf (dataframe) : DataFrame of results for the cluster of the form:
                        DateTime | Plane Fit Coefficients | Gradient | Direction | R2
        '''
    
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
        '''
        Analyzes all well clusters
                                
            Returns:
                    Results df (dataframe) : DataFrame of results for each cluster of the form:
                        Well Cluster | Results DataFrame
                    
                    Where Results DataFrame is of the form:
                    DateTime | Plane Fit Coefficients | Gradient | Direction | R2
        '''

        
        results = []
        for c in self.clusters:
            r = self.analyzeCluster(self.data,c,self.locnamecol,self.tcol,self.xcol,self.ycol,self.zcol)
            results.append(r)

        resultsdf = pd.DataFrame(data=[self.clusters,results])
        resultsdf = resultsdf.T
        resultsdf.columns = ['Well Cluster','Results DataFrame']

        return resultsdf