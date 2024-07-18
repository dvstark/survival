'''survival analysis routines...mostly R wrappers'''

import os

import numpy as np
import rpy2.robjects as ro
import fortranformat as ff
from rpy2.robjects.packages import importr

def setup():

    '''Run this routine once to install the R packages needed for the
    routines in this package. You may get a warning message about a library
    not being writable, and an option to use a private library. Go ahead and say 
    "yes" and this should all work'''

    # import standard packages
    utils = importr('utils')
    base = importr('base')      

    # import those packages needed for the survival analysis codes
    packages = ['survival', 'NADA', 'stats']
    for p in packages:
        utils.install_packages(p)
    

def run_ats_r(datafile):

    '''
    note: rename this to something 
    '''

    ro.r('tbl <- read.table("{}", head=T, fill=T)'.format(datafile))
    ro.r('dim(tbl) ; names(tbl) ; summary(tbl)')
              
    #read in various variables        
    ro.r('cen_gs <- seq(FALSE, FALSE, length=dim(tbl)[1])')
    ro.r('cen_gs[tbl$lim==1] <- TRUE')
    
    ro.r('library(survival)')
    ro.r('library(NADA)')
    ro.r('library(stats)')
    
    ro.r('cenken_out <- cenken(tbl$gs, cen_gs, tbl$xvar)')
    ro.r('print(cenken_out)')
    
    ro.r('yfit = tbl$xvar*cenken_out$slope') #set y-intercept to zero for now
    ro.r('resid = tbl$gs - yfit') #residuals
    
    #based on examples, I think I need to switch residuals to be right-censored
    ro.r('Xstatus <- seq(1,1,length.out=length(resid))')  # 1 = detected
    ro.r('Xstatus[cen_gs==TRUE] <- 0')  # 0 = right-censored
    ro.r('survobj <- Surv(-resid, Xstatus)')
    ro.r('KM.XLF<- survfit(survobj ~ 1, conf.int=0.68, conf.type="plain",conf.lower="modified")')
    ro.r('KM.output <- summary(survfit(survobj~1))')
    
    #interpolate to gte y intercept
    ro.r('int=approx(KM.XLF$surv,-(KM.XLF$time),0.5)') #interpolate to find where surv = 0.5
    
    #also get the 1-sigma scatter
    ro.r('slow=approx(KM.XLF$surv,-(KM.XLF$time),0.16)') #interpolate to find where surv = 0.5
    ro.r('shigh=approx(KM.XLF$surv,-(KM.XLF$time),0.84)') #interpolate to find where surv = 0.5
    
    par = np.asarray([ro.r('cenken_out$slope'),ro.r('int$y'), ro.r('cenken_out$tau'), ro.r('cenken_out$p'), ro.r('int$y-slow$y'), ro.r('shigh$y - int$y')])
    return par



def ats_fit(xvar,yvar,ylim,fmt='(I4, 2F10.3)'):
    
    '''Runs cenken from R (Alritas-Sen-Theil estimator) to determine
    correlation strength and line fit parameters in the presence of 
    upper limits. This is essentially a wrapper for run_ats_r
    
    Input:
    xvar : array-like
        The "x" variable
    yvar: array-like
        The "y" variable
    ylim: bool
        Flag indicating whether the corresponding y value is an upper limit
    fmt: str
        The format of the output data file that's sent to R. Just don't mess
        with this unless you know what you're doing.

    Output:
    
    out: dict
        A dictionary containing the results of the fit:
            slope = slope of fitted line
            int = intercept of fitted line
            tau = correlation coefficient
            p = statistical significance of correlation (e.g. < 0.05 would be
            significant at 2-sigma)
            sigma_up = scatter above the fitted line
            sigma_down = scatter below the fitted line
    
    '''
    
    tup = np.column_stack((ylim,xvar,yvar))
    header_line = ff.FortranRecordWriter(fmt)
    line = ff.FortranRecordWriter(fmt)

    datafile = os.getcwd() + 'rdata.dat'

    pfile=open(datafile,'w')
    pfile.write('lim xvar gs\n')
    for t in tup:
        line.write(t)
        pfile.write(line.write(t)+'\n')
    pfile.close()
    out = run_ats_r(datafile)
    
    #make this into a dictionary
    fit = {'slope':out[0],'intercept':out[1],'tau':out[2],'p':out[3],'sigma_up':out[4],'sigma_down':out[5]}
    
    return fit

def kaplan_meier(var,ylim):
    
    '''Runs the kaplan meier estimator to determine the CDF of the input data
    
    Input:

    var : array-like
        variable of interest
    ylim: array-like bool
        boolean array indicating if the corresponding y value is an upper limit
    
    Output:
    
    out: dict
        A dictionary containing:
            midpoint: where the CDF hits 0.5
            l68: lower 1-sigma confidence interval on midpoint
            h68: upper 1-sigma confidence interval on midpoint
            km: another dictionary with the results:
                x: the y variable value
                surv: the resulting cd

    
    '''
    



    tup = np.column_stack((ylim,var))
    header_line = ff.FortranRecordWriter('(I4, F10.3)')
    line = ff.FortranRecordWriter('(I4, F10.3)')

    pfile=open('km.dat','w')
    pfile.write('lim gs\n')
    for t in tup:
        line.write(t)
        pfile.write(line.write(t)+'\n')
    pfile.close()
    
    ro.r('tbl <- read.table("km.dat", head=T, fill=T)')
    ro.r('cen_gs <- seq(FALSE, FALSE, length=dim(tbl)[1])')
    ro.r('cen_gs[tbl$lim==1] <- TRUE')
    
    ro.r('print(tbl)')
    
    ro.r('library(survival)')
    ro.r('library(NADA)')
    ro.r('library(stats)')
    
    #based on examples, I think I need to switch residuals to be right-censored
    ro.r('Xstatus <- seq(1,1,length.out=length(tbl$gs))')  # 1 = detected
    ro.r('Xstatus[cen_gs==TRUE] <- 0')  # 0 = right-censored
    ro.r('survobj <- Surv(-tbl$gs, Xstatus)')
    ro.r('KM.XLF<- survfit(survobj ~ 1, conf.int=0.68, conf.type="plain",conf.lower="modified")')
    ro.r('KM.output <- summary(survfit(survobj~1))')
    
    #interpolate to gte y intercept
    ro.r('int=approx(KM.XLF$surv,-(KM.XLF$time),0.5)') #interpolate to find where surv = 0.5
    
    #also get the 1-sigma scatter
    ro.r('slow=approx(KM.XLF$surv,-(KM.XLF$time),0.16)') #interpolate to find where surv = 0.5
    ro.r('shigh=approx(KM.XLF$surv,-(KM.XLF$time),0.84)') #interpolate to find where surv = 0.5
 
    km={'surv':ro.r('KM.XLF$surv'),'x':ro.r('-(KM.XLF$time)')}
    par = np.asarray([ro.r('int$y'), ro.r('int$y-slow$y'), ro.r('shigh$y - int$y')])
    dict={'midpoint':ro.r('int$y'),'l68':ro.r('int$y-slow$y'),'h68':ro.r('shigh$y - int$y'),'km':km}
    return dict





    #thetaplot = radontrace['theta']*180/np.pi
    #rhoplot = radontrace['rho']
    #sel=thetaplot > 180
    #thetaplot[sel] = thetaplot[sel] - 180
    #rhoplot[sel] = -1*rhoplot[sel]
    #sel=thetaplot < 0
    #thetaplot[sel] = thetaplot[sel] + 180
    #rhoplot[sel] = -1*rhoplot[sel]
    
    #rhoplot[thetaplot<180]=-rhoplot[thetaplot<180]
