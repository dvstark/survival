'''survival analysis routines...mostly R wrappers'''

import rpy2.robjects as ro
import numpy as np
import fortranformat as ff
from rpy2.robjects.packages import importr

def setup():

    # import standard packages
    utils = importr('utils')
    base = importr('base')      

    # import those packages needed for the survival analysis codes
    packages = ['survival', 'NADA', 'stats']
    for p in packages:
        utils.install_packages(p)
    

def rfitter():

    ro.r('tbl <- read.table("/Users/dstark/python_packages/my_packages/survival/rdata.dat", head=T, fill=T)')
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
    upper limits. This is essentially a wrapper for rfitter'''
    
    tup = np.column_stack((ylim,xvar,yvar))
    header_line = ff.FortranRecordWriter(fmt)
    line = ff.FortranRecordWriter(fmt)

    pfile=open('/Users/dstark/python_packages/my_packages/survival/rdata.dat','w')
    pfile.write('lim xvar gs\n')
    for t in tup:
        line.write(t)
        pfile.write(line.write(t)+'\n')
    pfile.close()
    out = rfitter()
    
    #make this into a dictionary
    fit = {'slope':out[0],'intercept':out[1],'tau':out[2],'p':out[3],'sigma_out':out[4],'sigma_down':out[5]}
    
    return fit

def kaplan_meier(var,ylim):
    
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
