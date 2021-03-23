import sys
import os
import os.path
import numpy as np
from astropy.io import fits as pf
include_path='/Users/simon/common/python/include/'
sys.path.append(include_path)
import scipy.optimize as op
#from multiprocessing import Pool
#from iminuit import Minuit


import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import time
#from time import gmtime, strftime
t_i = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())



    
def pass_model(Mpass,OptimMpass):
    global M
    global OptimM
    M=Mpass
    OptimM=OptimMpass
    


def lnlike(theta):

    nvar=len(theta)
    names = list(map( (lambda x: x[0]),OptimM.domain))
    for iparam in range(nvar):
        # print("in lnlike settting attributes", names[iparam],theta[iparam])
        setattr(M,names[iparam],theta[iparam])


    chi2=M.polar_expansions()
        
    statusstring=''
        
    for iparam in range(nvar):
        statusstring=statusstring+names[iparam]+" "+str(theta[iparam])+" "

    statusstring=statusstring+" -> "+str(chi2)

    if (OptimM.PrintOptimStatus):
        print( statusstring)
    
    return -0.5*chi2



def lnprior(theta):
    inside=1
    bnds = list(map( (lambda x: x[1]),OptimM.domain))
    for iparam in list(range(len(theta))):
        if (bnds[iparam][0] < theta[iparam] < bnds[iparam][1]):
            inside *=1
        else:
            inside *=0
    if (inside): 
        return 0.0
    else:
        return -np.inf

    
#def lnprob(theta, bnds):
#    lp = lnprior(theta,bnds)
#    if not np.isfinite(lp):
#        return -np.inf
#    return lp + lnlike(theta)

def lnprob(theta):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta)


def run_scipy_optimize_minimize(M,OptimM,x,bnds):
    pass_model(M,OptimM)
    print( "starting op.minimize")
    start_time=time.time()
    nll = lambda *args: -lnlike(*args)
    print( "domain: ",OptimM.domain)
    ftol=1E-8 # 1e-10 too small leads to abnormal termination
    #eps=0.1*np.ones(len(x))
    eps=0.1*np.ones(len(x))
    for iparam in list(range(len(x))):
        fullrange=(bnds[iparam][1]-bnds[iparam][0])
        eps[iparam]=fullrange*1E-2
    print("step sizes:",eps)

    result = op.minimize(nll, x, tol=ftol,bounds=bnds,options={'eps':eps})
    print( "result",result)
    result_ml  = result["x"]
    print( "Optim done in (elapsed time):",   time.time()-start_time)
    print( "computing errors with Hessian")
    tmp_i = np.zeros(len(result_ml))
    errors_ml= np.zeros(len(result_ml))
    for i in list(range(len(result_ml))):
        tmp_i[i] = 1.0
        uncertainty_i = np.sqrt(result.hess_inv(tmp_i)[i])
        errors_ml[i]=uncertainty_i
        tmp_i[i] = 0.0
        print(('{0:12.4e} +- {1:.1e}'.format(result.x[i], uncertainty_i)))
    return (result_ml,errors_ml)




def exec_ConjGrad(M,OptimM):
    
    names = list(map( (lambda x: x[0]),OptimM.domain))
    bnds = list(map( (lambda x: x[1]),OptimM.domain))
    nvar=len(list(names))
    sample_theta=list(range(nvar))
    for iparam in list(range(nvar)):
        sample_theta[iparam]=getattr(M,names[iparam])

    x = np.array( sample_theta  ) 

    M.Verbose=False
    M.DumpAllFitsFiles=False
    M.Grid=True
    M.PlotAzimuthalProfile=False
    M.PlotRadialProfile=False
    M.XCheckInv=False

    M.prep_files()

    print("Init ConjGrad with params:")
    for iparam in list(range(nvar)):
        print("x",x[iparam],"bnds",bnds[iparam])

    (result_ml,errors_ml)=run_scipy_optimize_minimize(M,OptimM,x,bnds)
            
    np.save(M.workdir+'result_ml.dat',result_ml)
    np.save(M.workdir+'result_ml_errors.dat',errors_ml)


    statusstring=''
    for iparam in range(nvar):
        setattr(M,names[iparam],result_ml[iparam])
        statusstring=statusstring+names[iparam]+" %.3f " %(result_ml[iparam])

    print("Finished conjgrad and set M object to: "+statusstring)


    M.XCheckInv=True
    M.PlotRadialProfile=True
    M.PlotAzimuthalProfile=True
    M.DumpAllFitsFiles=True
    M.Verbose=True
    M.Grid=False

    print("running final polar expansion")
    M.polar_expansions()


    return  result_ml




def exec_emcee(M,result_ml,RunMCMC,OptimM):
    Nit=OptimM.Nit
    nwalkers=OptimM.nwalkers
    n_cores=OptimM.n_cores_MCMC
    burn_in=OptimM.burn_in #100

    M.DumpAllFitsFiles=False
    M.Verbose=False
    M.Grid=True
    M.PlotAzimuthalProfile=False
    M.PlotRadialProfile=False
    M.XCheckInv=False


    OptimM.PrintOptimStatus=False

    pass_model(M,OptimM)

    workdir=M.workdir
    names = list(map( (lambda x: x[0]),OptimM.domain))
    bnds = list(map( (lambda x: x[1]),OptimM.domain))

    ranges = list(map( (lambda x: x[1][1]-x[1][0]),OptimM.domain))

    allowed_ranges=np.array(ranges)
    print("allowed_ranges ",allowed_ranges)

    
    nvar = len(names)
    print( "mcmc with nvar=",nvar)
    
    ndim =nvar
    #ndim, nwalkers = nvar, 60
    #pos = [result_ml + 1e-1*np.random.randn(ndim) for i in list(range(nwalkers))]
    pos=[]
    for i in list(range(nwalkers)):
        if (np.any(allowed_ranges < 0.)):
            sys.exit("wrong order of bounds in domains")
        awalkerinit=result_ml+(1e-3*np.random.randn(ndim)*allowed_ranges)
        pos.append(awalkerinit)

    print("init for emcee :", result_ml)
    #print("init ball for emcee :", pos)
    
    os.environ["OMP_NUM_THREADS"] = "1"
    import emcee
    #nit=3000
    print( "in exec_emcee with RunMCMC=",RunMCMC)
    if RunMCMC:
        print( bnds)
        print( "now about to call run_mcmc with Nit",Nit,"and nmwalkers",nwalkers," and ncores",n_cores)
        #sampler = emcee.ensemblesampler(nwalkers, ndim, lnprob, args=(bnds))

        from multiprocessing import Pool
        with Pool(n_cores) as pool:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool)
            start = time.time()
            sampler.run_mcmc(pos, Nit, progress=True)
            end = time.time()
            multi_time = end - start
            print("Multiprocessing took {0:.1f} seconds".format(multi_time))
        
        print( "************ finish ***************")
        samples = sampler.chain  # chain= array(nwalkers,nit,ndim)
        lnprobs = sampler.lnprobability
        
        ######### save samples
        np.save(workdir+'samples.dat',samples)
        np.save(workdir+'lnprobs.dat',lnprobs)
        # end time
        t_f = time.strftime("%y-%m-%d %h:%m:%s", time.gmtime())
        print( "t_i = "+str(t_i))
        print( "t_f = "+str(t_f))
        
        print(("mean acceptance fraction: {0:.3f} "  .format(np.mean(sampler.acceptance_fraction))))
        f=open(workdir+'acceptance.dat', 'w')
        f.write(str(t_i)+' \n')
        f.write(str(t_f)+' \n')
        f.write("Nit = "+str(Nit)+' \n')
        f.write("nwalkers = "+str(nwalkers)+' \n')
        f.write("ndim = "+str(ndim)+' \n')
        f.write("mean acceptance fraction: {0:.3f}"  .format(np.mean(sampler.acceptance_fraction)) +' \n')
        f.close() 
        
        #autocorr=sampler.get_autocorr_time(c=1, low=1)
        #print( "autocorr\n",autocorr  )
        
    else:
        samples=np.load(workdir+'samples.dat.npy')
        lnprobs=np.load(workdir+'lnprobs.dat.npy')
        



        
    chains=np.zeros(((Nit-burn_in)*nwalkers,ndim))
    chains2=np.zeros((Nit-burn_in, nwalkers,ndim))
    lnpchain=np.zeros(((Nit-burn_in)*nwalkers))
    lnpchain2=np.zeros(((Nit-burn_in), nwalkers))
    


    chains[:,:]=samples[:,burn_in:,:].reshape((nwalkers*(Nit-burn_in), ndim),order='c')
    lnpchain[:]=lnprobs[:,burn_in:].reshape((nwalkers*(Nit-burn_in)),order='c')
    
    ibestparams=np.argmax(lnpchain)
    bestparams=chains[ibestparams,:]
    
    ######### save bestparams
    np.save(workdir+'bestparams.dat',bestparams)
    

    for j in list(range(nwalkers)):
        chains2[:,j,:]=samples[j,burn_in:,:].reshape((Nit-burn_in, ndim),order='c')
        lnpchain2[:,j]=lnprobs[j,burn_in:].reshape(((Nit-burn_in)),order='c')

        
    #fig=plt.figure(figsize=(10,8))
    #par_labels=names
    #for i in list( range(nwalkers)):
    #    for ip in list(range(ndim)):
    #        ax=fig.add_subplot(ndim+1,1,ip+1)
    #        ax.plot(chains2[:,i,ip],alpha=0.1)
    #        ax.set_ylabel(par_labels[ip])
    #        
    #        ax=fig.add_subplot(ndim+1,1,ndim+1)
    #        ax.plot(lnpchain2[:,i],alpha=0.1)
    #        ax.set_ylabel('ln(p)')


    fig=plt.figure(figsize=(10,8))
    par_labels=names
    ax_lnprob=fig.add_subplot(ndim+1,1,ndim+1)
    for ip in list(range(ndim)):
        ax_chain=fig.add_subplot(ndim+1,1,ip+1)
        for i in list( range(nwalkers)):
            ax_chain.plot(chains2[:,i,ip],alpha=0.1)
            ax_chain.set_ylabel(par_labels[ip])
            ax_lnprob.plot(lnpchain2[:,i],alpha=0.1)
            ax_lnprob.set_ylabel('ln(p)')



            
    #plt.show()
    plt.savefig(workdir+'chains.png', bbox_inches='tight')
    plt.close(fig)




    #samples = sampler.chain[:, burn_in:, :].reshape((-1, ndim))

    
    mcmc_results = list(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(chains, [16, 50, 84],
                                                axis=0))))

    mcmc_results_0 = np.zeros(nvar)

    print( "param     distrib     max ")
    for iparam in list(range(nvar)):
        print( names[iparam],mcmc_results[iparam],bestparams[iparam])
        mcmc_results_0[iparam]= mcmc_results[iparam][0]
        

    #print( "mcmc median values:")
    #model_median =  np.array(modelfunk(mcmc_results_0, m))


    
    import corner

    labels=names
    for ilabel,alabel in enumerate(labels):
        if (alabel eq 'dra_off'):
            labels[ilabel] = r"$\Delta \alpha$"
        if (alabel eq 'ddec_off'):
            labels[ilabel] = r"$\Delta \delta$"
        

    
    fig=corner.corner(chains,
                      labels=names,
                      quantiles=[0.16, 0.5,0.84],
                      bins=20, truths=bestparams,
                      levels=[0.68, 0.95, 0.997],
                      show_titles=True,
                      title_fmt=".3f",
                      title_kwards={"fontsize": 10}) #, smooth=1.0




    fig.savefig(workdir+OptimM.TriangleFile)

    print( "finished MCMC for region workdir",workdir)
    return [names,mcmc_results]



def exec_Grid(M,OptimM):
    
    names = list(map( (lambda x: x[0]),OptimM.domain))
    bnds = list(map( (lambda x: x[1]),OptimM.domain))
    nvar=len(list(names))

    if (nvar > 2):
        sys.exit("only 2D grid for now")

    sample_theta=list(range(nvar))
    for iparam in list(range(nvar)):
        sample_theta[iparam]=getattr(M,names[iparam])

    x = np.array( sample_theta  ) 

    M.Verbose=False
    M.DumpAllFitsFiles=False
    M.Grid=True
    M.PlotAzimuthalProfile=False
    M.PlotRadialProfile=False
    M.XCheckInv=False

    M.prep_files()

    print("Init Grid with params:")
    for iparam in list(range(nvar)):
        print("x",x[iparam],"bnds",bnds[iparam])


    nmesh=OptimM.nmesh_grid
    
    dims=nmesh*(np.ones(nvar,dtype=int))
    print("dims",dims)
    print("dims.tolist",dims.tolist())
    chi2map=np.zeros(dims.tolist())

    print("Chi2map shape",chi2map.shape)
    
    xs=bnds[0][0]+(bnds[0][1]-bnds[0][0])*np.arange(nmesh)/(nmesh-1)
    ys=bnds[1][0]+(bnds[1][1]-bnds[1][0])*np.arange(nmesh)/(nmesh-1)
    
    xxs, yys = np.meshgrid(xs, ys)
    
    for ix in list(range(nmesh)):
        for iy in list(range(nmesh)):
            theta=[xs[ix],ys[iy]]
            statusstring='ix '+str(ix)+' iy '+str(iy)+' '
            for iparam in range(nvar):
                setattr(M,names[iparam],theta[iparam])
                #statusstring=statusstring+names[iparam]+" "+str(theta[iparam])+" "
                statusstring=statusstring+names[iparam]+" %.3f " %(theta[iparam])
            chi2=M.polar_expansions()
            statusstring=statusstring+" -> "+str(chi2)
            chi2map[iy,ix]=chi2
            if (OptimM.PrintOptimStatus):
                print( statusstring)

    hdu=pf.PrimaryHDU()
    hdu.data=chi2map
    mapheader=hdu.header
    mapheader['CRPIX1']=1
    mapheader['CRPIX2']=1
    mapheader['CRVAL1']=bnds[0][0]
    mapheader['CRVAL2']=bnds[1][0]
    mapheader['NAXIS1']=nmesh
    mapheader['NAXIS2']=nmesh
    mapheader['CDELT1']=xs[1]-xs[0]
    mapheader['CDELT2']=ys[1]-ys[0]
    hdu.header=mapheader
    hdu.writeto(M.workdir+'chi2map.fits',overwrite=True)

    statusstring=''
    result_grid=np.zeros(nvar)

    indminimum = np.unravel_index(np.argmin(chi2map, axis=None), chi2map.shape)
    result_grid[0]=xxs[indminimum]
    result_grid[1]=yys[indminimum]
    print("result_grid",result_grid)
    #iminimum=np.argmin(chi2map)
    #result_grid[0]=xxs.flatten()[iminimum]
    #result_grid[1]=yys.flatten()[iminimum]
    #print("result_grid",result_grid)

    for iparam in range(nvar):
        if OptimM.SetOptim:
            setattr(M,names[iparam],result_grid[iparam])
        statusstring=statusstring+names[iparam]+" %.3f " %(result_grid[iparam])


    statusstring=statusstring+"\n"
    fout=open(M.workdir+"result_grid.txt","a+")
    fout.write(statusstring)
    fout.close()

    np.save(M.workdir+'result_grid.dat',result_grid)

    M.Verbose=False
    M.DumpAllFitsFiles=True
    M.Grid=False
    M.PlotAzimuthalProfile=True
    M.PlotRadialProfile=True
    M.XCheckInv=False

    fig=plt.figure(figsize=(10,8))

    ax = plt.subplot(1, 1, 1)

    #plt.gca().set_xlabel(names[0])
    #plt.gca().set_ylabel(names[1])
    ax.set_xlabel(names[0])
    ax.set_ylabel(names[1])

    x1=bnds[0][0]
    x2=bnds[0][1]
    y1=bnds[1][0]
    y2=bnds[1][1]
    cmap='RdBu_r'
    cmap='ocean'
    vmin=chi2map.min()
    vmax=chi2map.max()

    print(" display range should be:",vmin,vmax)

    rasterimage=plt.imshow(chi2map,origin='lower', cmap=cmap,vmin=vmin,vmax=vmax,extent=[x1,x2,y1,y2])

    #plt.xlim(x1,x2)
    #plt.ylim(y1,y2)
    #from mpl_toolkits.axes_grid1 import make_axes_locatable
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="3%", pad=0.2)
    #cax.yaxis.set_ticks_position('left')
    #cax.xaxis.set_ticks_position('top')
    #cax.xaxis.set_tick_params(labelsize=12, direction='in')
    #cax.yaxis.set_tick_params(labelsize=12, direction='in')
    #fmt='%.2f'
    #cb = plt.colorbar(rasterimage, cax=cax, format=fmt, ticks=clevs)
    cb = plt.colorbar(rasterimage)


    fileout = M.workdir+'fig_chi2map.pdf'

    print(fileout)

    plt.savefig(fileout)

    return

