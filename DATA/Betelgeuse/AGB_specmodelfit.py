thesteps = []
step_title = {0: 'Spectral fitting',
              1: 'velocity max moment',
              3: '2 component fit'}

try:
    print('List of steps to be executed ...', mysteps)
    thesteps = mysteps
except:
    print('global variable mysteps not set.')
  
if (thesteps==[]):
    thesteps = range(0,len(step_title))
    print('Executing all steps: ', thesteps)

import re

import os



mystep = 0
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  basename="Betelgeuse.28SiOv2"
#  basename="Betelgeuse.29SiOv0"
#  basename="Betelgeuse.28SiOv1"
#  basename="Betelgeuse.12CO"
  pest=[10,50,25,"lorentzian"]  # initial fit estimates
#  regs='circle[[253pix,254pix],15pix]' 
#  regs='circle[[125pix,125pix],15pix]' 
  regs='annulus[[253pix,254pix],[8pix,16pix]]'
#  regs='annulus[[125pix,125pix],[7pix,15pix]]' 
  amprange=[0.003,0.02]
  centerrange=[25,65]   # in pixels!
  

  image=basename+".image"
  fitbase=basename+".fit2"
  modelfile=fitbase+".model.image"
  residualfile=fitbase+".model.residual"
  ampfile=fitbase+".amp"
  amperrfile=fitbase+".amp.err"
  centerfile=fitbase+".center.vel"
  centererrfile=fitbase+".center.err"
  fwhmfile=fitbase+".fwhm"
  fwhmerrfile=fitbase+".fwhm.err"

  os.system("rm -rf "+fitbase+"*")

  specfit(imagename=image,
          box="",
          region=regs,
          chans="",
          stokes="",
          axis=-1,
          mask="",
          ngauss=1,
          poly=0,
          estimates="",
          minpts=1,
          multifit=True,
          model=modelfile,
          residual=residualfile,
          amp=ampfile,
          amperr=amperrfile,
          center=centerfile,
          centererr=centererrfile,
          fwhm=fwhmfile,
          fwhmerr=fwhmerrfile,
          integral="",
          integralerr="",
          wantreturn=False,
          stretch=False,
          logresults=True,
          pampest=pest[0],
          pcenterest=pest[1],
          pfwhmest=pest[2],
          pfix={},
          pfunc=pest[3],
          gmfix="",
          logfile="",
          append=True,
          goodamprange=amprange,
          goodcenterrange=centerrange,
          goodfwhmrange=[0.0],
          sigma="",
          outsigma="")


mystep = 1
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

#  basename="RDorB6_29SiOv0"
#  basename="RDorB6_SO2_1616"
#  basename="RDorB6_SO2_2626"
  basename="RDorB6_SO3S54"
#  basename="RDorB6_SiOv1"
  image=basename+".clean.image"
  maskpos=basename+".mask1"
  maskneg=basename+".mask2"
#  plev=0.01
#  mlev=-0.01
  plev=0.015
  mlev=-0.015

  os.system("rm -rf "+image+".mom*")

  immoments(imagename=image,
            moments=[8],
            excludepix = [-1000,plev],
            outfile=image+'.mom8')

  immoments(imagename=image,
            moments=[9],
            excludepix = [-1000,plev],
            outfile=image+'.mom9')

  immoments(imagename=image,
            moments=[10],
            excludepix = [mlev,1000],
            outfile=image+'.mom10')

  immoments(imagename=image,
            moments=[11],
            excludepix = [mlev,1000],
            outfile=image+'.mom11')


  
mystep = 3
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

  basename="RDorB6_29SiOv1"
  pest1=[-42,36,12,"gaussian"]
  pest2=[13,44,8,"gaussian"]
  masks='circle[[121pix,138pix],8pix]' 

  image=basename+".clean.image"
  fitbase=basename+".2comp.fit"
  modelfile=fitbase+".model.image"
  residualfile=fitbase+".model.residual"
  ampfile=fitbase+".amp"
  amperrfile=fitbase+".amp.err"
  centerfile=fitbase+".center.vel"
  centererrfile=fitbase+".center.err"
  fwhmfile=fitbase+".fwhm"
  fwhmerrfile=fitbase+".fwhm.err"

  os.system("rm -rf "+fitbase+"*")

  specfit(imagename=image,
          box="",
          region=masks,
          chans="",
          stokes="",
          axis=-1,
          mask="",
          ngauss=2,
          poly=0,
          estimates="",
          minpts=1,
          multifit=True,
          model=modelfile,
          residual=residualfile,
          amp=ampfile,
          amperr=amperrfile,
          center=centerfile,
          centererr=centererrfile,
          fwhm=fwhmfile,
          fwhmerr=fwhmerrfile,
          integral="",
          integralerr="",
          wantreturn=False,
          stretch=False,
          logresults=True,
          pampest=[pest1[0],pest2[0]],
          pcenterest=[pest1[1],pest2[1]],
          pfwhmest=[pest1[2],pest2[2]],
          pfix="",
          pfunc=[pest1[3],pest2[3]],
          gmncomps=0,
          gmampcon="",
          gmcentercon="",
          gmfwhmcon="",
          gmampest=[0.0],
          gmcenterest=[0.0],
          gmfwhmest=[0.0],
          gmfix="",
          logfile="",
          append=True,
          goodamprange=[0.0],
          goodcenterrange=[0.0],
          goodfwhmrange=[0.0],
          sigma="",
          outsigma="")


