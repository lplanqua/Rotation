from casatasks import specfit, exportfits
import glob, os, sys

args = sys.argv

starname = 'R_Dor' #or R_Dor or Betelgeuse or R_Leo

date = '2023-10-12'


if (len(args) > 1):
    if "-i" in args:
        starname = args[args.index("-i") + 1 ]

    if "-date" in args:
        date = args[args.index("-date") + 1 ]
    if "-line" in args:
        transition = args[args.index("-line") + 1 ]



if starname == 'R_Leo':
  date = '2025-07-15'
# transition = '29SiO_v=1_8-7'



# transition = 'CO_v=1_3-2'
filepath = 'DATA/'+ starname+ '/' + date + '/'



basename= filepath + starname + "_"+ transition +  ".clean"

image=basename+".image"
#basename="R_Dor_29SiO_v=1_8-7.clean"
#  basename="Betelgeuse.29SiOv0"
#  basename="Betelgeuse.28SiOv1"
#  basename="Betelgeuse.12CO"
pest=[10,50,25,"lorentzian"]  # initial fit estimates
pest = [0.2,36,12,"gaussian"]
#  regs='circle[[253pix,254pix],15pix]' 
if '-l' in args:
  # suffix = '.clean.large_scale.fit2'
  basename = basename+ '.large_scale'

if starname == 'R_Dor':

  if date == '2023-10-12':
    center = '[127pix,146pix]'
  elif date == '2023-10-06' or date == '2023-10-04' :
    center = '[121pix,141pix]'
  elif date == '2023-10-09':
    center = '[116pix,143pix]'

  regs='circle['+ center+ ',8pix]'
    #

elif starname == 'R_Leo':
  center = '[156pix,156pix]'
  regs='circle['+ center+ ',7pix]'


if 'large_scale' in basename:
  regs = 'annulus['+ center+ ',[7pix,12pix]]'
  regs='circle['+ center+ ',12pix]'

amprange=[-0,0.5]
centerrange=[0]#[25,65]   # in pixels!
fwhmrange = [2,25]


#name of the output file

#

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
# print(stop)
print('starting the fitting')

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
        goodfwhmrange=fwhmrange,
        sigma="",
        outsigma="")
print('AGB_specfit finished')

filenames = glob.glob('*'+ fitbase +'*')

for filename in filenames:
    nameout = filename+".fits"

    if not os.path.exists(nameout):
        print(nameout)
        exportfits(imagename=filename, fitsimage=nameout, dropstokes =True, velocity = True)
        os.system("rm -rf "+filename)



for log_file in glob.glob("casa-*.log"):
  os.remove(log_file)


"""
mystep = 3
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print('Step ', mystep, step_title[mystep])

#  basename="RDorB6_29SiOv0"
#  basename="RDorB6_SO2_1616"
#  basename="RDorB6_SO2_2626"
  basename="R_Dor_SiO_v=2_8-7"
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

  basename="R_Dor_SiO_v=2_8-7"
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

"""
