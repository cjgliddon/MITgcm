import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter
import MITgcmutils as mu
import cartopy.crs as ccrs
import string
import xarray as xr
import gsw as gsw
from scipy.interpolate import griddata
from scipy.interpolate import interp1d,interp2d,RectBivariateSpline
from scipy.interpolate import RegularGridInterpolator
import scipy.ndimage.filters as filters
from IPython.display import display
import pandas as pd
from os import path
import struct
import re
import glob
# haven't incorporated the complex ptide_ext, pcond_ext thing
class Experiment():
    def __init__(self, **kwargs):
        self.base_dir = kwargs.get('base_dir',
                                   "./data_")
        base_dir = self.base_dir
        self.name = kwargs.get("name",None)
        self.path = base_dir + self.name + "/"
        update=kwargs.get("update",False)
        runpath=kwargs.get("runpath",None)
        if runpath is None:
            runpath='/'
            if not path.exists(self.path+runpath+"XC.data"):
                runpath='coords/'
                if not path.exists(self.path+runpath+"XC.data"):
                    runpath='run0/'
                
        self.xc = mu.rdmds(self.path+runpath + "XC")
        self.yc = mu.rdmds(self.path+runpath + "YC")
        self.ny,self.nx = np.shape(self.yc)
        self.xg = mu.rdmds(self.path+runpath + "XG")
        self.yg = mu.rdmds(self.path+runpath + "YG")
        self.yg = np.append(self.yg,-self.yg[0]*np.ones((1,len(self.xc[0,:]))),axis=0)
        self.zc = mu.rdmds(self.path+runpath + "RC").squeeze()
        self.nz = len(self.zc)
        self.zf = mu.rdmds(self.path+runpath + "RF").squeeze()
        self.dzc = -np.diff(self.zf)
        self.dzf = -np.diff(self.zc)
        self.dxc = mu.rdmds(self.path+runpath + "DXC")
        self.dyc = mu.rdmds(self.path+runpath + "DYC")
        self.dxg = mu.rdmds(self.path+runpath + "DXG")
        self.dyg = mu.rdmds(self.path+runpath + "DYG")
        rac = mu.rdmds(self.path+runpath + "RAC")
        self.ras = mu.rdmds(self.path+runpath + "RAS").squeeze()
        self.raw = mu.rdmds(self.path+runpath + "RAW").squeeze()
        self.raz = mu.rdmds(self.path+runpath + "RAZ").squeeze()
        self.rhoprof0=mu.rdmds(self.path+runpath+'/RhoRef').squeeze()
        hfacc = mu.rdmds(self.path+runpath + "hFacC")
        hfacw = mu.rdmds(self.path+runpath + "hFacW")
        drf = mu.rdmds(self.path+runpath + "DRF").squeeze()
        self.drf=xr.DataArray(drf,dims=('z'),coords=dict(z=self.zc))
        self.spheric=True
        self.GM=False

        self.cp = 4000
        self.lf=334000
        self.beta=0.0
        self.topo=False
        self.Dfig = kwargs.get('Dfig', 500000)  # Used only for polar plots
        self.deepatm=False
        self.a=1
        self.omega=0
        self.Q0base=0
        self.difconv=0
        with open(self.path + "data") as namelist:
            for line in namelist:
                if '#' in line:
                    continue
                if 'dumpFreq=' in line:
                    self.dumpfreq = float(re.search('=(.*),',line).group(1))
                if 'gravity=' in line:
                    self.g = float(re.search('=(.*),',line).group(1))
                if 'rSphere=' in line:
                    self.a = float(re.search('=(.*),',line).group(1))
                if 'tAlpha=' in line:
                    self.alpha = float(re.search('=(.*),',line).group(1))
                if 'sBeta=' in line:
                    self.beta = float(re.search('=(.*),',line).group(1))
                if 'HeatCapacity_Cp=' in line:
                    try:
                        self.cp = float(re.search('=(.*),',line).group(1))
                    except:
                        self.cp=4000
                if 'rhoNil=' in line:
                    self.rhoref=float(re.search('=(.*),',line).group(1))
                if 'eosType=' in line:
                    self.eos=re.search('=(.*),',line).group(1)
                    self.eos=self.eos[1:-1]
                if 'viscAh=' in line:
                    self.ah=float(re.search('=(.*),',line).group(1))
                if 'viscAr=' in line or 'viscAz=' in line:
                    self.av=float(re.search('=(.*),',line).group(1))
                if 'viscArNr=' in line:
                    tmp=re.search('=(.*),',line).group(1)
                    av1=re.split(',',tmp)[0]
                    av2=re.split('\*',av1)[1]
                    self.av=float(av2)
                if 'diffKhT=' in line:
                    self.difh=float(re.search('=(.*),',line).group(1))
                if 'diffKrT=' in line or 'diffKzT=' in line:
                    self.difv=float(re.search('=(.*),',line).group(1))
                if '#' not in line and 'ivdc_kappa=' in line:
                    self.difconv=float(re.search('=(.*),',line).group(1))
                if 'bottomDragLinear' in line:
                    self.gammabot=float(re.search('=(.*),',line).group(1))
                if 'rotationPeriod=' in line:
                    self.omega=2*np.pi/float(re.search('=(.*),',line).group(1))
                if 'deepAtmosphere=' in line:
                    tmp=re.search('=.(.*).,',line).group(1)
                    if tmp=='False' or tmp=='FALSE':
                        self.deepatm=False
                    else:
                        self.deepatm=True
                if 'tRef=' in line:
                    try:
                        self.Tref=float(re.search('\*(.*),',line).group(1))
                    except:
                        tmp=re.search('\*(.*),',line).group(1)
                        self.Tref=float(re.search('(.*),',tmp).group(1))
                    self.Ts0=self.Tref
                if 'sRef=' in line:
                    self.Sref=float(re.search('\*(.*),',line).group(1))
                if 'deltaT=' in line:
                    self.dt=float(re.search('=(.*),',line).group(1))
                if 'nIter0=' in line:
                    tmp=re.search('=(.*)',line).group(1)
                    if ',' in tmp:
                        tmp=tmp[:-1]
                    self.itermax=int(tmp)
                if 'nTimeSteps=' in line:
                    self.eachiter=int(re.search('=(.*),',line).group(1))
                if 'dumpFreq=' in line:
                    self.dtoutput=round(float(re.search('=(.*),',line).group(1)))
                    self.iteroutput=self.dtoutput/self.dt
                    #self.iteroutput=self.iteroutput*10
                if 'geothermalFile' in line:
                    tmp=re.search('bottom_(.*)mW',line).group(1)
                    if 'W' in tmp:
                        tmp=re.search('(.*)mW_',tmp).group(1)
                    self.Q0base=float(tmp)/1e3
                if 'delX' in line:
                    self.delx=float(re.search('\*(.*),',line).group(1))
                    self.xg = np.append(self.xg,self.xg+self.delx*np.ones((len(self.yc[:,0]),1)),axis=1)
        self.runmax=int(self.itermax/self.eachiter)-1
        files=glob.glob(self.path+'run*[0-9]')
        self.iter0=0
        self.iter0s=[int(re.search('run(.*)',file).group(1)) for file in files]
        if len(self.iter0s)>0:
            self.iter0=min(self.iter0s)
        self.runlength=(self.itermax-self.iter0)*self.dt/86400/360
        self.difconv=self.difv

        if path.exists(self.path + "data.pkg"): 
            with open(self.path + "data.pkg") as pkgdata:
                for line in pkgdata:
                    if 'useGMRedi' in line and '#' not in line:
                        if 'TRUE' in line or 'True' in line:
                            self.GM=True
                        else:
                            self.GM=False

        self.Ntrc=0
        self.trcnames=[]
        if path.exists(self.path + "data.ptracers"): 
            with open(self.path + "data.ptracers") as ptracersdata:
                for line in ptracersdata:
                    if 'PTRACERS_numInUse=' in line and '#' not in line:
                        self.Ntrc=int(re.search('=(.*),',line).group(1))
                    if self.Ntrc!=0:
                        for ipt in range(1,self.Ntrc+1):
                            if 'PTRACERS_names('+str(ipt)+')=' in line:
                                self.trcnames.append(str(re.search("='(.*)',",line).group(1)))

        self.meridionalTs=0
        self.tide2d=False
        self.uniHtide=False
        self.pcond=-1.0
        self.realtopo=0
        self.Hice_P1=0.0
        self.Hice_P2=0.0
        self.Hice_P3=0.0
        self.Hice_k2=0.0
        self.mixbend=False
        self.addmixbend=False
        self.Hice_gaussA=0
        self.Hice_slope=0
        self.Htide0_portion=None
        self.icetopofile=None
        self.icemassfile=None
        self.rhoice=920
        
        if path.exists(self.path + "data.shelfice"): 
            self.haveice=True
            with open(self.path + "data.shelfice") as shelficedata:
                for line in shelficedata:
                    if '#' in line:
                        continue
                    if 'SHELFICEtopoFile=' in line:
                        self.icetopofile=re.search("='(.*)',",line).group(1)
                    if 'SHELFICEmassFile=' in line:
                        self.icemassfile=re.search("='(.*)',",line).group(1)
                    if 'rhoShelfIce=' in line:
                        self.rhoice=float(re.search('=(.*),',line).group(1))
                    if 'SHELFICEthetaSurface=' in line:
                        self.Ts0=float(re.search('=(.*),',line).group(1))+273.15
                        self.Ts=self.Ts0*np.ones_like(self.yc[:,0])
                    if 'SHELFICEheatTransCoeff' in line:
                        self.gammaT=float(re.search('=(.*),',line).group(1))
                    if 'SHELFICEDragLinear' in line:
                        self.gammatop=float(re.search('=(.*),',line).group(1))
                    if 'ptide' in line:
                        self.ptide=float(re.search('=(.*),',line).group(1))
                    if 'pcond' in line:
                        self.pcond=float(re.search('=(.*),',line).group(1))
                    if 'obliquity' in line:
                        self.obliquity=float(re.search('=(.*),',line).group(1))
                    if 'uniHtide' in line:
                        tmp=re.search('=.(.*).,',line).group(1)
                        if tmp == 'False' or tmp == 'FALSE':
                            self.uniHtide=False
                        else:
                            self.uniHtide=True
                    if 'tide2d' in line:
                        tmp=re.search('=.(.*).,',line).group(1)
                        if tmp == 'False' or tmp == 'FALSE':
                            self.tide2d=False
                        else:
                            self.tide2d=True
        else:
            self.haveice=False
            self.Hice0=0

        self.startevolve=0
        if path.exists(self.path + "startevolve"):
            with open(self.path + "startevolve") as startevolve:
                for line in startevolve:
                    self.startevolve=float(line)
            

        if path.exists(self.path + "gendata.m"):
            with open(self.path + "gendata.m") as gendata:
                for line in gendata:
                    if '%%%%%%' in line:
                        break
                    if line[0]=='%':
                        continue
                    if 'Htot=' in line:
                        self.Htot=float(re.search('=(.*);',line).group(1))
                    else:
                        self.Htot=np.max(self.zf) - np.min(self.zf)
                    if 'a0=' in line and 'kappa' not in line and 'fa0' not in line:
                        self.rsphere=float(re.search('=(.*);',line).group(1))
                    if 'kappa0=' in line:
                        self.kappa0=float(re.search('=(.*);',line).group(1))
                        self.fluxfac=self.kappa0/651
                    if 'Hice0=' in line or 'Hice=' in line:
                        self.Hice0=float(line[line.find('Hice0=')+len('Hice0='):line.find(';')])
                    if 'Htide0_portion=' in line:
                        tmp=re.search('=(.*)',line).group(1)
                        if ';' in tmp:
                            tmp=re.search('=(.*);',line).group(1)
                        self.Htide0_portion=float(tmp)
                    if 'Hice_P1=' in line and '==' not in line and line[0] != '%':
                        icm=line.find('%')
                        icl=line.find(';')
                        if not (icm>=0 and icm<icl):
                            tmp=line[line.find('Hice_P1=')+len('Hice_P1='):]
                            self.Hice_P1=float(tmp[:tmp.find(';')])
                    if 'Hice_P2=' in line and '==' not in line and line[0] != '%':
                        icm=line.find('%')
                        icl=line.find(';')
                        if not (icm>=0 and icm<icl):
                            tmp=line[line.find('Hice_P2=')+len('Hice_P2='):]
                            self.Hice_P2=float(tmp[:tmp.find(';')])
                    if 'Hice_P3=' in line and '==' not in line and line[0] != '%':
                        icm=line.find('%')
                        icl=line.find(';')
                        if not (icm>=0 and icm<icl):
                            tmp=line[line.find('Hice_P3=')+len('Hice_P3='):]
                            self.Hice_P3=float(tmp[:tmp.find(';')])
                    if 'Hice_gaussA=' in line and '==' not in line:
                        self.Hice_gaussA=float(re.search('=(.*);',line).group(1))
                    if 'Hice_k2=' in line and '==' not in line:
                        self.Hice_k2=float(re.search('=(.*);',line).group(1))
                    if 'Hice_slope=' in line and '==' not in line:
                        self.Hice_slope=float(re.search('=(.*);',line).group(1))
                    if 'qbotvary=' in line and 'qbotvary==' not in line:
                        tmp=float(re.search('=(.*);',line).group(1))
                        self.qbotvary=tmp
                    else:
                        self.qbotvary=0
                    if 'qbot_gausssigma=' in line:
                        self.qbot_gausssigma=float(re.search('=(.*);',line).group(1))
                    if 'spheric=' in line and '==' not in line and 'replace_string' not in line and 'fprint' not in line:
                        tmp=float(re.search('=(.*);',line).group(1))
                        if tmp==0:
                            self.spheric=False
                            self.yc = self.yc - np.mean(self.yc)
                            self.xc = self.xc - np.mean(self.xc)
                            self.yg = self.yg - np.mean(self.yg)
                            self.xg = self.xg - np.mean(self.xg)
                    if 'mixbend=' in line and '==' not in line and 'replace_string' not in line and 'fprint' not in line:
                        tmp=float(re.search('=(.*);',line).group(1))
                        if tmp==1:
                            self.mixbend=True
                    if 'addmixbend=' in line and '==' not in line and 'replace_string' not in line and 'fprint' not in line:
                        tmp=float(re.search('=(.*);',line).group(1))
                        if tmp==1:
                            self.addmixbend=True
                        else:
                            self.addmixbend=False
                    if 'Htidemode=' in line:
                        if (not hasattr(self, 'Htidemode')) or self.mixbend:
                            tmp=re.search('=\[(.*)\];*',line).group(1)
                            self.Htidemode = [float(item) for item in tmp.split(',')]
                    else:
                        self.Htidemode = np.zeros(5)
                    if 'Hmixbendmode=' in line:
                        tmp=re.search('=\[(.*)\];*',line).group(1)
                        self.Hmixbendmode = [float(item) for item in tmp.split(',')]
                    else:
                        self.Hmixbend = np.zeros(5)
                    if 'realtopo=' in line and 'realtopo==' not in line:
                        self.realtopo=float(re.search('=(.*);',line).group(1))
                    if 'realtopopath=' in line:
                            self.realtopopath=re.search("='(.*)';",line).group(1)
                    if 'meridionalTs' in line and 'meridionalTs==' not in line and 'replace_string' not in line and 'fprint' not in line and 'if' not in line:
                        self.meridionalTs=float(re.search('=(.*);',line).group(1))

        if self.haveice and (self.Hice_P1!=0 or self.Hice_P2!=0 or self.Hice_P3!=0 or self.realtopo!=0 or self.Hice_gaussA or self.Hice_k2 or self.Hice_slope):
            self.topo=True
        self.addmixbend = kwargs.get('addmixbend', self.addmixbend)
        self.hfacc = xr.DataArray( hfacc, dims=('z', 'lat', 'lon'),
            coords=dict(z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        self.hfacf=np.vstack((self.hfacc.values,self.hfacc.values[-1,np.newaxis,:,:]))
        self.rac= xr.DataArray( rac, dims=('lat', 'lon'),
            coords=dict(lat=self.yc[:, 0], lon=self.xc[0, :]))
        self.hfacw = xr.DataArray( hfacw, dims=('z', 'lat', 'lon'),
            coords=dict(z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        self.mask=self.hfacc*0+1
        self.mask=self.mask.where(self.hfacc>0.0)


        # just for scale estimation, except Q0 is taken as the global mean heat flux released that is normalized by the surface area
        if self.spheric:
            self.wgt=len(self.yc[:,0])*(np.sin(np.radians(self.yg[1:,0]))-np.sin(np.radians(self.yg[:-1,0])))/(np.sin(np.radians(self.yg[-1,0]))-np.sin(np.radians(self.yg[0,0])))
            self.ainter=self.a-(self.rhoice/self.rhoref)*self.Hice0
            self.Q0=self.Q0base*((self.a-self.Htot)/self.ainter)**(2*self.deepatm)
            self.widthfac=360/(self.xg[0,-1]-self.xg[0,0])
            self.deepfacc = 1 + self.zc / self.rsphere
            self.deepfacf = 1 + self.zf / self.rsphere
            self.deepfac2c = self.deepfacc**2
            self.deepfac2f = self.deepfacf**2
            self.totarea_zc=np.sum(self.rac.values[np.newaxis,:,:]*self.deepfac2c[:,np.newaxis,np.newaxis]*self.widthfac*self.hfacc,axis=1).squeeze()
            self.totarea_zf=np.sum(self.rac.values[np.newaxis,:,:]*self.deepfac2f[:,np.newaxis,np.newaxis]*self.widthfac*self.hfacf,axis=1).squeeze()
        else:
            self.wgt=self.yc[:,0]*0+1
            self.Q0=self.Q0base
            self.widthfac=1.
            self.deepfacc = 1 + self.zc*0
            self.deepfacf = 1 + self.zf*0
            self.deepfac2c=self.deepfacc
            self.deepfac2f=self.deepfacf
            self.totarea_zc=np.sum(self.rac.values[np.newaxis,:,:]*self.deepfac2c[:,np.newaxis,np.newaxis]*self.widthfac*self.hfacc,axis=1).squeeze()
            self.totarea_zf=np.sum(self.rac.values[np.newaxis,:,:]*self.deepfac2f[:,np.newaxis,np.newaxis]*self.widthfac*self.hfacf,axis=1).squeeze()
        dV=self.rac.values[np.newaxis,:,:]*self.deepfac2c[:,np.newaxis,np.newaxis]*self.widthfac*self.dzc[:,np.newaxis,np.newaxis]*self.hfacc
        self.dV=xr.DataArray( dV, dims=('z', 'lat', 'lon'),coords=dict(z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        slat1=np.sin(np.radians(self.yc[:,0]))
        P1=slat1
        P2=1.5*slat1**2-0.5
        P3=2.5*slat1**3-1.5*slat1
        P4=(35*slat1**4-30.0*slat1**2 + 3)/8.0
        P6=(231.0*slat1**6.0-315.0*slat1**4.0+105.0*slat1**2.0-5.0)/16.0
#            if self.realtopo!=0:
#                f=open(self.realtopopath,'rb')
#                fomat='>'+720*1440*'d'
#                Hice_=np.array(struct.unpack(fomat,f.read())).reshape((1440,720)).transpose()
#                Hice_=Hice_*1e3
#                if self.realtopo==1:
#                    y_=np.arange(-90+0.125,90,0.25)
#                    y_=np.hstack((-90,y_,90))
#                    Hice_=np.mean(Hice_,1)
#                    Hice_=np.hstack((Hice_[0],Hice_,Hice_[-1]))
#                    fHice=interp1d(y_,Hice_)
#                    Hice=fHice(self.yc[:,0])
#                    self.Hice=Hice[:,np.newaxis]
#                else:
#                    y_=np.arange(-90+0.125,90,0.25)
#                    y_=np.hstack((-90,y_,90))
#                    x_=np.arange(0,360.25,0.25)
#                    Hice_=np.vstack((Hice_[np.newaxis,0,:],Hice_,Hice_[np.newaxis,-1,:]))
#                    Hice_=np.hstack((Hice_,Hice_[:,0,np.newaxis]))
#                    fHice=RectBivariateSpline(y_,x_,Hice_)
#                    self.Hice=fHice(self.yc[:,0],self.xc[0,:])
#                self.Hice0=np.mean(self.Hice*self.wgt[:,np.newaxis])
#
#            else:
#                Hice=self.Hice0+self.Hice_P1*P1+self.Hice_P2*P2+self.Hice_P3*P3
#                self.Hice=Hice[:,np.newaxis]

#            self.Hunder=np.sum((1-self.hfacc)*self.drf,axis=0)
#            self.Hice=self.Hunder/self.rhoice*self.rhoref
#            self.Hice=self.Hice.values
        if path.exists(self.path+'gravity_myMoon.bin'):
            fgrav=open(self.path+'gravity_myMoon.bin', mode='rb')
            self.gprofile=np.fromfile(fgrav,'>f4').reshape((len(self.zc)))

        if update and self.haveice:
            self.Hice=self.real_Hice(runpath=runpath).mean('lon').squeeze().values[:,np.newaxis]
            self.Mice=self.Hice*self.rhoice
        else:
            if self.haveice:
                Micefile=self.icemassfile
                if path.exists(self.path+Micefile):
                    ficemass=open(self.path+Micefile, mode='rb')
                else:
                    ficemass=open(self.path+'../'+Micefile, mode='rb')

                self.Mice=np.fromfile(ficemass,'>f4')
                self.Mice=self.Mice.reshape((len(self.yc[:, 0]),len(self.xc[0, :])))
                self.Hice=self.Mice/self.rhoice

        ktop=np.sum(self.hfacc==0,axis=0)
        self.colmask=np.where(ktop==self.nz,0,1)
        self.ktop=np.minimum(ktop,self.nz-1)
        if self.haveice:
            if self.meridionalTs==1:
                ycr=np.radians(self.yc[:,0])
                self.Ts[abs(self.yc[:,0])<=90-self.obliquity]=self.Ts0*np.cos(ycr[abs(self.yc[:,0])<=90-self.obliquity])**0.25
                self.Ts[abs(self.yc[:,0])>90-self.obliquity]=self.Ts0*(((np.pi/2-abs(ycr[abs(self.yc[:,0])>90-self.obliquity]))**2+(self.obliquity*np.pi/180)**2)/2)**0.125
            if self.meridionalTs==3:
                cobl=np.cos(np.radians(self.obliquity))
                p2obl=1.5*cobl**2-0.5
                p4obl=(35*cobl**4-30.0*cobl**2 + 3)/8.0
                p6obl=(231.0*cobl**6.0-315.0*cobl**4.0+105.0*cobl**2.0-5.0)/16.0
                self.Ts=self.Ts0*(1.0-(5/8)*p2obl*P2-(9/64)*p4obl*P4-(65/1024)*p6obl*P6)**0.25
            if hasattr(self, 'kappa0'):
                self.Hcond=self.kappa0*np.log((self.Tref+273.15)/(self.Ts[:,np.newaxis]))/self.Hice0*(self.Hice/self.Hice0)**self.pcond
                self.Hcond=self.Hcond.mean(1)
                self.Hcond0=np.mean(self.Hcond*self.wgt)
                self.Htide0=self.Hcond0*self.Htide0_portion
                self.showname='core:{}W, shell:{}W, visr:{}, vish:{}'.format(self.Q0,self.Htide0,self.av,self.ah)
            else:
                self.Htide0=0.0
                self.Hcond0=0.0

        self.D=self.Htot-self.Hice0
        if self.spheric:
            self.ro=self.rsphere-self.Hice0
            self.ri=self.rsphere-self.Htot
            self.geomratio = self.ri / self.ro
            self.thetatc = np.degrees(np.arccos(self.geomratio))

        clat1=np.cos(np.radians(self.yc[:,0]))
        slat1=np.sin(np.radians(self.yc[:,0]))
        Y00=np.ones_like(clat1)/np.sqrt(4*np.pi)
        Y20=(1.5*slat1**2-0.5)/np.sqrt(4*np.pi/5)
        Y40=(35/8*slat1**4-30/8*slat1**2+3/8)/np.sqrt(4*np.pi/9)
        if self.spheric:
            self.thetaeq3=np.degrees(np.arcsin(slat1[-1]/3))
            self.thetaeq2=np.degrees(np.arcsin(slat1[-1]/2))
            if self.tide2d:
                c2lon=np.cos(np.radians(self.xc[0,:])*2)
                c4lon=np.cos(np.radians(self.xc[0,:])*4)
                Y22=(3*clat1[:,np.newaxis]**2)*(2*c2lon[np.newaxis,:])/np.sqrt(96*np.pi/5)
                Y42=(7.5*(7*slat1[:,np.newaxis]**2-1)*clat1[:,np.newaxis]**2)*(2*c2lon[np.newaxis,:])/np.sqrt(1440*np.pi/9)
                Y44=(105*clat1[:,np.newaxis]**4)*(2*c4lon[np.newaxis,:])/np.sqrt(40320*4*np.pi/9)
            else:
                Y22=0
                Y42=0
                Y44=0

        self.Q=self.Q0*np.ones_like(self.yc[:,0])
        if self.qbotvary==1:
            #qprofile=1-25/2/(60-25/2)*(2*clat1**2-1)
            qprofile=1.08449 + 0.252257*np.cos(np.radians(2*(90-self.yc[:,0]))) + 0.00599489*np.cos(np.radians(4*(90-self.yc[:,0])));
            qprofile=qprofile/np.mean(qprofile*self.wgt)
            self.Q=self.Q*qprofile
        elif self.qbotvary==2:
            qprofile=np.exp(-(self.yc[:,0]-np.mean(self.yc[:,0]))**2/2/self.qbot_gausssigma**2)
            qprofile=qprofile/np.mean(qprofile*self.wgt)
            self.Q=self.Q*qprofile
        elif self.qbotvary==3:
            qprofile=np.exp(-(self.yc[:,0]-np.mean(self.yc[:,0]))**2/2/self.qbot_gausssigma**2)
            qprofile=qprofile*np.exp(-(self.xc[0,:]-np.mean(self.xc[0,:]))**2/2/self.qbot_gausssigma**2)
            qprofile=qprofile/np.mean(qprofile*self.wgt)
            self.Q=self.Q*qprofile
        elif self.qbotvary==4:
            qprofile=1+cosd(2*(90-self.yc[:,0]))
            qprofile=qprofile/np.mean(qprofile*self.wgt)
            self.Q=self.Q*qprofile
        elif self.qbotvary==5:
            qprofile=np.exp(-(self.yc[:,0]-90)**2/(2*self.qbot_gausssigma**2))+np.exp(-(self.yc[:,0]+90)**2/(2*self.qbot_gausssigma**2))
            qprofile=qprofile/np.mean(qprofile*self.wgt)
            self.Q=self.Q*qprofile
        elif self.qbotvary==6:
            qprofile=np.exp(-(self.yc[:,0]+90)**2/(2*self.qbot_gausssigma**2))
            qprofile=qprofile/np.mean(qprofile*self.wgt)
            self.Q=self.Q*qprofile
        
        if self.haveice:
            if self.uniHtide:
                Htideprof=np.sqrt(4*np.pi)*Y00[:,np.newaxis]
            else:
                Htideprof=np.sqrt(4*np.pi)*(Y00[:,np.newaxis]+self.Htidemode[0]*Y20[:,np.newaxis]+self.Htidemode[1]*Y40[:,np.newaxis]+self.Htidemode[2]*Y22+self.Htidemode[3]*Y42+self.Htidemode[4]*Y44)

            self.Htideprof0=Htideprof
            Htideprof2=Htideprof*(self.Hice/self.Hice0)**(-2.0)
            Htideprof1=Htideprof*(self.Hice/self.Hice0)**(-1.0)
            Htideprof15=Htideprof*(self.Hice/self.Hice0)**(-1.5)
            Htideprof=Htideprof*(self.Hice/self.Hice0)**self.ptide
            if self.addmixbend:
                Hmixbendprof=np.sqrt(4*np.pi)*(self.Hmixbendmode[0]*Y00[:,np.newaxis]+self.Hmixbendmode[1]*Y20[:,np.newaxis]+self.Hmixbendmode[2]*Y40[:,np.newaxis]+self.Hmixbendmode[3]*Y22+self.Hmixbendmode[4]*Y42+self.Hmixbendmode[5]*Y44)
                self.Hmixbendprof=Hmixbendprof
                Htideprof1=Htideprof1+Hmixbendprof
                Htideprof2=Htideprof2+Hmixbendprof
                Htideprof15=Htideprof15+Hmixbendprof
                Htideprof=Htideprof+Hmixbendprof

            Htideprof=Htideprof/np.mean(Htideprof*self.wgt[:,np.newaxis])
            Htideprof1=Htideprof1/np.mean(Htideprof1*self.wgt[:,np.newaxis])
            Htideprof2=Htideprof2/np.mean(Htideprof2*self.wgt[:,np.newaxis])
            Htideprof15=Htideprof15/np.mean(Htideprof15*self.wgt[:,np.newaxis])
            self.Htide=self.Htide0*np.mean(Htideprof,axis=1)
            self.Htide1=self.Htide0*np.mean(Htideprof1,axis=1)
            self.Htide2=self.Htide0*np.mean(Htideprof2,axis=1)
            self.Htide15=self.Htide0*np.mean(Htideprof15,axis=1)

        # nondimensional numbers
        self.B = self.g * self.alpha * (self.Q0) / self.rhoref / self.cp
        self.unonrot = (self.B * (self.D))**(1 / 3)
        self.gprimenonrot = (self.B**2 / (self.D))**(1 / 3) / self.alpha / self.g
        if self.omega!=0:
            self.rostar = (self.B * (2 * self.omega)**-3)**(1 / 2) / (self.D)
            self.E = self.av / (2 * self.omega * (self.D)**2)
            self.taylor = 4 * self.omega**2 * (self.D)**4 / self.av**2
            self.lrot = (self.B * (2 * self.omega)**-3)**(1 / 2)
            self.urot = (self.B * (2 * self.omega)**-1)**(1 / 2)
            self.ucone = 2 * self.omega * (self.D) * self.rostar**(1 / 2)
            self.lcone = (self.D) * np.sqrt(self.rostar)
            self.gprime = np.sqrt(self.B * 2 * self.omega) / self.alpha / self.g
            self.rastarq = self.alpha * self.g * self.Q0 / self.rhoref / self.cp / (
                self.omega**3) / (self.D**2)

        self.umax = kwargs.get('umax', None)
        self.wmax = kwargs.get('wmax', None)
        self.Tmax = kwargs.get('Tmax', None)
        self.Tmin = kwargs.get('Tmin', None)
        self.ummax = kwargs.get('ummax', None)
        self.Tmmax = kwargs.get('Tmmax', None)


    def real_Hice(self,iteration=np.inf,runpath=None,runnum=np.inf):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        self.Mice=self.get_var(runpath+'shiDiag',iteration=iteration,dimnum=2,rec=9)
        #self.Hice=self.get_var(runpath+'shiDiag',iteration=iteration,dimnum=2,rec=10)
        realHice=self.Mice/self.rhoice
        print('check whether Mice/rhoice=Hice')
        return realHice

    def real_Hcond(self,iteration=np.inf,runpath=None,runnum=np.inf,returndim=1):
        realHice=self.real_Hice(iteration=iteration,runpath=runpath,runnum=runnum)
        realHcond=self.kappa0*np.log((self.Tref+273.15)/(self.Ts[:,np.newaxis]))/self.Hice0*(realHice/self.Hice0)**self.pcond
        realHcond_zm=realHcond.mean('lon')
        realHcond0=(realHcond_zm.squeeze()*self.wgt).mean('lat')
        if returndim==0:
            return realHcond0
        elif returndim==1:
            return realHcond_zm
        else:
            return realHcond
        
    def real_Htide(self,iteration=np.inf,runpath=None,runnum=np.inf,returndim=1):
        realHice=self.real_Hice(iteration=iteration,runpath=runpath,runnum=runnum)
        profile=self.Htideprof0*(realHice/self.Hice0)**self.ptide
        if self.addmixbend:
            profile=profile+self.Hmixbendprof
        profile=profile/np.mean(profile*self.wgt[np.newaxis,:,np.newaxis])
        realHtide=self.Htide0*profile
        realHtide_zm=realHtide.mean('lon')
        realHtide0=(realHtide_zm.squeeze()*self.wgt).mean('lat')
        if self.Htide0_portion is not None and self.Htide0 > 0:
            tmp=realHtide0
            realHtide0=self.real_Hcond(iteration=iteration,runpath=runpath,runnum=runnum,returndim=0)*self.Htide0_portion
            realHtide=realHtide*realHtide0/tmp
            realHtide_zm=realHtide_zm*realHtide0/tmp
        if returndim==0:
            return realHtide0
        elif returndim==1:
            return realHtide_zm
        else:
            return realHtide

    def eta(self,iteration=np.inf,runpath=None,runnum=np.inf):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        ETANdata=self.get_var(runpath+'surfDiag',rec=0,iteration=iteration,sample=sample,dimnum=2)
        ETANdata.attrs['units']='m'
        ETANdata.attrs['showname']='Eta (m)'
        ETANdata.attrs['longshowname']='free surface elevation (m)'
        return ETANdata

    def energetics(self,iteration=np.inf,runpath=None,runnum=np.inf,onlyz=True,rminbalance=True,decim=1,title=''):
        # onlyz: used if eos isn't linear. it determines whether alpha and beta is computed using horizontally averaged T,S or only longitudinally averaged T,S.
        BQ=self.energetics_Q(iteration=iteration,runpath=runpath,runnum=runnum,onlyz=onlyz,rminbalance=rminbalance).values.round(decimals=decim)
        BS=self.energetics_S(iteration=iteration,runpath=runpath,runnum=runnum,onlyz=onlyz,rminbalance=rminbalance).values.round(decimals=decim)
        #BdT=self.energetics_diffT(iteration=iteration,runpath=runpath,runnum=runnum,onlyz=onlyz).values.round(decimals=decim)
        #BdS=self.energetics_diffS(iteration=iteration,runpath=runpath,runnum=runnum,onlyz=onlyz).values.round(decimals=decim)
        BdT=self.energetics_diffT_flux(iteration=iteration,runpath=runpath,runnum=runnum,onlyz=onlyz).values.round(decimals=decim)
        BdS=self.energetics_diffS_flux(iteration=iteration,runpath=runpath,runnum=runnum,onlyz=onlyz).values.round(decimals=decim)
        Bfric=self.energetics_fric(iteration=iteration,runpath=runpath,runnum=runnum).values.round(decimals=decim)
        Bvisc=self.energetics_visc_tend(iteration=iteration,runpath=runpath,runnum=runnum).round(decimals=decim)
        #Bconv=self.energetics_conv(iteration=iteration,runpath=runpath,runnum=runnum).values.round(decimals=decim)
        #BconvTS=self.energetics_convTS(iteration=iteration,runpath=runpath,runnum=runnum).values.round(decimals=decim)
        #Bdiff=self.energetics_diff(iteration=iteration,runpath=runpath,runnum=runnum).values.round(decimals=decim)
        Badv=self.energetics_adv(iteration=iteration,runpath=runpath,runnum=runnum,onlyz=onlyz).values.round(decimals=decim)
        Bgen=BdT+BdS+BQ+BS
##        Btot=BQ+BS+Bfric+Bvisc+BconvTS+BdT+BdS
##        dict={'Energy terms':['heat flux','salinity flux','T diffusion','S diffusion','friction','dissip (zeta)','convection','sum'],
##              'Contribution (W)':[BQ,BS,BdT,BdS,Bfric,Bvisc,BconvTS,Btot]}
        Btot=BQ+BS+Bvisc+BdT+BdS
        dict={'Energy terms':['heat flux','salinity flux','diff/conv T', 'diff/conv S','dissip/fric','sum','adv','gen','boundary fric'],
              'Contribution (W)':[BQ,BS,BdT,BdS,Bvisc,Btot,Badv,Bgen,Bfric]}
        df = pd.DataFrame(dict)
        df.columns.name = title
        display(df)
        return

    def energetics_fric(self,iteration=np.inf,runpath=None,runnum=np.inf,U=None, V=None):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        if U is None:
            U=self.U(iteration=iteration,runpath=runpath,runnum=runnum)
        if V is None:
            V=self.V(iteration=iteration,runpath=runpath,runnum=runnum)
        Utop=U.isel(z=self.ktop,lat=np.arange(0,self.ny),lon=np.arange(0,self.nx))
        Ubot=U.isel(z=self.nz-1)
        Vtop=V.isel(z=self.ktop,lat=np.arange(0,self.ny),lon=np.arange(0,self.nx))
        Vbot=V.isel(z=self.nz-1)
        frictop=self.gammatop*(Utop**2+Vtop**2).mean('Time')
        fricbot=self.gammabot*(Ubot**2+Vbot**2).mean('Time')
        areatop=(self.rac.squeeze()*self.widthfac*((self.a-(self.rhoice/self.rhoref)*self.Hice)**2/self.a**2).squeeze())
        areabot=(self.rac.squeeze()*self.widthfac*((self.a-self.Htot)**2/self.a**2))
        B_fric=-np.sum(self.rhoref*(frictop*areatop+fricbot*areabot))
        return B_fric

    def energetics_Q(self,iteration=np.inf,runpath=None,runnum=np.inf,onlyz=True, rminbalance=True):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        # ref height: bottom, flux upward positive
        Qsurf=self.heatflux_surf(iteration=iteration,runpath=runpath,runnum=runnum).mean('Time')
        area=(self.rac.squeeze()*self.widthfac*((self.a-(self.rhoice/self.rhoref)*self.Hice)**2/self.a**2).squeeze())
        areabot=(self.rac.squeeze()*self.widthfac*((self.a-self.Htot)**2/self.a**2))
        if rminbalance:
            Qsurf=Qsurf-(Qsurf*area).mean()/area.mean()+(self.Q0base*areabot).mean()/area.mean()
        if self.eos.lower() == 'linear':
            alpha=self.alpha
        else:
            alpha=self.nonlinear_alpha(iteration=iteration,runpath=runpath,runnum=runnum, onlyz=onlyz)
            if onlyz:
                alpha=alpha[:,np.newaxis,np.newaxis]
            else:
                alpha=alpha[:,:,np.newaxis]
        alphagdz=(alpha*self.gprofile[:,np.newaxis,np.newaxis]*self.dzc[:,np.newaxis,np.newaxis]*self.hfacc)
        alphagint=np.sum(alphagdz,axis=0)
        B_Qtop=np.sum(Qsurf/(self.rhoref*self.cp)*alphagint*area*self.rhoref)

        return B_Qtop
        
    def energetics_S(self,iteration=np.inf,runpath=None,runnum=np.inf,onlyz=True, rminbalance=True):
        # ref height: bottom, flux upward positive
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        Ssurf=(-self.Sref*self.q(iteration=iteration,runpath=runpath,runnum=runnum,inmassunit=True).mean('Time')/self.rhoref)
        area=(self.rac.squeeze()*self.widthfac*((self.a-(self.rhoice/self.rhoref)*self.Hice)**2/self.a**2).squeeze())
        areabot=(self.rac.squeeze()*self.widthfac*((self.a-self.Htot)**2/self.a**2))
        if rminbalance:
            Ssurf=Ssurf-(Ssurf*area).mean()/area.mean()
        if self.eos.lower() == 'linear':
            beta=self.beta
        else:
            beta=self.nonlinear_beta(iteration=iteration,runpath=runpath,runnum=runnum,onlyz=onlyz)
            if onlyz:
                beta=beta[:,np.newaxis,np.newaxis]
            else:
                beta=beta[:,:,np.newaxis]
        betadz=(beta*self.gprofile[:,np.newaxis,np.newaxis]*self.dzc[:,np.newaxis,np.newaxis]*self.hfacc)
        betaint=np.sum(betadz,axis=0)
        B_Stop=np.sum(-Ssurf*betaint*area*self.rhoref)

        return B_Stop

    def energetics_convTS(self,iteration=np.inf,runpath=None,runnum=np.inf,onlyz=True, T=None, S=None):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
##        if T is None:
##            T=self.T(iteration=iteration,runpath=runpath,runnum=runnum).mean('Time').mean('lon')
##        if S is None:
##            S=self.S(iteration=iteration,runpath=runpath,runnum=runnum).mean('Time').mean('lon')
##        T_z=T.differentiate('z')
##        S_z=S.differentiate('z')
##        if self.eos.lower() == 'linear':
##            alpha=self.alpha
##            beta=self.beta
##        else:
##            alpha=self.nonlinear_alpha(iteration=iteration,runpath=runpath,runnum=runnum,T=T,S=S,onlyz=onlyz)
##            beta=self.nonlinear_beta(iteration=iteration,runpath=runpath,runnum=runnum,T=T,S=S,onlyz=onlyz)
##        Rho_z=-self.rhoref*alpha*T_z+self.rhoref*beta*S_z
        Rho_z=self.dRhodR(iteration=iteration,runpath=runpath,runnum=runnum).mean('Time').mean('lon')
        B_convTS=np.sum(-(self.difconv-self.difv)*Rho_z.where(Rho_z>1e-7)*self.dV.sum('lon')*self.gprofile[:,np.newaxis])
        return B_convTS

    def energetics_conv(self,iteration=np.inf,runpath=None,runnum=np.inf,onlyz=True, dRhodR=None):
        # meridional T S variation induced alpha change are not accounted, need to work on onlyz=False case 
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        if dRhodR is None:
            Rho_z=self.dRhodR(iteration=iteration,runpath=runpath,runnum=runnum).mean('Time').mean('lon')
        convdif=-(self.difconv-self.difv)*Rho_z.where(Rho_z>0)*self.gprofile[:,np.newaxis]*self.dV.sum('lon')
        #convdif.plot()
        B_conv=np.sum(convdif)

        return B_conv

    def energetics_adv(self,iteration=np.inf,runpath=None,runnum=np.inf,onlyz=True):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''

        ADVr_TH = self.get_var( runpath + 'flxDiag', iteration=iteration,rec=1).mean('Time')
        ADVr_SLT = self.get_var(runpath + 'flxDiag', iteration=iteration,rec=6).mean('Time') 
        if self.eos.lower() == 'linear':
            alphaf=self.alpha
        else:
            alpha=self.nonlinear_alpha(iteration=iteration,runpath=runpath,runnum=runnum,T=T,S=S,onlyz=onlyz)
            alphaf=(alpha[1:]+alpha[:-1])/2.0
            if onlyz:
                alphaf=alphaf[:,np.newaxis,np.newaxis]
            else:
                alphaf=alphaf[:,:,np.newaxis]

        if self.eos.lower() == 'linear':
            betaf=self.beta
        else:
            beta=self.nonlinear_beta(iteration=iteration,runpath=runpath,runnum=runnum,T=T,S=S,onlyz=onlyz)
            betaf=(beta[1:]+beta[:-1])/2.0
            if onlyz:
                betaf=betaf[:,np.newaxis,np.newaxis]
            else:
                betaf=betaf[:,:,np.newaxis]

        gf=(self.gprofile[1:]+self.gprofile[:-1])/2.0
        B_advT=np.sum((alphaf*gf[:,np.newaxis,np.newaxis])*ADVr_TH[1:]*self.dzf[:,np.newaxis,np.newaxis]*self.widthfac*self.rhoref)
        B_advS=np.sum((-betaf*gf[:,np.newaxis,np.newaxis])*ADVr_SLT[1:]*self.dzf[:,np.newaxis,np.newaxis]*self.widthfac*self.rhoref)

        B_adv=B_advS+B_advT
        return B_adv

    def energetics_diffT_flux(self,iteration=np.inf,runpath=None,runnum=np.inf,onlyz=True):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        DFrE_TH = self.get_var( runpath + 'flxDiag', iteration=iteration,rec=3).mean('Time')
        DFrI_TH = self.get_var(runpath + 'flxDiag', iteration=iteration,rec=4).mean('Time') 
        DFr_TH=DFrE_TH+DFrI_TH
        mask=self.hfacc*0+1
        #mask=self.mask.where(self.hfacf>0.01)
        DFr_TH=DFr_TH*mask

        # written on L grid (grid center - 1/2), higher in z 
        if self.eos.lower() == 'linear':
            alphaf=self.alpha
        else:
            alpha=self.nonlinear_alpha(iteration=iteration,runpath=runpath,runnum=runnum,T=T,S=S,onlyz=onlyz)
            alphaf=(alpha[1:]+alpha[:-1])/2.0
            if onlyz:
                alphaf=alphaf[:,np.newaxis,np.newaxis]
            else:
                alphaf=alphaf[:,:,np.newaxis]
        gf=(self.gprofile[1:]+self.gprofile[:-1])/2.0
        B_diffT=np.sum((-alphaf*gf[:,np.newaxis,np.newaxis])*DFr_TH[1:]*self.dzf[:,np.newaxis,np.newaxis]*self.widthfac*self.rhoref)
        return B_diffT

    def energetics_diffS_flux(self,iteration=np.inf,runpath=None,runnum=np.inf,onlyz=True):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        DFrE_SLT = self.get_var(runpath + 'flxDiag', iteration=iteration,rec=8).mean('Time')
        DFrI_SLT = self.get_var(runpath + 'flxDiag', iteration=iteration,rec=9).mean('Time') 
        DFr_SLT=DFrE_SLT+DFrI_SLT
        mask=self.hfacc*0+1
        #mask=self.mask.where(self.hfacc>0.01)
        DFr_SLT=DFr_SLT*mask

        # written on L grid (grid center - 1/2), higher in z 
        if self.eos.lower() == 'linear':
            betaf=self.beta
        else:
            beta=self.nonlinear_beta(iteration=iteration,runpath=runpath,runnum=runnum,T=T,S=S,onlyz=onlyz)
            betaf=(beta[1:]+beta[:-1])/2.0
            if onlyz:
                betaf=betaf[:,np.newaxis,np.newaxis]
            else:
                betaf=betaf[:,:,np.newaxis]
        gf=(self.gprofile[1:]+self.gprofile[:-1])/2.0
        B_diffS=np.sum((betaf*gf[:,np.newaxis,np.newaxis])*DFr_SLT[1:]*self.dzf[:,np.newaxis,np.newaxis]*self.widthfac*self.rhoref)
        return B_diffS

    def energetics_diff(self,iteration=np.inf,runpath=None,runnum=np.inf,onlyz=True, dRhodR=None):
        # meridional T S variation induced alpha change are not accounted, need to work on onlyz=False case 
        if dRhodR is None:
            Rho_z=self.dRhodR(iteration=iteration,runpath=runpath,runnum=runnum).mean('Time').mean('lon')
        B_diff=np.sum(-self.difv*Rho_z*self.gprofile[:,np.newaxis]*self.dV.sum('lon'))

        return B_diff

    def energetics_diffT(self,iteration=np.inf,runpath=None,runnum=np.inf,onlyz=True, T=None,S=None):
        # meridional T S variation induced alpha change are not accounted, need to work on onlyz=False case 
        if T is None:
            T=self.T(iteration=iteration,runpath=runpath,runnum=runnum).mean('Time')
        if S is None:
            S=self.S(iteration=iteration,runpath=runpath,runnum=runnum).mean('Time')
        T_z=T.differentiate('z')
        if self.eos.lower() == 'linear':
            alpha=self.alpha
        else:
            alpha=self.nonlinear_alpha(iteration=iteration,runpath=runpath,runnum=runnum,T=T,S=S,onlyz=onlyz)
            if onlyz:
                alpha=alpha[:,np.newaxis,np.newaxis]
            else:
                alpha=alpha[:,:,np.newaxis]
        B_diffT=np.sum(self.difv*alpha*T_z*self.dV*self.gprofile[:,np.newaxis,np.newaxis]*self.rhoref)

        return B_diffT

    def energetics_diffS(self,iteration=np.inf,runpath=None,runnum=np.inf,onlyz=True, T=None,S=None):
        # meridional T S variation induced beta change are not accounted, need to work on onlyz=False case 
        if T is None:
            T=self.T(iteration=iteration,runpath=runpath,runnum=runnum).mean('Time')
        if S is None:
            S=self.S(iteration=iteration,runpath=runpath,runnum=runnum).mean('Time')
        S_z=S.differentiate('z')
        if self.eos.lower() == 'linear':
            beta=self.beta
        else:
            beta=self.nonlinear_beta(iteration=iteration,runpath=runpath,runnum=runnum,T=T,S=S,onlyz=onlyz)
            if onlyz:
                beta=beta[:,np.newaxis,np.newaxis]
            else:
                beta=beta[:,:,np.newaxis]
        B_diffS=np.sum(-self.difv*beta*S_z*self.dV*self.gprofile[:,np.newaxis,np.newaxis]*self.rhoref)

        return B_diffS

    def energetics_visc_tend(self,iteration=np.inf,runpath=None,runnum=np.inf):
        if runpath is None:
            _,runpath=self.get_runpath(iteration,runnum=runnum)
            print(runpath)
        # read velocity
        _,its,meta=mu.rdmds(self.path + runpath + 'dynDiag', np.inf,returnmeta=True,rec=0)
        iUVEL=meta['fldlist'].index('UVEL')
        U=self.get_var(runpath+'dynDiag',iteration=iteration,rec=iUVEL).values 
        iVVEL=meta['fldlist'].index('VVEL')
        V=self.get_var(runpath+'dynDiag',iteration=iteration,rec=iVVEL).values 
        iWVEL=meta['fldlist'].index('WVEL')
        W=self.get_var(runpath+'dynDiag',iteration=iteration,rec=iWVEL).values 
        # read tendency
        _,its,meta=mu.rdmds(self.path + runpath + 'momDiag', np.inf,returnmeta=True,rec=0)
        iUm_Diss=meta['fldlist'].index('Um_Diss')
        UmD=self.get_var(runpath+'momDiag',iteration=iteration,rec=iUm_Diss).values
        iUm_ImplD=meta['fldlist'].index('Um_ImplD')
        UmI=self.get_var(runpath+'momDiag',iteration=iteration,rec=iUm_ImplD).values
        UmD=UmD+UmI
        iVm_Diss=meta['fldlist'].index('Vm_Diss')
        VmD=self.get_var(runpath+'momDiag',iteration=iteration,rec=iVm_Diss).values
        iVm_ImplD=meta['fldlist'].index('Vm_ImplD')
        VmI=self.get_var(runpath+'momDiag',iteration=iteration,rec=iVm_ImplD).values
        VmD=VmD+VmI
        iWm_Diss=meta['fldlist'].index('Wm_Diss')
        WmD=self.get_var(runpath+'momDiag',iteration=iteration,rec=iWm_Diss).values
#        if np.size(U,0)>1:
#            U=U[1:np.size(UmD,0)+1]
#            V=V[1:np.size(VmD,0)+1]
#            W=W[1:np.size(WmD,0)+1]

        # shift W grid
        WWmD=W*WmD
        WWmD[:,:-1,:,:]=0.5*(WWmD[:,1:,:,:]+WWmD[:,:-1,:,:])
        WWmD[:,-1,:,:]=0.5*WWmD[:,-1,:,:]
        
        # calculate dissipation
        viscdamp=np.mean(U*UmD+V*VmD+WWmD,0) # time average
        B_visc=np.nansum(self.rhoref*viscdamp*self.dV)
        plt.pcolormesh(self.yc[:,0],self.zc,np.sum(self.rhoref*viscdamp*self.dV,2).squeeze(),cmap='hot')
        plt.colorbar()
        return B_visc

    def energetics_visc(self,iteration=np.inf,runpath=None,runnum=np.inf):
        U=self.U(iteration=iteration,runpath=runpath,runnum=runnum)
        V=self.V(iteration=iteration,runpath=runpath,runnum=runnum)
        W=self.W(iteration=iteration,runpath=runpath,runnum=runnum)
        # fill in values
        U[:,self.ktop-1,np.arange(0,self.ny),np.arange(0,self.nx)]=U.isel(z=self.ktop,lat=np.arange(0,self.ny),lon=np.arange(0,self.nx))
        V[:,self.ktop-1,np.arange(0,self.ny),np.arange(0,self.nx)]=V.isel(z=self.ktop,lat=np.arange(0,self.ny),lon=np.arange(0,self.nx))
        W[:,self.ktop-1,np.arange(0,self.ny),np.arange(0,self.nx)]=W.isel(z=self.ktop,lat=np.arange(0,self.ny),lon=np.arange(0,self.nx))
        U[:,np.maximum(self.ktop-2,0),np.arange(0,self.ny),np.arange(0,self.nx)]=U.isel(z=self.ktop,lat=np.arange(0,self.ny),lon=np.arange(0,self.nx))
        V[:,np.maximum(self.ktop-2,0),np.arange(0,self.ny),np.arange(0,self.nx)]=V.isel(z=self.ktop,lat=np.arange(0,self.ny),lon=np.arange(0,self.nx))
        W[:,np.maximum(self.ktop-2,0),np.arange(0,self.ny),np.arange(0,self.nx)]=W.isel(z=self.ktop,lat=np.arange(0,self.ny),lon=np.arange(0,self.nx))
        # shift W grid
        W[:,:-1,:,:]=0.5*(W[:,1:,:,:].values+W[:,:-1,:,:].values)
        W[:,-1,:,:]=0.5*W[:,-1,:,:].values

        # compute vorticity
        vortz=self.vorticity_z(iteration=iteration,runpath=runpath,runnum=runnum,U=U,V=V)
        vorty=self.vorticity_y(iteration=iteration,runpath=runpath,runnum=runnum,U=U,W=W)
        vortx=self.vorticity_x(iteration=iteration,runpath=runpath,runnum=runnum,V=V,W=W)

        # assign coord
        U['lat']=np.radians(U.lat)
        V['lat']=np.radians(V.lat)
        W['lat']=np.radians(W.lat)
        U['lon']=np.radians(U.lon)
        V['lon']=np.radians(V.lon)
        W['lon']=np.radians(W.lon)

        # compute derivatives
        r=(self.a+self.zc)[np.newaxis,:,np.newaxis,np.newaxis]
        Uzz=(U.differentiate('z')*r**2).differentiate('z')/r**2
        Vzz=(V.differentiate('z')*r**2).differentiate('z')/r**2
        Wzz=(W.differentiate('z')*r**2).differentiate('z')/r**2
        clat=np.cos(np.radians(self.yc))[np.newaxis,np.newaxis,:,:]
        Uyy=(U.differentiate('lat')*clat).differentiate('lat')/clat/r**2
        Vyy=(V.differentiate('lat')*clat).differentiate('lat')/clat/r**2
        Wyy=(W.differentiate('lat')*clat).differentiate('lat')/clat/r**2
        if self.nx>1:
            Uxx=(U.differentiate('lon')).differentiate('lon')/clat**2/r**2
            Vxx=(V.differentiate('lon')).differentiate('lon')/clat**2/r**2
            Wxx=(W.differentiate('lon')).differentiate('lon')/clat**2/r**2
        else:
            Uxx=0
            Vxx=0
            Wxx=0
        # smooth
        Uxx=gaussian_filter(Uxx.values,sigma=2)
        Uyy=gaussian_filter(Uyy.values,sigma=2)
        Uzz=gaussian_filter(Uzz.values,sigma=2)
        Vxx=gaussian_filter(Vxx.values,sigma=2)
        Vyy=gaussian_filter(Vyy.values,sigma=2)
        Vzz=gaussian_filter(Vzz.values,sigma=2)
        Wxx=gaussian_filter(Wxx.values,sigma=2)
        Wyy=gaussian_filter(Wyy.values,sigma=2)
        Wzz=gaussian_filter(Wzz.values,sigma=2)

        #viscdamp=np.mean(-self.ah*(vortz.values**2+vortx.values**2*1+vorty.values**2*1)+1*(self.av-self.ah)*(U.values*Uzz+V.values*Vzz+W.values*Wzz),0)
        #viscdamp=-self.ah*np.mean(vortz.values**2,0)
        #viscdamp=-self.ah*np.mean(vortz.values**2+vortx.values**2*(self.av/self.ah)+vorty.values**2*(self.av/self.ah),0)
        viscdamp=np.mean(self.ah*(U.values*Uxx+V.values*Vxx+W.values*Wxx+U.values*Uyy+V.values*Vyy+W.values*Wyy)+self.av*(U.values*Uzz+V.values*Vzz+W.values*Wzz),0)
        dV=self.rac.values[np.newaxis,:,:]*self.deepfac2c[:,np.newaxis,np.newaxis]*self.widthfac*self.dzc[:,np.newaxis,np.newaxis]
        #(U*Uzz+V*Vzz+W*Wzz).mean('Time').squeeze().plot(vmin=-2e-12,vmax=2e-12,cmap='RdBu_r') 
        hfacc1=self.hfacc
        #hfacc1[1:]=hfacc1[0:-1].values
        B_visc=np.nansum(self.rhoref*viscdamp*dV*(hfacc1))
        plt.pcolormesh(self.yc[:,0],self.zc,np.sum(self.rhoref*viscdamp*dV*(hfacc1),2).squeeze(),cmap='hot')
        plt.colorbar()
        return B_visc

    def vorticity_x(self,iteration=np.inf,runpath=None,runnum=np.inf,V=None,W=None,zfac=1):
        if V is None:
            V=self.V(iteration=iteration,runpath=runpath,runnum=runnum)
        if W is None:
            W=self.W(iteration=iteration,runpath=runpath,runnum=runnum)
        clat=np.cos(np.radians(self.yc[:,0]))
        r=self.a+self.zc
        Vr=V*r[:,np.newaxis,np.newaxis]
        Vr['z']=r
        Vz=Vr.differentiate('z')/r[:,np.newaxis,np.newaxis]
        vortx=-Vz*zfac
        if self.nx>1:
            W['lat']=np.radians(self.yc[:,0])
            Wy=W.differentiate('lat')/r[:,np.newaxis,np.newaxis]
            vortx=vortx-Wy.values

        return vortx

    def vorticity_y(self,iteration=np.inf,runpath=None,runnum=np.inf,U=None,W=None,zfac=1):
        if U is None:
            U=self.U(iteration=iteration,runpath=runpath,runnum=runnum)
        if W is None:
            W=self.W(iteration=iteration,runpath=runpath,runnum=runnum)
        r=self.a+self.zc
        Ur=U*r[:,np.newaxis,np.newaxis]
        Ur['z']=r
        Uz=Ur.differentiate('z')/r[:,np.newaxis,np.newaxis]
        vorty=Uz*zfac
        if self.nx>1:
            W['lon']=np.radians(self.xc[0,:])
            Wx=W.differentiate('lon')/np.cos(np.radians(self.yc))/r[:,np.newaxis,np.newaxis]
            vorty=vorty-Wx.values
        return vorty

    def vorticity_z(self,iteration=np.inf,runpath=None,runnum=np.inf,U=None,V=None):
        if U is None:
            U=self.U(iteration=iteration,runpath=runpath,runnum=runnum)
        if V is None:
            V=self.V(iteration=iteration,runpath=runpath,runnum=runnum)
        Uc=U*np.cos(np.radians(self.yc))
        Uc['lat']=np.sin(np.radians(self.yc[:,0]))
        r=self.a+self.zc
        Uy=Uc.differentiate('lat')/r[:,np.newaxis,np.newaxis]
        vortz=-Uy
        if self.nx>1:
            V['lon']=np.radians(self.xc[0,:])
            Vx=V.differentiate('lon')/np.cos(np.radians(self.yc))/r[:,np.newaxis,np.newaxis]
            vortz=vortz+Vx.values
        return vortz

    def nonlinear_alpha(self,iteration=np.inf,runpath=None,runnum=np.inf,onlyz=True, T=None, S=None):
        if T is None:
            T=self.T(iteration=iteration,runpath=runpath,runnum=runnum).mean('Time').mean('lon')
        if S is None:
            S=self.S(iteration=iteration,runpath=runpath,runnum=runnum).mean('Time').mean('lon')
        if onlyz and 'lat' in T.dims:
            T=(T*self.wgt[np.newaxis,:]).mean('lat')/(self.mask.mean('lon')*self.wgt[np.newaxis,:]).mean('lat')
        if onlyz and 'lat' in S.dims:
            S=(S*self.wgt[np.newaxis,:]).mean('lat')/(self.mask.mean('lon')*self.wgt[np.newaxis,:]).mean('lat')
        P=self.rhoref*np.cumsum(self.gprofile*self.dzc)/1e4 # dbar
        if not onlyz:
            P=P[:,np.newaxis]
        alpha=gsw.alpha(S,gsw.CT_from_pt(S,T),P)*self.mask
        return alpha.mean('lon')

    def nonlinear_beta(self,iteration=np.inf,runpath=None,runnum=np.inf,onlyz=True, T=None, S=None):
        if T is None:
            T=self.T(iteration=iteration,runpath=runpath,runnum=runnum).mean('Time').mean('lon')
        if S is None:
            S=self.S(iteration=iteration,runpath=runpath,runnum=runnum).mean('Time').mean('lon')
        if onlyz and 'lat' in T.dims:
            T=(T*self.wgt[np.newaxis,:]).mean('lat')/(self.mask.mean('lon')*self.wgt[np.newaxis,:]).mean('lat')
        if onlyz and 'lat' in S.dims:
            S=(S*self.wgt[np.newaxis,:]).mean('lat')/(self.mask.mean('lon')*self.wgt[np.newaxis,:]).mean('lat')
        P=self.rhoref*np.cumsum(self.gprofile*self.dzc)/1e4 # dbar
        if not onlyz:
            P=P[:,np.newaxis]
        beta=gsw.beta(S,gsw.CT_from_pt(S,T),P)*self.mask
        return beta.mean('lon')

    def update(self,iteration=np.inf,runpath=None,runnum=np.inf):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        self.__init__(name=self.name,runpath=runpath,update=True)
#        hfacc1=mu.rdmds(self.path + runpath+ 'hFacC')
#        self.hfacc = xr.DataArray( hfacc1, dims=('z', 'lat', 'lon'),
#            coords=dict(z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
#        hfacw1=mu.rdmds(self.path + runpath+ 'hFacW')
#        self.hfacw = xr.DataArray( hfacw1, dims=('z', 'lat', 'lon'),
#            coords=dict(z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
#        self.mask=self.hfacc*0+1
#        self.mask=self.mask.where(self.hfacc>0.0)
#        self.totarea_zc=np.sum(self.rac.values[np.newaxis,:,:]*self.deepfac2c[:,np.newaxis,np.newaxis]*self.widthfac*self.hfacc,axis=1).squeeze()
    
    def get_Tminmax(self,iteration=np.inf,runpath=None,runnum=np.inf,constrain=None):
        T=self.T(iteration=iteration,runpath=runpath,runnum=runnum).mean('Time').mean('lon')
        if constrain is 'NH':
            T=T.where(T.lat>=0,drop=True)
        if constrain is 'SH':
            T=T.where(T.lat<=0,drop=True)

        self.Tmin=float(T.min().values)
        self.Tmax=float(T.max().values)
        return (self.Tmin,self.Tmax)

    def get_Sminmax(self,iteration=np.inf,runpath=None,runnum=np.inf,constrain=None):
        S=self.S(iteration=iteration,runpath=runpath,runnum=runnum).mean('Time').mean('lon')
        if constrain is 'NH':
            S=S.where(S.lat>=0,drop=True)
        if constrain is 'SH':
            S=S.where(S.lat<=0,drop=True)
        self.Smin=float(S.min().values)
        self.Smax=float(S.max().values)
        return (self.Smin,self.Smax)


    def get_non_dims(self, final_tsteps=50):
        return [
            self.D / 1000, self.Q0,
            np.mean(self.ra().values[-final_tsteps:]),
            np.mean(self.rastar().values[-final_tsteps:]), self.rastarq,
            self.E, self.Ev,
            np.mean(self.nustar().values[-final_tsteps:]),
            np.mean(self.nu().values[-final_tsteps:]),
            np.mean(self.roloc().values[-final_tsteps:])
        ]

    def get_deltaT(self):
        T = self.T(iteration=np.inf)
        vmeanT = (T * self.rac*self.mask).sum(('lat', 'lon')) / (self.rac*self.mask).sum(('lat', 'lon'))
        vmeanT=vmeanT.dropna(dim='z')
        deltaT = vmeanT.isel(z=-1) - vmeanT.isel(z=0)
        return deltaT

    def rastar(self):
        return self.alpha * self.g / (self.omega**2 *
                                      self.D) * self.get_deltaT()

    def ra(self):
        return self.alpha * self.g * self.D**3 / (
            self.av**2) * self.get_deltaT()  # assumed Pr = 1

    def nustar(self):
        return self.rastarq / self.rastar()

    def nu(self):
        return self.rastarq / self.ra() / self.E**3  # assumed Pr = 1

    def roloc(self):
        return self.ra()**(5 / 4) * self.E**2

    def get_iter(self,Time,unit='year'):
        if unit=='day':
            Time=Time/360
        iter1=np.int(np.floor(self.iter0*self.dt/self.dtoutput+Time*360*86400/self.dtoutput))*self.dtoutput/self.dt
        return iter1

    def get_runpath(self,iteration,runnum=None):
        it=iteration
        if runnum is None:
            if self.itermax>0:
                if it>=0 and it<=self.itermax:
                    runnum=np.mod(it,self.eachiter)
                    runnum=it-runnum-int(runnum==0 and it!=0)*self.eachiter
                else:
                    runnum=self.itermax-self.eachiter
                    if it>=0:
                        it=self.itermax
                runpath='run{}/'.format(int(runnum))
            else:
                if not np.isnan(it):
                    it=np.inf
                runpath=''
        else:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
            it=None
        return it,runpath;

    def get_var(self,var,iteration=np.nan,sample=1,dimnum=3,rec=None):
        # tuple iteration will cause the script to read all available records and select those within the range
        #it,runpath=self.get_runpath(iteration)
        if type(iteration)==tuple:
            it=np.nan
        else:
            it=iteration
        Vardata,its,meta = mu.rdmds(self.path + var, it,returnmeta=True,rec=rec)
        its=np.array(its)
        if ((sample>1 and np.isnan(it)) or  ((type(iteration)==tuple) and len(iteration)>1)):
            itsid=np.arange(len(its))
            if (type(iteration)==tuple) and len(iteration)>1:
                ibeg=iteration[0]
                iend=iteration[1]
            else:
                ibeg=0
                iend=np.inf
            ikeep=np.array(np.where((np.mod(itsid,sample)==0) & (its>=ibeg) & (its<=iend))).squeeze()
            if dimnum==3:
                Vardata=Vardata[ikeep,:,:,:]
            if dimnum==2:
                Vardata=Vardata[ikeep,:,:]
                
            its=its[ikeep]
        Time = np.array(its)*self.dt/86400

        if np.ndim(Vardata)==dimnum:
            Vardata=np.reshape(Vardata,np.concatenate([[1],np.array(Vardata.shape)]))

        if dimnum==3:
            #Vardata=np.multiply(Vardata,np.tile(self.mask,[np.size(its),1,1,1]))
            Varxdata= xr.DataArray(
                Vardata,
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=Time, z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
            Varxdata=Varxdata.where(self.hfacc>0.0)

        if dimnum==2:
            Varxdata= xr.DataArray(
                Vardata,
                dims=('Time', 'lat', 'lon'),
                coords=dict(
                    Time=Time, lat=self.yc[:, 0], lon=self.xc[0, :]))

        return Varxdata
    
    def U(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        Udata=self.get_var(runpath+'U',iteration=iteration,sample=sample)
        Udata.name='U'
        Udata.attrs['units']='m/s'
        Udata.attrs['showname']='U (m/s)'
        Udata.attrs['longshowname']='Zonal flow speed (m/s)'
        return Udata

    def T(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        Tdata=self.get_var(runpath+'T',iteration=iteration,sample=sample)
        Tdata.name='T'
        Tdata.attrs['units']='degC'
        Tdata.attrs['showname']='T (degC)'
        Tdata.attrs['longshowname']='Temperature (degC)'
        return Tdata

    def V(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        Vdata=self.get_var(runpath+'V',iteration=iteration,sample=sample)
        Vdata.name='V'
        Vdata.attrs['units']='m/s'
        Vdata.attrs['showname']='V (m/s)'
        Vdata.attrs['longshowname']='Meridional flow speed (m/s)'
        return Vdata

    def W(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        Wdata=self.get_var(runpath+'W',iteration=iteration,sample=sample)
        Wdata.name='W'
        Wdata.attrs['units']='m/s'
        Wdata.attrs['showname']='W (m/s)'
        Wdata.attrs['longshowname']='Vertical flow speed (m/s)'
        return Wdata

    def SHIRshel(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        SHIRsheldata=self.get_var(runpath+'surfDiag',rec=7,iteration=iteration,sample=sample,dimnum=2)
        SHIRsheldata.attrs['units']='m'
        SHIRsheldata.attrs['showname']='SHIRshel (m)'
        SHIRsheldata.attrs['longshowname']='ice shelf thickness (m)'
        return SHIRsheldata

    def SHI_mass(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        SHI_massdata=self.get_var(runpath+'surfDiag',rec=6,iteration=iteration,sample=sample,dimnum=2)
        SHI_massdata.name='SHI_mass'
        SHI_massdata.attrs['units']='kg/m2'
        SHI_massdata.attrs['showname']='SHI_mass (kg/m2)'
        SHI_massdata.attrs['longshowname']='ice shelf mass (kg/m2)'
        return SHI_massdata

    def PNH(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        PNHdata=self.get_var(runpath+'PNH',iteration=iteration,sample=sample)
        PNHdata.values=self.rhoref*PNHdata.values
        PNHdata.name='PNH'
        PNHdata.attrs['units']='Pa'
        PNHdata.attrs['showname']='PNH (Pa)'
        PNHdata.attrs['longshowname']='Non-hydrostatic pressure (Pa)'
        return PNHdata

    def PH(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        PHdata=self.get_var(runpath+'PH',iteration=iteration,sample=sample)
        PHdata.values=self.rhoref*PHdata.values
        PHdata.name='PH'
        PHdata.attrs['units']='Pa'
        PHdata.attrs['showname']='PH (Pa)'
        PHdata.attrs['longshowname']='Hydrostatic pressure (Pa)'
        return PHdata

    def q(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None,inmassunit=False):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        iteration=iteration
        qdata=self.get_var(runpath+'SHICE_fwFlux',iteration=iteration,sample=sample,dimnum=2)*self.colmask[np.newaxis,:]
        if inmassunit:
            qdata.name='q'
            qdata.attrs['units']='kg/s/m^2'
            qdata.attrs['showname']='q (kg/s/m^2)'
            qdata.attrs['longshowname']='freezing rate (kg/s/m^2)'
        else:
            qdata=qdata/1000*(86400*360*1e6)/1000
            qdata.name='q'
            qdata.attrs['units']='km/Myr'
            qdata.attrs['showname']='q (km/Myr)'
            qdata.attrs['longshowname']='freezing rate (km/Myr)'
        return qdata

    def eta(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None,inmassunit=False):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        iteration=iteration
        etadata=self.get_var(runpath+'Eta',iteration=iteration,sample=sample,dimnum=2)*self.colmask[np.newaxis,:]
        etadata.name='eta'
        etadata.attrs['units']='m'
        etadata.attrs['showname']='eta (m)'
        etadata.attrs['longshowname']='surface ele (m)'
        return etadata

    def S(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        Sdata=self.get_var(runpath+'S',iteration=iteration,sample=sample)
        Sdata.name='S'
        Sdata.attrs['units']='psu'
        Sdata.attrs['showname']='S (psu)'
        Sdata.attrs['longshowname']='Salinity (psu)'
        return Sdata

    def Psi_tot(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None):
        def fill_array(arr):
            arr_=arr.copy()
            arr_[:,0]=arr[:,1].copy()
            arr_[:,-1]=arr[:,-2].copy()
            for ilat in range(self.ny):
                imiss=np.where(np.isnan(arr_[:,ilat]))[0]
                arr_[imiss,ilat]=arr_[len(imiss),ilat].copy()
            return arr_
        PsiGM=self.Psi_GM(iteration=iteration,runnum=runnum,sample=sample,runpath=runpath)
        PsiEu=self.Psi(iteration=iteration,runnum=runnum,sample=sample,runpath=runpath)
        Psieddy=self.Psi_eddy(iteration=iteration,runnum=runnum,sample=sample,runpath=runpath)
        Psi=PsiEu
        Psi.values=PsiGM.values+PsiEu.values
        if self.nx>1:
            Psi=Psi.mean('Time')+Psieddy
        missmask=np.where(np.isnan(Psi),0,1)
        Psi=fill_array(Psi)
        Psismt=Psi.copy()
        Psismt.values=gaussian_filter(Psi,sigma=[3,9],truncate=3,mode='nearest')
        Psi=fill_array(Psismt)*missmask
        Psi[:,0]=0
        Psi[:,-1]=0
        Psi.attrs['units']=PsiEu.attrs['units']
        Psi.attrs['longshowname']='Total meridional streamfunction ('+Psi.attrs['units']+')'
        Psi.attrs['showname']='Psi_tot ('+Psi.attrs['units']+')'
        return Psi

    def Psi_eddy(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None):
        if self.nx>1:
            if runpath==None:
                if np.isnan(runnum):
                    runpath='run*/'
                elif runnum<0:
                    runpath=''
                else:
                    if self.runmax>=0:
                        runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                        runpath='run{}/'.format(runid)
                    else:
                        runpath=''
            def fill_array(arr,fill='const'):
                arr_=arr.copy()
                if fill=='const':
                    arr_[:,0]=arr[:,1].copy()
                    arr_[:,-1]=arr[:,-2].copy()
                elif fill=='grad':
                    arr_[:,0]=2*arr[:,1]-arr[:,2]
                    arr_[:,-1]=2*arr[:,-2]-arr[:,-3]

                for ilat in range(self.ny):
                    imiss=np.where(np.isnan(arr_[:,ilat]))[0]
                    #print(ilat,imiss)
                    arr_[imiss,ilat]=arr_[len(imiss),ilat].copy()
                return arr_
            #ADVr_TH=self.get_var(runpath+'flxDiag',iteration=iteration,sample=sample,rec=1).mean('Time').mean('lon')
            _,its,meta=mu.rdmds(self.path + runpath + 'dynDiag', np.inf,returnmeta=True,rec=0)
            iWVELMASS=meta['fldlist'].index('WVELMASS')
            WVELMASS=self.get_var(runpath+'dynDiag',iteration=iteration,sample=sample,rec=iWVELMASS).mean('Time').mean('lon')
            iWTHMASS=meta['fldlist'].index('WTHMASS')
            WTHMASS=self.get_var(runpath+'dynDiag',iteration=iteration,sample=sample,rec=iWTHMASS).mean('Time').mean('lon')
            iTHETA=meta['fldlist'].index('THETA')
            T=self.get_var(runpath+'dynDiag',iteration=iteration,sample=sample,rec=iTHETA).mean('Time').mean('lon')
            wt=WTHMASS-WVELMASS*T
            #wt=(ADVr_TH)/self.rac.mean('lon')/self.deepfac2c[:,np.newaxis]-W*T
            wt=fill_array(wt,fill='const')
            T=fill_array(T, fill='grad')
            Ty=T.differentiate('lat')
            if self.spheric:
                Ty=Ty.values/self.a*180.0/np.pi
            Tymag=np.nanmean(np.abs(Ty).flatten())
            Ty=np.maximum(np.abs(Ty),Tymag/2)*np.sign(Ty+1e-15)
            if self.spheric:
                coslat=(np.cos(np.radians(self.yc[:,0]))).squeeze()
                Psieddy=-2*np.pi*self.a*coslat[np.newaxis,:]*wt/(Ty)*self.rhoref
            else:
                Psieddy=-wt/(Ty)*self.rhoref
            Psieddysmt=Psieddy.copy()
            Psieddysmt.values=gaussian_filter(Psieddy,sigma=4,truncate=2,mode='nearest')
            Psieddy=fill_array(Psieddysmt)
        else:
            Psieddy=np.zeros(np.shape(self.hfacc.values[np.newaxis,:,:,:]))
            Psieddy=xr.DataArray( Psieddy, dims=('z', 'lat'),
                coords=dict(z=self.zc, lat=self.yc[:, 0]))

        Psieddy.name='Psi_eddy'
        if self.spheric:
            Psieddy.attrs['units']='kg/s'
            Psieddy.attrs['showname']='Psi_eddy (kg/s)'
            Psieddy.attrs['longshowname']='eddy-induced meridional stream function (kg/s)'
        else:
            Psieddy.attrs['units']='kg/s/m'
            Psieddy.attrs['showname']='Psi_eddy (kg/s/m)'
            Psieddy.attrs['longshowname']='eddy-induced meridional stream function (kg/s/m)'
        
        return Psieddy

    def Psi_GM(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None):
        if self.GM:
            if runpath==None:
                if np.isnan(runnum):
                    runpath='run*/'
                elif runnum<0:
                    runpath=''
                else:
                    if self.runmax>=0:
                        runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                        runpath='run{}/'.format(runid)
                    else:
                        runpath=''
            _,its,meta=mu.rdmds(self.path + runpath + 'GMDiag', np.inf,returnmeta=True,rec=0)
            iGM_PsiY=meta['fldlist'].index('GM_PsiY')
            PsiGM=self.get_var(runpath+'GMDiag',iteration=iteration,sample=sample,rec=iGM_PsiY)
            PsiGM=PsiGM.where(self.hfacc>0.0)
            tmp=PsiGM.values.squeeze()
            tmp[np.isnan(tmp)]=0
            PsiGM.values=filters.gaussian_filter(tmp,[0.8,0.8],mode='constant',cval=0)[np.newaxis,:,:,np.newaxis]
            PsiGM=PsiGM.where(self.hfacc>0.0)
            PsiGM=np.mean(PsiGM,axis=-1)
            if self.spheric:
                PsiGM=PsiGM*self.rhoref*(2*np.pi*self.a*np.cos(np.radians(self.yg[np.newaxis,np.newaxis,:-1,0])))
                PsiGM.attrs['units']='kg/s'
                PsiGM.attrs['showname']='Psi_GM (kg/s)'
                PsiGM.attrs['longshowname']='GM-induced meridional stream function (kg/s)'
            else:
                PsiGM=PsiGM*self.rhoref
                PsiGM.attrs['units']='kg/s/m'
                PsiGM.attrs['showname']='Psi (kg/s/m)'
                PsiGM.attrs['longshowname']='GM-induced meridional stream function (kg/s/m)'

            PsiGM=PsiGM.where(self.hfacc>0.0)
            PsiGM.name='Psi_GM'
        else:
            PsiGM=np.zeros(np.array([1,self.nz,self.ny]))
            PsiGM=xr.DataArray( PsiGM, dims=('Time','z', 'lat'),
                    coords=dict(Time=[0],z=self.zc, lat=self.yc[:, 0]))
            PsiGM.name='Psi_GM'
            PsiGM.attrs['units']='kg/s'
            PsiGM.attrs['showname']='Psi_GM (kg/s)'
            PsiGM.attrs['longshowname']='GM-induced meridional stream function (kg/s)'

        return PsiGM


    def Rho(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None,S=None,T=None):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        if self.eos.lower() == 'linear':
            if not 'rray' in type(S).__name__:
                S=self.get_var(runpath+'S',iteration=iteration,sample=sample)
            if not 'rray' in type(T).__name__:
                T=self.get_var(runpath+'T',iteration=iteration,sample=sample)
            if 'Time' in S.dims:
                T['Time']=S.Time
            density=self.rhoref*(self.beta*(S-self.Sref)-self.alpha*(T-self.Tref))
        else:
            rhotmp=self.get_var(runpath+'oceDiag',iteration=iteration,sample=sample,rec=1)
            rho_bot=rhotmp.isel(z=-1).values
            rho_bot=rho_bot[:,np.newaxis,:,:]
            drhodr=self.get_var(runpath+'oceDiag',iteration=iteration,sample=sample,rec=0)
            drhodrmid=(drhodr.values[:,1:,:,:]+drhodr.values[:,:-1,:,:])/2.0
            intdrhodr=np.cumsum(drhodrmid[:,::-1,:,:]*self.dzf[np.newaxis,::-1,np.newaxis,np.newaxis],axis=1)
            intdrhodr=np.concatenate((np.zeros_like(rho_bot),intdrhodr),axis=1)[:,::-1,:,:]
            density=rhotmp
            rho=intdrhodr+rho_bot
            density.values=rho-np.nanmean(rho)

###            density=density-self.rhoprof0[:,np.newaxis,np.newaxis]
###            density=density-density.mean()
###            rho_zvar=density.mean('lon').mean('lat').mean('Time')
###            rho_zvar=rho_zvar.fillna(rho_zvar.dropna(dim='z')[0])
###            rho_zvarsmt=savgol_filter(rho_zvar.values,int(len(self.zc)/4),1)
###            rho_corr=np.array(rho_zvarsmt-rho_zvar)
###            #rho_corr=np.where(np.abs(rho_corr)>1e-4,1e-4*np.sign(rho_corr),rho_corr)
###            density=density+rho_corr[np.newaxis,:,np.newaxis,np.newaxis]
        density.name='Rho'
        density.attrs['units']='kg/m3'
        density.attrs['showname']='Rho (kg/m3)'
        density.attrs['longshowname']='density (kg/m3)'
        return density

    def dRhodR(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        drhodr=self.get_var(runpath+'oceDiag',iteration=iteration,sample=sample,rec=0)
        drhodr.name='Rho_z'
        drhodr.attrs['units']='kg/m3/m'
        drhodr.attrs['showname']='Rho_z (kg/m4)'
        drhodr.attrs['longshowname']='Stratification (kg/m4)'
        return drhodr


    def Psi_(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None,V=None):
        if not 'rray' in type(V).__name__:
            if runpath==None:
                if np.isnan(runnum):
                    runpath='run*/'
                elif runnum<0:
                    runpath=''
                else:
                    if self.runmax>=0:
                        runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                        runpath='run{}/'.format(runid)
                    else:
                        runpath=''
            V=self.get_var(runpath+'V',iteration=iteration,sample=sample)

        if not 'lon' in V.dims:
            vdr=V*self.hfacc.mean(dim='lon')*self.drf*2*np.pi*self.a*np.cos(np.radians(self.yg[:-1]))*self.rhoref
            if self.deepatm:
                vdr=vdr*self.deepfacc[:, np.newaxis]
        else:
            #vdr=V*self.hfacc*self.drf*2*np.pi*self.a*np.cos(np.radians(self.yg[:-1]))*self.rhoref
            vdr=V*self.hfacc*self.drf*self.dxg*self.rhoref
            if self.deepatm:
                vdr=vdr*self.deepfacc[:, np.newaxis,np.newaxis]
        Psidata=-(vdr.reindex(z=vdr.z[::-1])).cumsum(dim='z').reindex(z=vdr.z).sum('lon')*(360/(self.xg[0,-1]-self.xg[0,0]))
        Psidata=Psidata.where(self.hfacc>0.0)
        Psidata.name='Psi'
        Psidata.attrs['units']='kg/s'
        Psidata.attrs['showname']='Psi (kg/s)'
        Psidata.attrs['longshowname']='Meridional stream function (kg/s)'
        return Psidata 

    def Psi(self,iteration=np.nan,runnum=np.inf,sample=1,runpath=None,V=None):
        if not 'rray' in type(V).__name__:
            if runpath==None:
                if np.isnan(runnum):
                    runpath='run*/'
                elif runnum<0:
                    runpath=''
                else:
                    if self.runmax>=0:
                        runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                        runpath='run{}/'.format(runid)
                    else:
                        runpath=''
            V=self.get_var(runpath+'V',iteration=iteration,sample=sample)
        if not 'lon' in V.dims:
            vdr=V*self.hfacc.mean(dim='lon')*self.drf*self.rhoref
            if self.deepatm:
                vdr=vdr*self.deepfacc[:, np.newaxis]
        else:
            #vdr=V*self.hfacc*self.drf*2*np.pi*self.a*np.cos(np.radians(self.yg[:-1]))*self.rhoref
            vdr=V*self.hfacc*self.drf*self.rhoref
            if self.deepatm:
                vdr=vdr*self.deepfacc[:, np.newaxis,np.newaxis]
            vdr=vdr.mean('lon')

        if self.spheric:
#            print((2*np.pi*self.a*np.cos(np.radians(self.yg[np.newaxis,np.newaxis,:-1,0]))))
#            print(self.dxg.sum(1)*(360/(self.xg[0,-1]-self.xg[0,0])))
            Psidata=-(vdr.reindex(z=vdr.z[::-1])).cumsum(dim='z').reindex(z=vdr.z)*(2*np.pi*self.a*np.cos(np.radians(self.yg[np.newaxis,np.newaxis,:-1,0])))
            Psidata.attrs['units']='kg/s'
            Psidata.attrs['showname']='Psi (kg/s)'
            Psidata.attrs['longshowname']='Meridional stream function (kg/s)'
        else:
            Psidata=-(vdr.reindex(z=vdr.z[::-1])).cumsum(dim='z').reindex(z=vdr.z)
            Psidata.attrs['units']='kg/s/m'
            Psidata.attrs['showname']='Psi (kg/s/m)'
            Psidata.attrs['longshowname']='Meridional stream function (kg/s/m)'
        Psidata=Psidata.where(self.hfacc.mean('lon')>0.0)
        tmp=Psidata.values.squeeze()
        tmp[np.isnan(tmp)]=0

        if len(Psidata.Time)>1:
            Psidata.values=filters.gaussian_filter(tmp,[0,0.8,0.8],mode='constant',cval=0)
        else:
            Psidata.values=filters.gaussian_filter(tmp,[0.8,0.8],mode='constant',cval=0)[np.newaxis,:,:]
        Psidata=Psidata.where(self.hfacc.mean('lon')>0.0)
        Psidata.name='Psi'
        return Psidata

    def TRC(self,trc=1,iteration=np.nan,runnum=np.inf,runpath=None,sample=1):
        #only one tracer can be read in this function, starting from 1
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        TRCdata=self.get_var(runpath+'PTRACER{:02}'.format(trc),iteration=iteration,sample=sample)
        TRCdata.name='TRC{:02}'.format(trc)
        TRCdata.attrs['units']='1'
        TRCdata.attrs['showname']='TRC{:02}'.format(trc)
        TRCdata.attrs['longshowname']='TRC '+self.trcnames[trc-1]
        return TRCdata

    def TRC01(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1):
        TRC01data=self.TRC(trc=1,iteration=iteration,runnum=runnum,runpath=runpath,sample=sample)
        return TRC01data
    def TRC02(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1):
        TRC02data=self.TRC(trc=2,iteration=iteration,runnum=runnum,runpath=runpath,sample=sample)
        return TRC02data
    def TRC03(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1):
        TRC03data=self.TRC(trc=3,iteration=iteration,runnum=runnum,runpath=runpath,sample=sample)
        return TRC03data
    def TRC04(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1):
        TRC04data=self.TRC(trc=4,iteration=iteration,runnum=runnum,runpath=runpath,sample=sample)
        return TRC04data
    def TRC05(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1):
        TRC05data=self.TRC(trc=5,iteration=iteration,runnum=runnum,runpath=runpath,sample=sample)
        return TRC05data

    def get(self,varlist,iteration=np.nan,runnum=np.inf,runpathstr=None,sample=1,dimmeths={},shift=False):
        dic={}
        codetemplate1="""dic['varname']=self.varname(iteration=it_,runnum=rn_,runpath=rp_,sample=sp_)"""
        for var in varlist:
            print('reading {}'.format(var))
            code=codetemplate1.replace('varname',var).replace('it_',str(iteration)).replace('rn_',str(runnum)).replace('nan','np.nan').replace('inf','np.inf').replace('sp_',str(sample)).replace('rp_',str(runpathstr))
            if var=='heatflux_surf':
                code=code.replace(")",",shift="+str(shift)+")")
            if var=='Rho' and 'T' in dic:
                code=code.replace(")",",T=dic['T'])")
            if var=='Rho' and 'S' in dic:
                code=code.replace(")",",S=dic['S'])")
            if var=='Psi' and 'V' in dic:
                code=code.replace(")",",V=dic['V'])")

            exec(code)
            dic[var]=dimcomp(dic[var],key=var,dimmeths=dimmeths)

        return dic


    def monitor(self,iteration=np.inf,runpathstr=None,runnum=np.inf,pltarray=(3,2),wdirs=[0,0,0,1,1,1],pltvar=['T','S','Rho','U','W','Psi'],dimmeths={'Time':'mean','lon':'mean','lat':None,'z':None},figsize=(12,12),projection=None,pltcontour={},flip=True,pltdir='F',labelorientation=None,xlabelpad=0,ylabelpad=20,labelpad=0.04,labelshrink=1,labelaspect=20,linetype=['k-']*20,sharex=True,sharey=True,xlims=None,ylims=None,vminmax={},cms={},savefig=0,sample=3,alwaysrefmean=False, returnimhandle=False, refcase=None,refrunpathstr=None,refrunnum=None,sharec=False,cutcm=[0,1],symmetrify=False,shift=False,xfactor=1,interp=False):
        d=self.get(pltvar,runpathstr=runpathstr,iteration=iteration,sample=sample,runnum=runnum,shift=shift,dimmeths=dimmeths)
        if refcase is not None:
            dr=self.get(pltvar,runpathstr="'../data_"+refcase+"/"+refrunpathstr[1:],iteration=iteration,sample=sample,shift=shift,dimmeths=dimmeths)
            for var in pltvar:
                d[var]=d[var]-dr[var].values
        if refrunnum is not None:
            dr=self.get(pltvar,iteration=iteration,sample=sample,runnum=refrunnum,shift=shift,dimmeths=dimmeths)
            for var in pltvar:
                d[var]=d[var]-dr[var].values

        #its=d[pltvar[0]].Time/self.dt*86400.0
        #its=np.round(its.values)
        #print('its={}-{}'.format(its[0],its[-1]))
        fig,axf,da,figname,imh=my2dplt(self,d,pltarray,wdirs=wdirs,dimmeths=dimmeths,pltvar=pltvar,figsize=figsize,projection=projection,pltcontour=pltcontour,flip=flip,pltdir=pltdir,labelorientation=labelorientation,xlabelpad=xlabelpad,ylabelpad=ylabelpad,labelpad=labelpad,labelshrink=labelshrink,labelaspect=labelaspect,linetype=linetype,sharex=sharex,sharey=sharey,xlims=xlims,ylims=ylims,vminmax=vminmax,cms=cms,savefig=savefig,alwaysrefmean=alwaysrefmean,sharec=sharec,cutcm=cutcm,symmetrify=symmetrify,xfactor=xfactor,interp=interp)
        if returnimhandle:
            return fig,axf,da,imh
        else:
            return fig,axf,da

    def diffusivity_eq(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=5,doplt=False):
        T_bar=self.T(iteration=iteration,runnum=runnum,runpath=runpath,sample=sample).mean(dim=('lon','Time')).values
        dTv=T_bar[-5,:]-T_bar[4,:]
        H_surf=self.heatflux_surf(iteration=iteration,runnum=runnum,runpath=runpath,sample=sample).mean(dim=('Time','lon')).values
        kappav_eq=H_surf/dTv/(self.cp*self.rhoref)*self.D
        
        if doplt:
            with plt.style.context(('labelsize15')):
                fig, ax = plt.subplots(1,1)
                im=ax.plot(self.yc[:,0],kappav_eq,color='black')
                ax.set_title('vertical diffusivity',loc='right')
                plt.xlabel('lat')

        return(kappav_eq)

    def thermalwind_uz(self,iteration=np.nan,runnum=np.inf,runpath=None, **kwargs):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''

        u=self.U(iteration=iteration,runnum=runnum,runpath=runpath).mean('lon')
        tmean=np.mean(u.Time.values)
        u=u.mean('Time')
        u_z=xr.concat([u.isel(z=0), u, u.isel(z=-1)], "z").transpose()
        du = u_z.diff("z").values
        drc = mu.rdmds(self.path + runpath+ 'DRC').squeeze()
        dudr = du / drc[:, np.newaxis]
        dudr = 0.5 * (dudr[:-1] + dudr[1:])
        u_y=xr.concat([u.isel(lat=1), u], "lat")
        dtheta = np.radians(np.diff(self.yc[:, 0])[0])
        dy = (self.rsphere + self.zc) * dtheta
        durdtheta = u_y.diff('lat') / dy[:, np.newaxis]
        theta = np.radians(self.yc[:, 0])
        dudz = -np.sin(theta) * dudr + np.cos(theta) * durdtheta
        omega_dudz = (2 * self.omega)*dudz.values

        omega_dudz = xr.DataArray(
                omega_dudz[np.newaxis,:,:], dims=('Time','z', 'lat'), coords=dict(Time=tmean[np.newaxis],z=self.zc, lat=self.yc[:, 0]))
        return omega_dudz
    
    def thermalwind_by(self,iteration=np.nan,runnum=np.inf,runpath=None, **kwargs):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''

        rho=self.Rho(iteration=iteration,runnum=runnum,runpath=runpath).mean('lon')
        b=self.g*(rho-self.rhoref)/self.rhoref
        tmean=np.mean(b.Time.values)
        b=b.mean('Time')
        dtheta = np.radians(np.diff(self.yc[:, 0])[0])
        r = (self.rsphere + self.zc)
        dy = (self.rsphere + self.zc) * dtheta

        b_=xr.concat([b[:,1],b[:,1], b[:,1:-1],b[:,-2]], "lat")
        dbdy = b_.diff('lat') / dy[:, np.newaxis]
        dbdy = dbdy.values

        dbdy = xr.DataArray(
                dbdy[np.newaxis,:,:], dims=('Time','z', 'lat'), coords=dict(Time=tmean[np.newaxis],z=self.zc, lat=self.yc[:, 0]))
        return dbdy


    def reconstruct_b(self,iteration=np.nan,runnum=np.inf,runpath=None, **kwargs):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''

        # calculate dudz 
        u=self.U(iteration=iteration,runnum=runnum,runpath=runpath).mean('lon')
        tmean=np.mean(u.Time.values)
        u=u.mean('Time')
        u_z=xr.concat([u.isel(z=0), u, u.isel(z=-1)], "z").transpose()
        du = u_z.diff("z").values
        drc = mu.rdmds(self.path + runpath+ 'DRC').squeeze()
        dudr = du / drc[:, np.newaxis]
        dudr = 0.5 * (dudr[:-1] + dudr[1:])
        u_y=xr.concat([u.isel(lat=1), u], "lat")
        dtheta = np.radians(np.diff(self.yc[:, 0])[0])
        dy = (self.rsphere + self.zc) * dtheta
        durdtheta = u_y.diff('lat') / dy[:, np.newaxis]
        theta = np.radians(self.yc[:, 0])
        dudz = -np.sin(theta) * dudr + np.cos(theta) * durdtheta
        # integrate in y
        rhs = (2 * self.omega)*dudz
        r = (self.rsphere + self.zc)
        rhs *= r[:, np.newaxis]
        reconst_b = np.zeros_like(rhs)
        reconst_b[:, theta > 0] = np.cumsum(rhs[:, theta > 0], axis=1) * dtheta
        reconst_b[:, theta < 0] = -np.cumsum(
            rhs[:, theta < 0][:, ::-1], axis=1)[:, ::-1] * dtheta
        b=self.Rho(iteration=iteration,runnum=runnum,runpath=runpath).mean('lon').mean('Time')
        ymin=np.min(np.abs(self.yc[:,0]))+0.1
        beq=b.sel(lat=slice(-ymin, ymin)).mean('lat').values
        reconst_b=reconst_b + beq[:,np.newaxis]

        reconst_b = xr.DataArray(
                reconst_b[np.newaxis,:,:], dims=('Time','z', 'lat'), coords=dict(Time=tmean[np.newaxis],z=self.zc, lat=self.yc[:, 0]))

        return reconst_b

    def heatbudget(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=5):

        H_surf=self.heatflux_surf(iteration=iteration,runnum=runnum,runpath=runpath,sample=sample).mean(dim=('Time','lon')).values
        H_ocny=self.bulk_heatflux_y(iteration=iteration,runnum=runnum,runpath=runpath,sample=sample).values
        
        # devide the whole globe into 3 equal pieces
        #boollow=(abs(self.yc[:,0])<=self.thetaeq3)
        #boolhigh=((self.yc[:,0])>self.thetaeq3)
        #indbndN=np.argmin(abs(self.yc[:,0]-self.thetaeq3))
        #indbndS=np.argmin(abs(-self.yc[:,0]-self.thetaeq3))

        # devide the NH into 2 equal pieces
        boollow=(self.yc[:,0]<=self.thetaeq2) & (self.yc[:,0]>=0)
        boolhigh=(self.yc[:,0]>self.thetaeq2)
        indbndN=np.argmin(abs(self.yc[:,0]-self.thetaeq2))
        indbndS=np.argmin(abs(-self.yc[:,0]-self.thetaeq2))

        zonalbandarea=np.sum(self.rac*((self.a-(self.rhoice/self.rhoref)*self.Hice)/self.a)**(2*self.deepatm),axis=1)*(360/(self.xg[0,-1]-self.xg[0,0]))
        arealow=np.sum(zonalbandarea[boollow])
        areahigh=np.sum(zonalbandarea[boolhigh])
        zonalbandarea_inter=np.sum(self.rac*(self.ainter/self.a)**(2*self.deepatm),axis=1)*(360/(self.xg[0,-1]-self.xg[0,0]))
        arealow_inter=np.sum(zonalbandarea_inter[boollow])
        areahigh_inter=np.sum(zonalbandarea_inter[boolhigh])

        scale=1e3 # W to mW
        dheat=dict(H_ocny=H_ocny[[indbndS,indbndN]]*1e6/np.array([arealow_inter,areahigh_inter])*scale,
                H_surf=[np.sum(H_surf[boollow]*zonalbandarea[boollow])/arealow_inter*scale,np.sum(H_surf[boolhigh]*zonalbandarea[boolhigh])/areahigh_inter*scale],
                Qbot=[np.sum(self.Q[boollow]*zonalbandarea_inter[boollow])/arealow_inter*scale,np.sum(self.Q[boolhigh]*zonalbandarea_inter[boolhigh])/areahigh_inter*scale],
                Htide=[np.sum(self.Htide[boollow]*zonalbandarea[boollow])/arealow_inter*scale,np.sum(self.Htide[boolhigh]*zonalbandarea[boolhigh])/areahigh_inter*scale],
                Hcond=[np.sum(self.Hcond[boollow]*zonalbandarea[boollow])/arealow_inter*scale,np.sum(self.Hcond[boolhigh]*zonalbandarea[boolhigh])/areahigh_inter*scale]
                )
        return dheat

    def bulk_heatflux_r_adv(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1,doplt=False):
        advflux_r=self.bulk_heatflux_r(iteration=iteration,runnum=runnum,runpath=runpath,sample=sample,doplt=doplt,component='adv')
        return advflux_r

    def bulk_heatflux_r_dif(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1,doplt=False):
        difflux_r=self.bulk_heatflux_r(iteration=iteration,runnum=runnum,runpath=runpath,sample=sample,doplt=doplt,component='dif')
        return difflux_r


    def bulk_heatflux_r(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1,component='advdif',doplt=False,tmean=True,return3d=False):
        if runpath is None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''

        if 'adv' in component:
            print('reading ADVr_TH...')
            Vardata,its,meta = mu.rdmds(self.path + runpath + 'flxDiag', iteration,returnmeta=True,rec=1)
            print(str(its))
            print('its={}-{}'.format(its[0],its[-1]))
            if len(its)==1:
                Vardata=Vardata[np.newaxis,:,:,:]
            Time = np.array(its)*self.dt/86400
            ADVr_TH= xr.DataArray(
                Vardata,
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=Time, z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        else:
            ADVr_TH=xr.DataArray(
                    np.zeros((1,len(self.zc),len(self.yc[:,0]),len(self.xc[0,:]))),
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=[0], z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        if 'dif' in component:
            print('reading DFrE_TH...')
            Vardata,its,meta = mu.rdmds(self.path + runpath + 'flxDiag', iteration,returnmeta=True,rec=3)
            print('its={}-{}'.format(its[0],its[-1]))
            if len(its)==1:
                Vardata=Vardata[np.newaxis,:,:,:]
            Time = np.array(its)*self.dt/86400
            DFr_TH= xr.DataArray(
                Vardata,
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=Time, z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
            print('reading DFrI_TH...')
            Vardata,its,meta = mu.rdmds(self.path + runpath + 'flxDiag', iteration,returnmeta=True,rec=4)
            if len(its)==1:
                Vardata=Vardata[np.newaxis,:,:,:]
            DFr_TH.data=DFr_TH.data+Vardata

        else:
            DFr_TH=xr.DataArray(
                    np.zeros((1,len(self.zc),len(self.yc[:,0]),len(self.xc[0,:]))),
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=[0], z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))

        FLXr_TH=DFr_TH+ADVr_TH

        Heatflux_r=FLXr_TH*self.rhoref*self.cp/self.dxg[np.newaxis,np.newaxis,:,:]/self.dyg[np.newaxis,np.newaxis,:,:]/self.deepfac2f[np.newaxis,1:,np.newaxis,np.newaxis]
        Heatflux_r.attrs['units']='W/m2'
        Heatflux_r.attrs['showname']='F_Heat_r (W/m2)'
        Heatflux_r.attrs['longshowname']='Vertical Heat Flux (W/m2)'
        Heatflux_r.name='Heatflux_r'
        if doplt:
            with plt.style.context(('labelsize15')):
                fig, ax = plt.subplots(1,1)
                colors = plt.cm.jet(np.linspace(0,1,self.nz))
                for iz in range(self.nz):
                    im=ax.plot(Heatflux_r.lat,Heatflux_r.mean('Time').mean('lon').squeeze()[iz,:],color=colors[iz])
                ax.plot(Heatflux_r.lat,Heatflux_r.lat*0.,'k--')
                ax.set_title(Heatflux_r.longshowname,loc='right')
                plt.xlabel('lat')
        if tmean:
            Heatflux_r=Heatflux_r.mean('Time')
            DFr_TH=DFr_TH.mean('Time')
            ADVr_TH=ADVr_TH.mean('Time')

        if return3d:
            #DFr_TH=DFr_TH*self.cp*self.rhoref
            #ADVr_TH=ADVr_TH*self.cp*self.rhoref
            return (Heatflux_r,ADVr_TH,DFr_TH)
        else:
            return Heatflux_r

    def bulk_heatflux_y_adv(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1,doplt=False,layers=None,intx=True):
        advflux_y=self.bulk_heatflux_y(iteration=iteration,runnum=runnum,runpath=runpath,sample=sample,doplt=doplt,component='adv',layers=layers,intx=intx)
        return advflux_y

    def bulk_heatflux_y_dif(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1,doplt=False,layers=None,intx=True):
        difflux_y=self.bulk_heatflux_y(iteration=iteration,runnum=runnum,runpath=runpath,sample=sample,doplt=doplt,component='dif',layers=layers,intx=intx)
        return difflux_y

    def bulk_heatflux_y(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1,component='advdif',doplt=False,layers=None,tmean=True,intx=True,return3d=False):
        if runpath is None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''

        if 'adv' in component:
            print('reading ADVy_TH...')
            Vardata,its,meta = mu.rdmds(self.path + runpath + 'flxDiag', iteration,returnmeta=True,rec=0)
            print(str(its))
            print('its={}-{}'.format(its[0],its[-1]))
            if len(its)==1:
                Vardata=Vardata[np.newaxis,:,:,:]
            Time = np.array(its)*self.dt/86400
            ADVy_TH= xr.DataArray(
                Vardata,
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=Time, z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        else:
            ADVy_TH=xr.DataArray(
                    np.zeros((1,len(self.zc),len(self.yc[:,0]),len(self.xc[0,:]))),
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=[0], z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        if 'dif' in component:
            print('reading DFyE_TH...')
            Vardata,its,meta = mu.rdmds(self.path + runpath + 'flxDiag', iteration,returnmeta=True,rec=2)
            print('its={}-{}'.format(its[0],its[-1]))
            if len(its)==1:
                Vardata=Vardata[np.newaxis,:,:,:]
            Time = np.array(its)*self.dt/86400
            DFyE_TH= xr.DataArray(
                Vardata,
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=Time, z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        else:
            DFyE_TH=xr.DataArray(
                    np.zeros((1,len(self.zc),len(self.yc[:,0]),len(self.xc[0,:]))),
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=[0], z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))

        if tmean:
            DFyE_TH=DFyE_TH.mean(dim='Time')
            ADVy_TH=ADVy_TH.mean(dim='Time')
        FLXy_TH=DFyE_TH+ADVy_TH

        if layers is None:
            layers=range(len(self.zc))
        if intx:
            if self.spheric:
                Heatflux_y=FLXy_TH.isel(z=layers).sum(dim=['lon','z']).squeeze()*self.rhoref*self.cp/1e6/((self.xg[0,-1]-self.xg[0,0])/360)
                Heatflux_y.attrs['units']='MW'
                Heatflux_y.attrs['showname']='F_Heat_y (MW)'
                Heatflux_y.attrs['longshowname']='Meridional Heat Flux (MW)'
            else:
                Heatflux_y=FLXy_TH.isel(z=layers).sum(dim=['z']).mean(dim='lon').squeeze()*self.rhoref*self.cp/self.dxg[0,0]/1e6
                Heatflux_y.attrs['units']='MW/m'
                Heatflux_y.attrs['showname']='F_Heat_y (MW/m)'
                Heatflux_y.attrs['longshowname']='Meridional Heat Flux (MW/m)'
            
            if doplt:
                with plt.style.context(('labelsize15')):
                    fig, ax = plt.subplots(1,1)
                    im=ax.plot(Heatflux_y.lat,Heatflux_y,color='black')
                    ax.plot(Heatflux_y.lat,Heatflux_y*0.,'k--')
                    ax.set_title(Heatflux_y.longshowname,loc='right')
                    plt.xlabel('lat')
        else:
            Heatflux_y=FLXy_TH.isel(z=layers).sum(dim=['z']).squeeze()*self.rhoref*self.cp/self.dxg[0,0]/1e6
            Heatflux_y.attrs['units']='MW/m'
            Heatflux_y.attrs['showname']='F_Heat_y (MW/m)'
            Heatflux_y.attrs['longshowname']='Meridional Heat Flux (MW/m)'

        
        Heatflux_y.name='Heatflux_y'
        if return3d:
            #DFyE_TH=DFyE_TH*self.cp*self.rhoref
            #ADVy_TH=ADVy_TH*self.cp*self.rhoref
            return (Heatflux_y,ADVy_TH,DFyE_TH)
        else:
            return Heatflux_y
        

    def bulk_saltflux_y_adv(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1,doplt=False,layers=None,intx=True):
        advflux_y=self.bulk_saltflux_y(iteration=iteration,runnum=runnum,runpath=runpath,sample=sample,doplt=doplt,component='adv',layers=layers,intx=intx)
        return advflux_y

    def bulk_saltflux_y_dif(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1,doplt=False,layers=None,intx=True):
        difflux_y=self.bulk_saltflux_y(iteration=iteration,runnum=runnum,runpath=runpath,sample=sample,doplt=doplt,component='dif',layers=layers,intx=intx)
        return difflux_y

    def bulk_saltflux_y(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1,component='advdif',doplt=False,layers=None,intx=True):
        if runpath is None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''

        if 'adv' in component:
            print('reading ADVy_SLT...')
            Vardata,its,meta = mu.rdmds(self.path + runpath + 'flxDiag', iteration,returnmeta=True,rec=5)
            print('its={}-{}'.format(its[0],its[-1]))
            if len(its)==1:
                Vardata=Vardata[np.newaxis,:,:,:]
            Time = np.array(its)*self.dt/86400
            ADVy_SLT= xr.DataArray(
                Vardata,
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=Time, z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        else:
            ADVy_SLT=xr.DataArray(
                    np.zeros((1,len(self.zc),len(self.yc[:,0]),len(self.xc[0,:]))),
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=[0], z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        if 'dif' in component:
            print('reading DFyE_SLT...')
            Vardata,its,meta = mu.rdmds(self.path + runpath + 'flxDiag', iteration,returnmeta=True,rec=7)
            print('its={}-{}'.format(its[0],its[-1]))
            if len(its)==1:
                Vardata=Vardata[np.newaxis,:,:,:]
            Time = np.array(its)*self.dt/86400
            DFyE_SLT= xr.DataArray(
                Vardata,
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=Time, z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        else:
            DFyE_SLT=xr.DataArray(
                    np.zeros((1,len(self.zc),len(self.yc[:,0]),len(self.xc[0,:]))),
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=[0], z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))

        FLXy_SLT=DFyE_SLT.mean(dim='Time')+ADVy_SLT.mean(dim='Time')

        if layers is None:
            layers=range(len(self.zc))
        if intx:
            if self.spheric:
                Saltflux_y=FLXy_SLT.isel(z=layers).sum(dim=['lon','z']).squeeze()*self.rhoref/1e3/((self.xg[0,-1]-self.xg[0,0])/360)
                Saltflux_y.attrs['units']='kg/s'
                Saltflux_y.attrs['showname']='F_Salt_y (kg/s)'
                Saltflux_y.attrs['longshowname']='Meridional Salt Flux (kg/s)'
            else:
                Saltflux_y=FLXy_SLT.isel(z=layers).sum(dim=['z']).squeeze()*self.rhoref/self.dxg[0,0]/1e3
                Saltflux_y.attrs['units']='kg/s/m'
                Saltflux_y.attrs['showname']='F_Salt_y (kg/s/m)'
                Saltflux_y.attrs['longshowname']='Meridional Salt Flux (kg/s/m)'
            if doplt:
                with plt.style.context(('labelsize15')):
                    fig, ax = plt.subplots(1,1)
                    im=ax.plot(Saltflux_y.lat,Saltflux_y,color='black')
                    ax.set_title(Saltflux_y.longshowname,loc='right')
                    plt.xlabel('lat')
        else:
            Saltflux_y=FLXy_SLT.isel(z=layers).sum(dim=['z']).squeeze()*self.rhoref/self.dxg[0,0]/1e3
            Saltflux_y.attrs['units']='kg/s/m'
            Saltflux_y.attrs['showname']='F_Salt_y (kg/s/m)'
            Saltflux_y.attrs['longshowname']='Meridional Salt Flux (kg/s/m)'

        Saltflux_y.name='Saltflux_y'

        return Saltflux_y


    def bulk_heatflux_x_adv(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1,doplt=False,layers=None,inty=True):
        advflux_x=self.bulk_heatflux_x(iteration=iteration,runnum=runnum,runpath=runpath,sample=sample,doplt=doplt,component='adv',layers=layers,inty=inty)
        return advflux_x

    def bulk_heatflux_x_dif(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1,doplt=False,layers=None,inty=True):
        difflux_x=self.bulk_heatflux_x(iteration=iteration,runnum=runnum,runpath=runpath,sample=sample,doplt=doplt,component='dif',layers=layers,inty=inty)
        return difflux_x

    def bulk_heatflux_x(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1,component='advdif',doplt=False,layers=None,tmean=True,inty=True):
        if runpath is None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''

        if 'adv' in component:
            print('reading ADVx_TH...')
            Vardata,its,meta = mu.rdmds(self.path + runpath + 'flxDiag', iteration,returnmeta=True,rec=10)
            print(str(its))
            print('its={}-{}'.format(its[0],its[-1]))
            if len(its)==1:
                Vardata=Vardata[np.newaxis,:,:,:]
            Time = np.array(its)*self.dt/86400
            ADVy_TH= xr.DataArray(
                Vardata,
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=Time, z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        else:
            ADVy_TH=xr.DataArray(
                    np.zeros((1,len(self.zc),len(self.yc[:,0]),len(self.xc[0,:]))),
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=[0], z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        if 'dif' in component:
            print('reading DFyE_TH...')
            Vardata,its,meta = mu.rdmds(self.path + runpath + 'flxDiag', iteration,returnmeta=True,rec=12)
            print('its={}-{}'.format(its[0],its[-1]))
            if len(its)==1:
                Vardata=Vardata[np.newaxis,:,:,:]
            Time = np.array(its)*self.dt/86400
            DFyE_TH= xr.DataArray(
                Vardata,
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=Time, z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        else:
            DFyE_TH=xr.DataArray(
                    np.zeros((1,len(self.zc),len(self.yc[:,0]),len(self.xc[0,:]))),
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=[0], z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))

        if tmean:
            FLXy_TH=DFyE_TH.mean(dim='Time')+ADVy_TH.mean(dim='Time')
        else:
            FLXy_TH=DFyE_TH+ADVy_TH

        if layers is None:
            layers=range(len(self.zc))
        if inty:
            if self.spheric:
                Heatflux_x=FLXy_TH.isel(z=layers).sum(dim=['lat','z']).squeeze()*self.rhoref*self.cp/1e6
                Heatflux_x.attrs['units']='MW'
                Heatflux_x.attrs['showname']='F_Heat_x (MW)'
                Heatflux_x.attrs['longshowname']='Zonal Heat Flux (MW)'
            else:
                Heatflux_x=FLXy_TH.isel(z=layers).sum(dim=['z']).mean(dim='lat').squeeze()*self.rhoref*self.cp/self.dyg[0,0]/1e6
                Heatflux_x.attrs['units']='MW/m'
                Heatflux_x.attrs['showname']='F_Heat_x (MW/m)'
                Heatflux_x.attrs['longshowname']='Zonal Heat Flux (MW/m)'
            
            if doplt:
                with plt.style.context(('labelsize15')):
                    fig, ax = plt.subplots(1,1)
                    im=ax.plot(Heatflux_x.lat,Heatflux_x,color='black')
                    ax.plot(Heatflux_x.lat,Heatflux_x*0.,'k--')
                    ax.set_title(Heatflux_x.longshowname,loc='left')
                    plt.xlabel('lat')
        else:
            Heatflux_x=FLXy_TH.isel(z=layers).sum(dim=['z']).squeeze()*self.rhoref*self.cp/self.dyg[0,0]/1e6
            Heatflux_x.attrs['units']='MW/m'
            Heatflux_x.attrs['showname']='F_Heat_x (MW/m)'
            Heatflux_x.attrs['longshowname']='Zonal Heat Flux (MW/m)'

        
        Heatflux_x.name='Heatflux_x'
        return Heatflux_x
        

    def bulk_saltflux_x_adv(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1,doplt=False,layers=None,inty=True):
        advflux_x=self.bulk_saltflux_x(iteration=iteration,runnum=runnum,runpath=runpath,sample=sample,doplt=doplt,component='adv',layers=layers,inty=inty)
        return advflux_x

    def bulk_saltflux_x_dif(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1,doplt=False,layers=None,inty=True):
        difflux_x=self.bulk_saltflux_x(iteration=iteration,runnum=runnum,runpath=runpath,sample=sample,doplt=doplt,component='dif',layers=layers,inty=inty)
        return difflux_x

    def bulk_saltflux_x(self,iteration=np.nan,runnum=np.inf,runpath=None,sample=1,component='advdif',doplt=False,layers=None,inty=True):
        if runpath is None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''

        if 'adv' in component:
            print('reading ADVx_SLT...')
            Vardata,its,meta = mu.rdmds(self.path + runpath + 'flxDiag', iteration,returnmeta=True,rec=11)
            print('its={}-{}'.format(its[0],its[-1]))
            if len(its)==1:
                Vardata=Vardata[np.newaxis,:,:,:]
            Time = np.array(its)*self.dt/86400
            ADVy_SLT= xr.DataArray(
                Vardata,
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=Time, z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        else:
            ADVy_SLT=xr.DataArray(
                    np.zeros((1,len(self.zc),len(self.yc[:,0]),len(self.xc[0,:]))),
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=[0], z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        if 'dif' in component:
            print('reading DFyE_SLT...')
            Vardata,its,meta = mu.rdmds(self.path + runpath + 'flxDiag', iteration,returnmeta=True,rec=13)
            print('its={}-{}'.format(its[0],its[-1]))
            if len(its)==1:
                Vardata=Vardata[np.newaxis,:,:,:]
            Time = np.array(its)*self.dt/86400
            DFyE_SLT= xr.DataArray(
                Vardata,
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=Time, z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        else:
            DFyE_SLT=xr.DataArray(
                    np.zeros((1,len(self.zc),len(self.yc[:,0]),len(self.xc[0,:]))),
                dims=('Time', 'z', 'lat', 'lon'),
                coords=dict(
                    Time=[0], z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))

        FLXy_SLT=DFyE_SLT.mean(dim='Time')+ADVy_SLT.mean(dim='Time')

        if layers is None:
            layers=range(len(self.zc))
        if inty:
            if self.spheric:
                Saltflux_x=FLXy_SLT.isel(z=layers).sum(dim=['lat','z']).squeeze()*self.rhoref/1e3
                Saltflux_x.attrs['units']='kg/s'
                Saltflux_x.attrs['showname']='F_Salt_x (kg/s)'
                Saltflux_x.attrs['longshowname']='Zonal Salt Flux (kg/s)'
            else:
                Saltflux_x=FLXy_SLT.isel(z=layers).sum(dim=['z']).mean(dim='lat').squeeze()*self.rhoref/self.dxg[0,0]/1e3
                Saltflux_x.attrs['units']='kg/s/m'
                Saltflux_x.attrs['showname']='F_Salt_x (kg/s/m)'
                Saltflux_x.attrs['longshowname']='Zonal Salt Flux (kg/s/m)'
            if doplt:
                with plt.style.context(('labelsize15')):
                    fig, ax = plt.subplots(1,1)
                    im=ax.plot(Saltflux_x.lat,Saltflux_x,color='black')
                    ax.set_title(Saltflux_x.longshowname,loc='left')
                    plt.xlabel('lat')
        else:
            Saltflux_x=FLXy_SLT.isel(z=layers).sum(dim=['z']).squeeze()*self.rhoref/self.dxg[0,0]/1e3
            Saltflux_x.attrs['units']='kg/s/m'
            Saltflux_x.attrs['showname']='F_Salt_x (kg/s/m)'
            Saltflux_x.attrs['longshowname']='Zonal Salt Flux (kg/s/m)'

        Saltflux_x.name='Saltflux_x'

        return Saltflux_x

    def xmom_budget(self,runpath=None,runnum=np.inf,iteration=np.nan,sample=1,shift=True,angular=True,title='',makeplot=True,viscdiag=False):
        taux_surf=self.xmomflux_surf(runpath=runpath,runnum=runnum,iteration=iteration,sample=sample,shift=shift).mean('Time')
        taux_bot=self.xmomflux_bot(runpath=runpath,runnum=runnum,iteration=iteration,sample=sample,shift=shift).mean('Time')
        #Hice=self.real_Hice(runpath=runpath,runnum=runnum,iteration=iteration).mean('Time') # 2D
        Hice=self.Hice
        Hunder=Hice*self.rhoice/self.rhoref
        coslat=np.cos(np.radians(self.yc))
        ## geometry based on desired topo
        #eta=self.eta(runpath=runpath,runnum=runnum,iteration=iteration).mean('Time') # 2D
        #acore=self.rsphere-self.Htot
        #asurf=self.rsphere-Hunder+eta

        ## geometry based on modelled topo
        zcfull=self.hfacc
        zcfull.values=np.tile(self.zc[:,np.newaxis,np.newaxis],[1,self.ny,self.nx])
        zsurf=zcfull.isel(z=self.ktop,lat=np.arange(0,self.ny),lon=np.arange(0,self.nx))
        asurf=zsurf+self.rsphere
        acore=self.rsphere+self.zc[-1]

        deepfac_bot=(acore/self.rsphere)**2
        deepfac_surf=(asurf/self.rsphere)**2
        maskoutTC=np.where(acore<asurf*coslat,1,0)
        taux_botmean=taux_bot.mean('lon')
        taux_surfmean=taux_surf.mean('lon')
        if makeplot:
            fig, ax = plt.subplots(1,1)
            ax.plot(self.yc[:,0],taux_bot.mean('lon'),'k-')
            ax.plot(self.yc[:,0],taux_surf.mean('lon'),'k--')
            plt.legend(['bot','surf'],fontsize=13)
            ax.ticklabel_format(axis='y',scilimits=(-1,2))
        if angular:
            tottaux_bot=np.sum(taux_bot*acore*coslat*self.rac*deepfac_bot)*self.widthfac
            tottaux_surf=np.sum(taux_surf*asurf*coslat*self.rac*deepfac_surf)*self.widthfac
            tottaux_surfout=np.sum(taux_surf*asurf*coslat*self.rac*deepfac_surf*maskoutTC)*self.widthfac
            tottaux_surfin=tottaux_surf-tottaux_surfout
        else:
            tottaux_bot=np.sum(taux_bot*self.rac*deepfac_bot)*self.widthfac
            tottaux_surf=np.sum(taux_surf*self.rac*deepfac_surf)*self.widthfac
            tottaux_surfout=np.sum(taux_surf*self.rac*deepfac_surf*maskoutTC)*self.widthfac
            tottaux_surfin=tottaux_surf-tottaux_surfout

        if viscdiag:
            tottaux_visc=self.xmom_visc(runpath=runpath,runnum=runnum,iteration=iteration)

        if angular:
            if viscdiag:
                dic={'xmom torque':['bottom','shell-highlat','shell-lowlat','total','visc'],'(PN*m)':[tottaux_bot.values/1e15, tottaux_surfin.values/1e15, tottaux_surfout.values/1e15, (tottaux_bot+tottaux_surf).values/1e15, tottaux_visc/1e15]}
            else:
                dic={'xmom torque':['bottom','shell-highlat','shell-lowlat','total'],'(PN*m)':[tottaux_bot.values/1e15, tottaux_surfin.values/1e15, tottaux_surfout.values/1e15, (tottaux_bot+tottaux_surf).values/1e15]}
        else:
            dic={'xmom torque':['bottom','shell-highlat','shell-lowlat','total'],'(GN)':[tottaux_bot.values/1e9, tottaux_surfin.values/1e9, tottaux_surfout.values/1e9, (tottaux_bot+tottaux_surf).values/1e9]}
        df = pd.DataFrame(dic)
        df.columns.name = title
        display(df)
        return dic,taux_botmean,taux_surfmean


    def ymomflux_surf(self,runpath=None,runnum=np.inf,iteration=np.nan,sample=1,shift=True):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        #shift=false
        if shift and np.isfinite(iteration):
            iteration=iteration+self.iteroutput/2
        tauydata=self.get_var(runpath+'shiDiag',iteration=iteration,sample=sample,dimnum=2,rec=6)
        tauydata=-tauydata
        if shift and np.isfinite(iteration):
            iteration=iteration-self.iteroutput/2
        tauydata.name='core-ocean meri mom exchange'
        tauydata.attrs['units']='n/m^2'
        tauydata.attrs['showname']='drag_bot (n/m^2)'
        tauydata.attrs['longshowname']='core-ocean meri mom exchange (n/m^2)'
        return tauydata

    def ymomflux_bot(self,runpath=None,runnum=np.inf,iteration=np.nan,sample=1,shift=True):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        #shift=false
        if shift and np.isfinite(iteration):
            iteration=iteration+self.iteroutput/2
        tauydata=self.get_var(runpath+'shiDiag',iteration=iteration,sample=sample,dimnum=2,rec=8)
        tauydata=-tauydata
        if shift and np.isfinite(iteration):
            iteration=iteration-self.iteroutput/2
        tauydata.name='core-ocean meri mom exchange'
        tauydata.attrs['units']='N/m^2'
        tauydata.attrs['showname']='drag_bot (N/m^2)'
        tauydata.attrs['longshowname']='core-ocean meri mom exchange (N/m^2)'
        return tauydata


    def xmomflux_surf(self,runpath=None,runnum=np.inf,iteration=np.nan,sample=1,shift=True):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        #shift=false
        if shift and np.isfinite(iteration):
            iteration=iteration+self.iteroutput/2
        tauxdata=self.get_var(runpath+'shiDiag',iteration=iteration,sample=sample,dimnum=2,rec=5)
        tauxdata=-tauxdata
        if shift and np.isfinite(iteration):
            iteration=iteration-self.iteroutput/2
        tauxdata.name='ice-ocean zonal mom exchange'
        tauxdata.attrs['units']='N/m^2'
        tauxdata.attrs['showname']='drag_ice (N/m^2)'
        tauxdata.attrs['longshowname']='ice-ocean zonal mom exchange (N/m^2)'
        return tauxdata

    def xmomflux_bot(self,runpath=None,runnum=np.inf,iteration=np.nan,sample=1,shift=True):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        #shift=false
        if shift and np.isfinite(iteration):
            iteration=iteration+self.iteroutput/2
        tauxdata=self.get_var(runpath+'shiDiag',iteration=iteration,sample=sample,dimnum=2,rec=7)
        tauxdata=-tauxdata
        if shift and np.isfinite(iteration):
            iteration=iteration-self.iteroutput/2
        tauxdata.name='core-ocean zonal mom exchange'
        tauxdata.attrs['units']='N/m^2'
        tauxdata.attrs['showname']='drag_bot (N/m^2)'
        tauxdata.attrs['longshowname']='core-ocean zonal mom exchange (N/m^2)'
        return tauxdata


    def xmom_visc(self,runpath=None,runnum=np.inf,iteration=np.nan,sample=1,return3D=False,shift=True):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        #shift=false
        if shift and np.isfinite(iteration):
            iteration=iteration+self.iteroutput/2
        _,its,meta=mu.rdmds(self.path + runpath + 'momDiag', np.inf,returnmeta=True,rec=0)
        iVISCy_Um=meta['fldlist'].index('VISCy_Um')
        FUy=self.get_var(runpath+'momDiag',iteration=iteration,sample=sample,rec=iVISCy_Um).mean('Time') #VISCy_Um
        iVISrE_Um=meta['fldlist'].index('VISrE_Um')
        iVISrI_Um=meta['fldlist'].index('VISrI_Um')
        FUr=self.get_var(runpath+'momDiag',iteration=iteration,sample=sample,rec=iVISrE_Um).mean('Time')+self.get_var(runpath+'momDiag',iteration=iteration,sample=sample,rec=iVISrI_Um).mean('Time') #VISrE_Um, VISrI_Um
        UtenddV=(FUy.values[:-1,:-1,:]-FUy.values[:-1,1:,:])+(FUr.values[1:,:-1,:]-FUr.values[:-1,:-1,:])
        Rhoanom=self.get_var(runpath+'oceDiag',iteration=iteration,sample=sample,rec=1).mean('Time')
        AMtend=UtenddV*(self.rhoref+Rhoanom[:-1,:-1,:])*self.widthfac*(self.deepfacc[:,np.newaxis,np.newaxis]*self.rsphere*np.cos(np.radians(self.yc))[np.newaxis,:,:])[:-1,:-1,:]
        TOTAMtend=np.sum(AMtend)
        if return3D:
            AMtend= xr.DataArray( AMtend, dims=('z', 'lat', 'lon'),
                    coords=dict(z=self.zc[:-1], lat=self.yc[:-1, 0], lon=self.xc[0, :]))
            return TOTAMtend.values,AMtend
        else:
            return TOTAMtend.values


    def heatflux_surf(self,runpath=None,runnum=np.inf,iteration=np.nan,sample=1,removeq=True,q=None,shift=True):
        if runpath==None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        #shift=False
        if shift and np.isfinite(iteration):
            iteration=iteration+self.iteroutput/2
        ForcTdata=self.get_var(runpath+'shiDiag',iteration=iteration,sample=sample,dimnum=2,rec=1)
        ForcTdata=-ForcTdata
        if shift and np.isfinite(iteration):
            iteration=iteration-self.iteroutput/2

        if removeq:
            if q is None:
                q=self.q(runpath=runpath,runnum=runnum,iteration=iteration,sample=sample,inmassunit=True)
            else:
                if q.units=='km/Myr':
                    q=q*1e3/1e6/360/86400*self.rhoref
            if 'Time' in q.dims:
                q=q.mean('Time').values
                q=q[np.newaxis,:]
            ForcTdata=ForcTdata*(self.rhoref*self.gammaT/(self.rhoref*self.gammaT-q))*self.colmask[np.newaxis,:]

        ForcTdata.name='Ice-ocean heat exchange'
        ForcTdata.attrs['units']='W/m^2'
        ForcTdata.attrs['showname']='H_ocn (W/m^2)'
        ForcTdata.attrs['longshowname']='Ice-ocean heat exchange (W/m^2)'
        return ForcTdata
        


    def get_tprime(self, lats, latn, runpath=None, runnum=np.inf, iteration=np.inf,sample=3):
        T = self.T(iteration=iteration,runnum=runnum,sample=sample,runpath=runpath)
        T=T.mean(dim='Time')
        T = T.sel(lat=slice(lats, latn))
        self.Tprime = np.sqrt(
            np.mean(T**2, axis=(1, 2)) - np.mean(T, axis=(
                1, 2)) * np.mean(T, axis=(1, 2)))
        return self.Tprime


    def get_uprime(self, lats, latn, runpath=None, runnum=np.inf, iteration=np.inf,sample=3):
        u = self.U(iteration=iteration,runnum=runnum,sample=sample,runpath=runpath)
        v = self.V(iteration=iteration,runnum=runnum,sample=sample,runpath=runpath)
        u=u.mean(dim='Time')
        v=v.mean(dim='Time')
        u = u.sel(lat=slice(lats, latn))
        v = v.sel(lat=slice(lats, latn))
        u_meanofsquares = np.mean(u**2 + v**2, axis=(1, 2))
        u_mean = np.mean(u, axis=(1, 2))
        v_mean = np.mean(v, axis=(1, 2))
        u_squareofmeans = u_mean**2 + v_mean**2
        self.Uprime = np.sqrt(u_meanofsquares - u_squareofmeans)
        return self.Uprime

    def umom(self, iteration=np.nan, runnum=np.inf, runpath=None):
        if runpath is None:
            if np.isnan(runnum):
                runpath='run*/'
            elif runnum<0:
                runpath=''
            else:
                if self.runmax>=0:
                    runid=int(np.minimum(runnum,self.runmax)*self.eachiter)
                    runpath='run{}/'.format(runid)
                else:
                    runpath=''
        umcori,its,meta = mu.rdmds(self.path + runpath + "momDiag", iteration, rec=0,returnmeta=True)
        if len(its)==1:
            umcori=umcori[np.newaxis,:,:,:]
        Time = np.array(its)*self.dt/86400
        dims = ('Time', 'z', 'lat', 'lon')
        coords = dict(Time=Time,
                      z=self.zc,
                      lat=self.yc[:, 0],
                      lon=self.xc[0, :])
        umcori = xr.DataArray(umcori, dims=dims, coords=coords)
        umadvec = mu.rdmds(self.path + runpath + "momDiag", iteration, rec=1)
        if len(its)==1:
            umadvec=umadvec[np.newaxis,:,:,:]
        umadvec = xr.DataArray(umadvec, dims=dims, coords=coords)
        umadvre = mu.rdmds(self.path + runpath + "momDiag", iteration, rec=2)
        if len(its)==1:
            umadvre=umadvre[np.newaxis,:,:,:]
        umadvre = xr.DataArray(umadvre, dims=dims, coords=coords)
        umdiss = mu.rdmds(self.path + runpath + "momDiag", iteration, rec=3)
        if len(its)==1:
            umdiss=umdiss[np.newaxis,:,:,:]
        umdiss = xr.DataArray(umdiss, dims=dims, coords=coords)
        umimpld = mu.rdmds(self.path + runpath + "momDiag", iteration, rec=4)
        if len(its)==1:
            umimpld=umimpld[np.newaxis,:,:,:]
        umimpld = xr.DataArray(umimpld, dims=dims, coords=coords)
        umdphix = mu.rdmds(self.path + runpath + "momDiag", iteration, rec=5)
        if len(its)==1:
            umdphix=umdphix[np.newaxis,:,:,:]
        umdphix = xr.DataArray(umdphix, dims=dims, coords=coords)
        return xr.Dataset(
            dict(umcori=umcori,
                 umadvec=umadvec,
                 umadvre=umadvre,
                 umdiss=umdiss,
                 umimpld=umimpld,
                 umdphix=umdphix))

    def midlatitude_scales(self, iteration=np.inf, lats=35, latn=45, z=10):
        it,runpath=self.get_runpath(iteration)
        u = mu.rdmds(self.path +runpath+ "U", it)
        u = xr.DataArray(
            u,
            dims=('z', 'lat', 'lon'),
            coords=dict(z=self.zc, lat=self.yc[:, 0], lon=self.xc[0, :]))
        u = u.sel(y=slice(lats, latn)).isel(z=z).mean('lat')
        n = len(u)
        ps = np.fft.rfft(u.values)
        # sample_rate = np.radians(1) * self.radius / 1000 * np.cos(
        #     np.radians((lats + latn) / 2))
        wl = np.fft.rfftfreq(n, 1)
        return (1 / wl, ps)

def cal_density(expt,T,S):
    density=expt.rhoref+expt.beta*(S-expt.Sref)-expt.alpha*(T-expt.Tref)
    density.name='Rho'
    density.attrs['units']='kg/m3'
    density.attrs['showname']='Rho (kg/m3)'
    density.attrs['longshowname']='Density (kg/m3)'
    return density

def cal_psi(expt,V):
    if not 'lon' in V.dims:
        vdr=V*expt.hfacc.mean(dim='lon')*expt.drf
    else:
        vdr=V*expt.hfacc*expt.drf
    Psi=-(vdr.reindex(z=vdr.z[::-1])).cumsum(dim='z').reindex(z=vdr.z)
    Psi.name='Psi'
    Psi.attrs['units']='m/s'
    Psi.attrs['showname']='Psi (m/s)'
    Psi.attrs['longshowname']='Meridional flow speed (m/s)'
    return Psi

def dimcomp(var,key=None,dimmeths={}):
    allvardim=var.dims
    # deal with dimensions
    for dim in allvardim:
        if not dim in dimmeths:
            continue
        dimmeth=dimmeths[dim]
        if dimmeth==None:
            continue
        elif dimmeth=='mean':
            var=var.mean(dim=dim,keep_attrs=True)
        elif dimmeth=='anom':
            var=var-var.mean(dim=dim,keep_attrs=True)
        elif type(dimmeth) is not string:
            if type(dimmeth) is float:
                if dimmeth is np.inf:
                    var=var.where((xr.DataArray(range(len(var[dim])),dims=dim)==len(var[dim])-1),drop=True).mean(dim=dim,keep_attrs=True)
                else:
                    var=var.where(abs(var[dim]-dimmeth)==min(abs(var[dim]-dimmeth)),drop=True).mean(dim=dim,keep_attrs=True)
            elif type(dimmeth) is int:
                var=var.where((xr.DataArray(range(len(var[dim])),dims=dim)==dimmeth),drop=True).mean(dim=dim,keep_attrs=True)
            elif type(dimmeth) is tuple or type(dimmeth) is list:
                if type(dimmeth[0]) is int:
                    var=var.where((xr.DataArray(range(len(var[dim])),dims=dim)>=dimmeth[0])&(xr.DataArray(range(len(var[dim])),dims=dim)<=dimmeth[1]),drop=True)
                elif type(dimmeth[0]) is float:
                    var=var.where((var[dim]>=dimmeth[0])&(var[dim]<=dimmeth[1]),drop=True)
                else:
                    raise NameError('dimmeths format wrong.')
                if type(dimmeth) is list:
                    var=var.mean(dim=dim,keep_attrs=True)
            else:
                raise NameError('dimmeths format wrong.')
        else:
            raise NameError('dimmeths format wrong.')
    return var

def my2dplt(expt,d,pltarray,wdirs=np.zeros(30),pltvar=None,dimmeths={'Time':'mean','lon':'mean','lat':None,'z':None},figsize=(12,12),projection=None,pltcontour={},flip=True,pltdir='F',labelorientation=None,xlabelpad=0,ylabelpad=20,labelaspect=20,labelshrink=1,labelpad=0.04,linetype=['k-']*20,sharex=True,sharey=True,xlims=None,ylims=None,vminmax={},cms={},savefig=0,alwaysrefmean=False,sharec=False,cutcm=[0,1],symmetrify=False,xfactor=1,interp=False):
    nplt=pltarray[0]*pltarray[1]
    dout={}
    if pltvar==None:
        pltvar=d.keys()
    if labelorientation==None:
        if pltarray[0]==1 and pltarray[1]!=1:
            labelorientation='horizontal'
        else:
            labelorientation='vertical'

    # create figure
    with plt.style.context(('labelsize15')):
        if projection=="sphere":
            fig, ax = plt.subplots(pltarray[0],pltarray[1], sharex=sharex, sharey=sharey, figsize=figsize, subplot_kw=dict(projection=ccrs.Orthographic(central_latitude=0)))
        else:
            fig, ax = plt.subplots(pltarray[0],pltarray[1], sharex=sharex, sharey=sharey, figsize=figsize)
        imh=[]

        if nplt>1:
            axf=ax.flatten(pltdir)
        else:
            axf=[ax]
        iplt=0
        nds=[]
        # loop over all vars
        for key in pltvar:
            if wdirs[iplt]:
                cm='RdBu_r'
            else:
                cm='hot' 
            if key=="S" or key=='salinity' or key=='q' or 'SHICE' in key:
                cm='viridis'
            if 'TRC' in key:
                tracercm=truncate_colormap(cmapIn='ocean_r',minval=0.05)
                cm=tracercm
            if key in cms:
                cm=cms[key]
                if cm=='tracercm':
                    tracercm=truncate_colormap(cmapIn='ocean_r',minval=0.05)
                    cm=tracercm

            #var=d[key]
            if cutcm != [0,1]:
                cm=truncate_colormap(cmapIn=cm, minval=cutcm[0], maxval=cutcm[1])
            var=dimcomp(d[key],key,dimmeths)

            newvardim=var.squeeze(drop=True).dims  
            if symmetrify and 'lat' in newvardim:
                if key=='V' or key =='Psi':
                    var.values=(var.values-var.reindex(lat=var.lat[::-1]).values)/2
                else:
                    var.values=(var.values+var.reindex(lat=var.lat[::-1]).values)/2
            labels=list(newvardim)
            if len(newvardim)!=2 and len(newvardim)!=1:
                raise NameError('{} dimension doesn''t equal to 1 or 2'.format(key))

            if len(newvardim)==1:

                dout[key]=var
                xcoord=var[newvardim[0]]
                if newvardim[0]=='lat':
                    var[0]=np.nan
                    var[-1]=np.nan
                if newvardim[0]=='z':
                    xcoord=xcoord/1e3
                    xcoord=xcoord+expt.Hice0/1e3*(not expt.topo)
                    if xlims==None and (not expt.topo):
                        xlims=[(-expt.Htot+expt.Hice0*(not expt.topo))/1e3,(expt.Hice0*(not expt.topo)-expt.Hice0*(not expt.topo))/1e3]
                if newvardim[0]=='Time':
                    xcoord=xcoord/360
                    if xlims==None:
                        xlims=[var.Time.min()/360,var.Time.max()/360]

                if not expt.spheric:
                    if newvardim[0]=='lat':
                        xcoord=xcoord/1e3
                        labels[0]="y (km)"
                        if xlims==None:
                            xlims=[np.min(xcoord),np.max(xcoord)]
                    if newvardim[0]=='lon':
                        xcoord=xcoord/1e3
                        labels[0]='x (km)'
                        if xlims==None:
                            xlims=[np.min(xcoord),np.max(xcoord)]

                if flip:
                    im=axf[iplt].plot(xcoord,var,linetype[iplt])
                    if wdirs[iplt]:
                        axf[iplt].plot(var[newvardim[0]],var*0.,'k--')
                    if hasattr(var,'showname'):
                        title=var.showname
                    else:
                        title=key
                    axf[iplt].set_title(title,loc='left')
                    plt.xlabel(labels[0])
                else:
                    im=axf[iplt].plot(var,xcoord,linetype[iplt])
                    if wdirs[iplt]:
                        axf[iplt].plot(var*0.,var[newvardim[0]],'k--')
                    if hasattr(var,'showname'):
                        title=var.showname
                    else:
                        title=key
                    axf[iplt].set_title(title,loc='left')
                    plt.ylabel(labels[0])

                if not ylims==None:
                    axf[iplt].set_ylim(ylims)
                if not xlims==None:
                    axf[iplt].set_xlim(xlims)

                nds=nds+[1]
                iplt=iplt+1

            else: # dim==2
                # flip direction if required
                if flip:
                    xdim=newvardim[1]
                    ydim=newvardim[0]
                else:
                    xdim=newvardim[0]
                    ydim=newvardim[1]
                    var=var.transpose()

                dout[key]=var

                xcoord=var[xdim]
                ycoord=var[ydim]
                if xdim=='z':
                    xcoord=xcoord/1e3
                    xcoord=xcoord+expt.Hice0/1e3*(not expt.topo)
                    if xlims==None and (not expt.topo):
                        xlims=[(-expt.Htot+expt.Hice0*(not expt.topo))/1e3,(expt.Hice0*(not expt.topo)-expt.Hice0*(not expt.topo))/1e3]
                if ydim=='z':
                    ycoord=ycoord/1e3
                    ycoord=ycoord+expt.Hice0/1e3*(not expt.topo)
                    if ylims==None and (not expt.topo):
                        ylims=[(-expt.Htot+expt.Hice0*(not expt.topo))/1e3,(expt.Hice0*(not expt.topo)-expt.Hice0*(not expt.topo))/1e3]
                if xdim=='Time':
                    xcoord=xcoord/360
                    if xlims==None:
                        xlims=[var.Time.min()/360,var.Time.max()/360]
                if ydim=='Time':
                    ycoord=ycoord/360
                    if ylims==None:
                        ylims=[var.Time.min()/360,var.Time.max()/360]
                if not xlims==None:
                    xcoord=xcoord.where((xcoord>=xlims[0])&(xcoord<=xlims[1]),drop=True)
                    var=var.where((xcoord>=xlims[0])&(xcoord<=xlims[1]),drop=True)
                if not ylims==None:
                    ycoord=ycoord.where((ycoord>=ylims[0])&(ycoord<=ylims[1]),drop=True)
                    var=var.where((ycoord>=ylims[0])&(ycoord<=ylims[1]),drop=True)

                if not expt.spheric:
                    if xdim=='lat':
                        xcoord=xcoord/1e3
                    if ydim=='lat':
                        ycoord=ycoord/1e3
                    if xdim=='lon':
                        xcoord=xcoord/1e3
                    if ydim=='lon':
                        ycoord=ycoord/1e3

                # offset
                offset=0
                if wdirs[iplt]==0:
                    if key=='S' and not alwaysrefmean:
                        offset=expt.Sref
                    elif key=='T' and not alwaysrefmean:
                        offset=expt.Tref
                    else:
                        varmean=var.mean(dim=[xdim,ydim]).values
                        print('{} mean:{}'.format(key,varmean))
                        if np.isnan(varmean):
                            print(var)
                        if (not np.isnan(varmean)) and varmean!=0:
                            offset=np.round(varmean,int(5-np.log10(np.abs(varmean))))

                # vmin vmax
                if key in vminmax:
                    v_min=vminmax[key][0]
                    v_max=vminmax[key][1]
                elif wdirs[iplt]==1:
                    vamp=abs(var).max()
                    v_max=vamp
                    v_min=-vamp
                else:
                    v_max=None
                    v_min=None
                xcoord=xcoord-xcoord.mean()
                
                # interpolate
                vari=var
                if interp:
                    xorder=np.sign(xcoord[1]-xcoord[0]).values.astype(int)
                    yorder=np.sign(ycoord[1]-ycoord[0]).values.astype(int)
                    interp = RegularGridInterpolator((ycoord[::yorder], xcoord[::xorder]), var.values[::yorder,::xorder])
                    xi = np.linspace(xcoord.min(), xcoord.max(), 500)
                    yi = np.linspace(ycoord.min(), ycoord.max(), 500)
                    Xi, Yi = np.meshgrid(xi, yi)
                    vari = interp((Yi, Xi))
                    xcoord=xi
                    ycoord=yi

                # plot
                if projection=='sphere':
                    im=axf[iplt].pcolormesh((xcoord)*xfactor,ycoord, vari-offset,transform=ccrs.PlateCarree(),cmap=cm,vmin=v_min,vmax=v_max)
                    gl = axf[iplt].gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='k', alpha=0.5)
                    if key in pltcontour:
                        if wdirs[iplt]==0:
                            linestyles='solid'
                        else:
                            linestyles=None
                        axf[iplt].contour((xcoord)*xfactor,ycoord, vari-offset,pltcontour[key],linestyles=linestyles,colors='gray',transform=ccrs.PlateCarree(),linewidths=1)
                else:
                    im=axf[iplt].pcolormesh((xcoord)*xfactor,ycoord, vari-offset,cmap=cm,vmin=v_min,vmax=v_max)
                    if key in pltcontour:
                        if wdirs[iplt]==0:
                            linestyles='solid'
                        else:
                            linestyles=None
                        axf[iplt].contour((xcoord)*xfactor,ycoord, vari-offset,pltcontour[key],linestyles=linestyles,colors='gray',linewidths=2)

                if offset!=0:
                    if hasattr(var,'showname'):
                        title=var.showname.replace(')',',ref:{:g})'.format(round(offset,3)))
                    else:
                        title='{} (ref:{:g})'.format(key,round(offset,3))

                else:
                    if hasattr(var,'showname'):
                        title=var.showname
                    else:
                        title=key

                axf[iplt].set_title(title,loc='left')
                if not sharec:
                    cb = fig.colorbar(im, ax=axf[iplt],orientation=labelorientation, shrink=labelshrink,aspect=labelaspect, pad=labelpad)
                   # cb.set_label(key)
                    cb.formatter.set_powerlimits((-1,2))
                    cb.formatter.set_useOffset(False)
                    cb.update_ticks()

                if not expt.spheric:
                    axf[iplt].set_ylim([np.min(ycoord),np.max(ycoord)])
                    axf[iplt].set_xlim([np.min(xcoord),np.max(xcoord)])
                #if not ylims==None:
                #    axf[iplt].set_ylim(ylims)
                #if not xlims==None:
                #    axf[iplt].set_xlim(xlims)

                nds=nds+[2]
                iplt=iplt+1
            
            imh.append(im)
            # common axis label
            if all(np.array(nds)==2):
                xlab=xdim
                ylab=ydim
                if xdim=='z':
                    xlab=xlab+' (km)'
                if ydim=='z':
                    ylab=ylab+' (km)'
                if xdim=='Time':
                    xlab=xlab+' (yr)'
                if ydim=='Time':
                    ylab=ylab+' (yr)'
                if not expt.spheric:
                    if xdim=='lat':
                        xlab='y (km)'
                    if ydim=='lat':
                        ylab='y (km)'
                    if xdim=='lon':
                        xlab='x (km)'
                    if ydim=='lon':
                        ylab='x (km)'
                if projection==None:
                    fig.add_subplot(111, frameon=False)
                    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
                    plt.xlabel(xlab,labelpad=xlabelpad)
                    plt.ylabel(ylab,labelpad=ylabelpad)
        
            if sharec:
                plt.tight_layout()
                if labelorientation is 'vertical':
                    #fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.02, hspace=0.02)
                    cx=fig.add_axes([1.0,0.24,0.015,0.65])
                    #cx=fig.add_axes([0.5,0.1,0.6,0.8])
                    #cx.set_axis_off()
                    cb = fig.colorbar(im,cax=cx,orientation=labelorientation, shrink=labelshrink,aspect=labelaspect, pad=labelpad)
                    cb.formatter.set_powerlimits((-1,2))
                    cb.formatter.set_useOffset(False)
                    cb.update_ticks()

                if labelorientation is 'horizontal':
                    #fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.02, hspace=0.02)
                    cx=fig.add_axes([0.2,-0.02,0.6,0.015])
                    #cx=fig.add_axes([0.5,0.1,0.6,0.8])
                    #cx.set_axis_off()
                    cb = fig.colorbar(im,cax=cx,orientation=labelorientation, shrink=labelshrink,aspect=labelaspect, pad=labelpad)
                    cb.formatter.set_powerlimits((-1,2))
                    cb.formatter.set_useOffset(False)
                    cb.update_ticks()
    # save figure
    figname=None
    if savefig:
        figname='{}_clim.png'.format(expt.name)
        fig.savefig(figname, bbox_inches='tight', dpi=150)
    return fig,axf,dout,figname,imh
        



def get_dudr(expt, iteration=np.nan):
    it,runpath=expt.get_runpath(expt,iteration)
    _,its,meta=mu.rdmds(self.path + runpath + 'dynDiag', np.inf,returnmeta=True,rec=0)
    iUVEL=meta['fldlist'].index('UVEL')
    u = self.get_var(runpath+'dynDiag',iteration=iteration,rec=iUVEL)
    if np.isnan(it):
        u = np.mean(u, axis=0)
    udy = u * expt.dyg[np.newaxis, :, :] * expt.deepfacc[:, np.newaxis,
                                                         np.newaxis]
    udy = np.concatenate((-udy[[0]], udy, -udy[[-1]]), axis=0)
    du = np.diff(udy, axis=0)
    drc = mu.rdmds(self.path + 'DRC').squeeze()
    dudr = du / (drc[:, np.newaxis, np.newaxis] * expt.dyg[np.newaxis, :, :] *
                 expt.deepfacf[:, np.newaxis, np.newaxis])
    dudr = 0.5 * (dudr[:-1] + dudr[1:])
    dudr = np.mean(dudr, axis=2)
    return dudr


def get_dudr_old(expt, iteration=np.nan):
    it,runpath=expt.get_runpath(iteration)
    _,its,meta=mu.rdmds(self.path + runpath + 'dynDiag', np.inf,returnmeta=True,rec=0)
    iUVEL=meta['fldlist'].index('UVEL')
    u = self.get_var(runpath+'dynDiag',iteration=iteration,rec=iUVEL)
    u = np.mean(u, axis=-1)
    if np.isnan(it):
        u = np.mean(u, axis=0)
    u = np.concatenate((-u[[0], :], u, -u[[-1], :]), axis=0)
    du = np.diff(u, axis=0)
    drc = mu.rdmds(self.path + 'DRC').squeeze()
    dudr = du / drc[:, np.newaxis]
    dudr = 0.5 * (dudr[:-1] + dudr[1:])
    return dudr


def get_1byrdudtheta(expt, iteration=np.nan):
    it,runpath=expt.get_runpath(iteration)
    _,its,meta=mu.rdmds(self.path + runpath + 'dynDiag', np.inf,returnmeta=True,rec=0)
    iUVEL=meta['fldlist'].index('UVEL')
    u = self.get_var(runpath+'dynDiag',iteration=iteration,rec=iUVEL)
    if np.isnan(it):
        u = np.mean(u, axis=0)
    u = np.concatenate(
        (-u[:, [0]], u),
        axis=1)  # Assumption of no-slip at southern and northen edge
    dudrdtheta = np.diff(
        u, axis=1
    ) / (expt.dyc[np.newaxis, :, :] * expt.deepfacc[:, np.newaxis, np.newaxis])
    dudrdtheta = np.mean(dudrdtheta, axis=2)
    return dudrdtheta


def get_1byrdudtheta_old(expt, iteration=np.nan):
    it,runpath=expt.get_runpath(iteration)
    _,its,meta=mu.rdmds(self.path + runpath + 'dynDiag', np.inf,returnmeta=True,rec=0)
    iUVEL=meta['fldlist'].index('UVEL')
    u = self.get_var(runpath+'dynDiag',iteration=iteration,rec=iUVEL)
    if np.isnan(it):
        u = np.mean(u, axis=0)
    u = np.concatenate(
        (-u[:, [0]], u),
        axis=1)  # Assumption of no-slip at southern and northen edge
    dtheta = np.radians(np.diff(expt.yc[:, 0])[0])
    dy = (expt.rsphere + expt.zc) * dtheta
    durdtheta = np.diff(u, axis=1) / dy[:, np.newaxis, np.newaxis]
    durdtheta = np.mean(durdtheta, axis=2)
    return durdtheta


def get_dudz(expt, iteration=np.nan):
    it,runpath=expt.get_runpath(iteration)
    dudr = get_dudr_old(expt, it)
    durdtheta = get_1byrdudtheta_old(expt, it)

    theta = np.radians(expt.yc[:, 0])
    dudz = -np.sin(theta) * dudr + np.cos(theta) * durdtheta
    #dudz = np.cos(theta) * durdtheta
    return xr.DataArray(
        dudz, dims=('z', 'lat'), coords=dict(z=expt.zc, lat=expt.yc[:, 0]))


def get_db_rdtheta_old(expt, iteration=np.nan):
    it,runpath=expt.get_runpath(iteration)
    _,its,meta=mu.rdmds(self.path + runpath + 'dynDiag', np.inf,returnmeta=True,rec=0)
    iTHETA=meta['fldlist'].index('THETA')
    T = self.get_var(runpath+'dynDiag',iteration=iteration,rec=iTHETA).mean('time').mean('lon')
    if np.isnan(it):
        T = np.mean(T, axis=0)
    b = expt.g * expt.alpha * T
    b = np.concatenate((b[:, [-1]], b), axis=1)
    dtheta = np.radians(np.diff(expt.yc[:, 0])[0])
    dy = (expt.rsphere + expt.zc) * dtheta
    dbdy = np.diff(b, axis=1) / dy[:, np.newaxis, np.newaxis]
    dbrdtheta = dbdy / 2 / expt.omega
    dbrdtheta = np.mean(dbrdtheta, axis=2)
    return xr.DataArray(
        dbrdtheta, dims=('z', 'lat'), coords=dict(z=expt.zc, lat=expt.yc[:, 0]))



def get_db_rdtheta(expt, iteration=np.nan):
    it,runpath=expt.get_runpath(iteration)
    _,its,meta=mu.rdmds(self.path + runpath + 'dynDiag', np.inf,returnmeta=True,rec=0)
    iTHETA=meta['fldlist'].index('THETA')
    T = mu.rdmds(self.path + 'dynDiag', it, rec=5)
    if np.isnan(it):
        T = np.mean(T, axis=0)
    b = expt.g * expt.alpha * T
    bdx = b * expt.dxc * expt.deepfacc[:, np.newaxis, np.newaxis]
    bdx = np.concatenate((bdx[:, [-1]], bdx), axis=1)
    dbrdtheta = np.diff(
        bdx, axis=1) / (expt.ras[np.newaxis, :, :] * expt.
                        deepfacc[:, np.newaxis, np.newaxis]**2)
    dbrdtheta /= (2 * expt.omega)
    dbrdtheta = np.mean(dbrdtheta, axis=2)
    return xr.DataArray(
        dbrdtheta, dims=('z', 'lat'), coords=dict(z=expt.zc, lat=expt.yc[:, 0]))


def get_thermal_wind(expt, iteration=np.nan):
    it,runpath=expt.get_runpath(iteration)
    dudz = get_dudz(expt, it)
    dbrdtheta = get_db_rdtheta_old(expt, it)
    return dudz, dbrdtheta


def integrate_dudz_in_z(expt):
    dudz, rhs = get_thermal_wind(expt)

    dr = np.diff(expt.zf)[:, np.newaxis]
    phi = np.radians(expt.yc[np.newaxis, :, 0])
    sinphirhsdr = rhs * np.sin(phi) * dr
    intsinphirhsdr = np.cumsum(sinphirhsdr, axis=0)

    Rdphi = np.mean(
        expt.dyc[np.newaxis, :, :] * expt.deepfacc[:, np.newaxis, np.newaxis],
        axis=2)
    Rcosphidphi = Rdphi * np.cos(phi)
    rhsRcosphidphi = rhs * Rcosphidphi
    intrhsRcosphidphi = np.cumsum(rhsRcosphidphi, axis=1)
    return intsinphirhsdr + intrhsRcosphidphi
    return Rdphi


def remap_from_rphi_toxz(var, expt, npoints=1000, method='linear'):
    phi = np.radians(expt.yc[:, 0])[np.newaxis, :]
    r = expt.rsphere + expt.zc[:, np.newaxis]
    #phi = np.radians(expt.yc[:, 0])
    #r = expt.rsphere + expt.zc
    xp, zp = r * np.cos(phi), r * np.sin(phi)
    xi, zi = np.linspace(0, expt.rsphere, npoints), np.linspace(
        -expt.rsphere, expt.rsphere, 2 * npoints)
    x, z = np.meshgrid(xi, zi)
    var_new = griddata(
        (zp.ravel(), xp.ravel()), var.ravel(), (z, x), method=method)
    var_new[x**2 + z**2 < r[-1]**2] = np.inf
    return xr.DataArray(var_new, dims=('z', 'x'), coords=dict(z=zi, x=xi))


def integrate_interpolated_dudz_in_z(expt):
    it,runpath=expt.get_runpath(np.nan)
    u = np.mean(
        mu.rdmds(self.path +runpath+ "dynDiag", np.nan, rec=0).squeeze(), axis=(0, -1))
    uinterp = remap_from_rphi_toxz(u, expt)

    rhs = get_db_rdtheta(expt)
    rhs_interpolated = remap_from_rphi_toxz(rhs.values, expt, npoints=1000)

    ro = expt.rsphere
    ri = expt.zc[-1] + ro

    x = rhs_interpolated['x'].values
    z = rhs_interpolated['z'].values
    dz = np.diff(z)[0]
    xx, zz = np.meshgrid(x, z)

    # phi = np.arctan(zz/xx)

    # region1 = (phi < - np.acos(r[-1]/expt.rsphere)) & (zz<0)
    # region3 = (phi >   np.acos(r[-1]/expt.rsphere)) & (zz>0)
    # region2 = (phi <= np.acos(r[-1]/expt.rsphere)) & (phi >= -np.acos(r[-1]/expt.rsphere))

    ufromth = np.nancumsum(rhs_interpolated, axis=0)
    # ufromth[xx**2 + zz**2 < ri**2] = np.inf
    # ufromth[xx**2 + zz**2 > ro**2] = np.inf
    u_fromth = xr.DataArray(
        ufromth * dz, dims=('z', 'x'), coords=dict(z=z, x=x))
    u_interp = xr.DataArray(uinterp, dims=('z', 'x'), coords=dict(z=z, x=x))
    cI = (u_fromth - u_interp).sel(
        z=0, x=(ro + ri) / 2, method='nearest').values
    #u_fromth += cI
    return u_fromth, u_interp


def integrate_dbdr(expt, iteration=np.nan):
    rhs = get_dudz(expt, iteration)
    rhs *= (2 * expt.omega)
    r = (expt.rsphere + expt.zc)
    theta = expt.yc[:, 0]
    dtheta = np.radians(np.diff(theta)[0])
    rhs *= r[:, np.newaxis]
    b = np.zeros_like(rhs)
    b[:, theta > 0] = np.cumsum(rhs[:, theta > 0], axis=1) * dtheta
    b[:, theta < 0] = -np.cumsum(
        rhs[:, theta < 0][:, ::-1], axis=1)[:, ::-1] * dtheta
    #b = np.cumsum(rhs, axis=1) * dtheta
    return b





# #suyash_enc1=Experiment(base_dir="/net/fs08/d0/bire/mitgcm/enceladus_sector/run_1/")
# #suyash_enc2=Experiment(base_dir="/net/fs08/d0/bire/mitgcm/enceladus_sector/run_2/")
# v12_lowgamma0=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_lowgamma0")
# v12_lowgamma0_weakic=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_lowgamma0_weakic")
# v12_lowgamma0_lowdiff=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_lowgamma0_lowdiff")
# #v12_lowgamma0_allshell=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_lowgamma0_allshell")
# v12_lowgamma0_flux10=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_lowgamma0_flux10")
# v12_lowgamma0_flux100=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_lowgamma0_flux100")
# v12_lowgamma0_flux100_weakic=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_lowgamma0_flux100_weakic")
# #v12_lowgamma0_lowvert=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_lowgamma0_lowvert")
# v12_lowgamma0_lowvert_weakic=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_lowgamma0_lowvert_weakic")
# v12_lowgamma0_lowvert_weakic_adv33=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_lowgamma0_lowvert_weakic_adv33")
# #v12_lowgamma0_lowlowvert_weakic=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_lowgamma0_lowlowvert_weakic")
# v12_lowgamma0_lowlowvert_weakweakic=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_lowgamma0_lowlowvert_weakweakic")
# v12_lowgamma0_lowvert_weakic_allshell=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_lowgamma0_lowvert_weakic_allshell")
# v12_lowgamma0_lowvert_weakic_halfhalf=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_lowgamma0_lowvert_weakic_halfhalf")
# #v12_topo_lowgamma0_lowvert_weakic=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_topo_lowgamma0_lowvert_weakic")
# #v12_topo_lowgamma0_lowvert_weakic_allshell=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_topo_lowgamma0_lowvert_weakic_allshell")
# #v12_topo_lowgamma0_lowvert_weakic_halfhalf=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_topo_lowgamma0_lowvert_weakic_halfhalf")
# v12_lowgamma0_lowvert_weakic_4psu=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_lowgamma0_lowvert_weakic_4psu")
# v12_lowgamma0_lowvert_weakic_allshell_4psu=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_lowgamma0_lowvert_weakic_allshell_4psu")
# v12_lowgamma0_lowvert_weakic_halfhalf_4psu=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_lowgamma0_lowvert_weakic_halfhalf_4psu")
# #v12_topo_lowgamma0_lowvert_weakic_4psu=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="enceladus_z40x90y324_v12_topo_lowgamma0_lowvert_weakic_4psu")
# #v12_topo_lowgamma0_lowvert_weakic_allshell_4psu=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="topo_lowgamma0_lowvert_weakic_allshell_4psu")
# #v12_topo_lowgamma0_lowvert_weakic_halfhalf_4psu=Experiment(base_dir="/net/fs09/d0/wanying/data_",name="topo_lowgamma0_lowvert_weakic_halfhalf_4psu")
#  
#  
# expts_lowres = [v12_lowgamma0,v12_lowgamma0_weakic,v12_lowgamma0_lowdiff,v12_lowgamma0_flux10,v12_lowgamma0_flux100,v12_lowgamma0_flux100_weakic,v12_lowgamma0_lowvert_weakic,v12_lowgamma0_lowvert_weakic_adv33,v12_lowgamma0_lowlowvert_weakweakic,v12_lowgamma0_lowvert_weakic_allshell,v12_lowgamma0_lowvert_weakic_halfhalf]
# 
# v13_40psu_allcore_strongic=Experiment(base_dir="/net/fs08/d0/wanying/data_",name="enceladus_z59x90y648_v13_40psu_allcore_strongic")
# v13_40psu_halfhalf=Experiment(base_dir="/net/fs08/d0/wanying/data_",name="enceladus_z59x90y648_v13_40psu_halfhalf")
# v13_40psu_allshell=Experiment(base_dir="/net/fs08/d0/wanying/data_",name="enceladus_z59x90y648_v13_40psu_allshell")
# v13_40psu_allcore_strongic_merid=Experiment(base_dir="/net/fs08/d0/wanying/data_",name="enceladus_z59x90y648_v13_40psu_allcore_strongic_merid")
# v13_40psu_allcore_strongic_merid_asym_1e_3vr=Experiment(base_dir="/net/fs08/d0/wanying/data_",name="enceladus_z59x90y648_v13_40psu_allcore_strongic_merid_asym_1e-3vr")
# v13_40psu_allcore_strongic_merid_asym_1e_5vr=Experiment(base_dir="/net/fs08/d0/wanying/data_",name="enceladus_z59x90y648_v13_40psu_allcore_strongic_merid_asym_1e-5vr")
# v13_40psu_allcore_strongic_merid_asym_1e_5vr_1e_2vh=Experiment(base_dir="/net/fs08/d0/wanying/data_",name="enceladus_z59x90y648_v13_40psu_allcore_strongic_merid_asym_1e-5vr_1e-2vh")
# v13_40psu_halfhalf_strongic_merid=Experiment(base_dir="/net/fs08/d0/wanying/data_",name="enceladus_z59x90y648_v13_40psu_halfhalf_strongic_merid")
# v13_40psu_allshell_strongic_merid=Experiment(base_dir="/net/fs08/d0/wanying/data_",name="enceladus_z59x90y648_v13_40psu_allshell_strongic_merid")
# 
# expts_highres=[v13_40psu_allcore_strongic,v13_40psu_halfhalf,v13_40psu_allshell, v13_40psu_allcore_strongic_merid, v13_40psu_allcore_strongic_merid_asym_1e_3vr,v13_40psu_allcore_strongic_merid_asym_1e_5vr,v13_40psu_allcore_strongic_merid_asym_1e_5vr_1e_2vh,v13_40psu_halfhalf_strongic_merid,v13_40psu_allshell_strongic_merid]

def expt_grid_meanvelocity(
        expts,
        iteration=np.inf,
        runnum=np.inf,
        runpath=None,
        sample=3,
        # ncols=4,
        lats=45,
        latn=50,
        func=np.mean,
        sharex=False,
        figsize=(9, 4)):
    labs = iter(string.ascii_lowercase)
    nax = len(expts)
    # fig, ax = plt.subplots(int(nax / ncols), ncols, figsize=figsize)
    fig, ax = plt.subplots(1, 4, figsize=figsize, sharex=sharex, sharey=True)
    meanvels = []
    ucones = []
    urots = []
    cols = []
    #styledict = {29000: '--', 100000: '-.', 500000: '-'}
    #coldict = {0.1: 'r', 10: 'b', 1000: 'k', 10000: 'g'}
    for expt in expts:
        u_prime = expt.get_uprime(lats, latn, iteration=iteration,runnum=runnum,runpath=runpath,sample=sample)

        z = u_prime.z / expt.Htot
        ls = expt.lstyle  #styledict[expt.htot]
        cs = expt.c  #coldict[expt.q]
        ax[0].plot(
            u_prime.values, u_prime.z / expt.Htot, ls=ls, c=cs, label=expt.name)
        ax[1].plot(
            u_prime.values / expt.ucone,
            u_prime.z / expt.Htot,
            ls=ls,
            c=cs,
            label=expt.name)
        ax[2].plot(
            u_prime.values / expt.unonrot,
            u_prime.z / expt.Htot,
            ls=ls,
            c=cs,
            label=expt.name)
        ax[3].plot(
            u_prime.values / expt.urot,
            u_prime.z / expt.Htot,
            ls=ls,
            c=cs,
            label=expt.name)
        meanvels.append(func(u_prime))
        urots.append(expt.urot)
        ucones.append(expt.ucone)
        #cols.append(coldict[expt.Htot])
        cols.append(expt.c)
    return fig, meanvels, urots, ucones, cols

def truncate_colormap(cmapIn='jet', minval=0.0, maxval=1.0, n=100):
    '''truncate_colormap(cmapIn='jet', minval=0.0, maxval=1.0, n=100)'''
    cmapIn = plt.get_cmap(cmapIn)

    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmapIn.name, a=minval, b=maxval),
        cmapIn(np.linspace(minval, maxval, n)))

    return new_cmap

def expt_grid_rmst(
        expts,
        iteration=np.inf,
        runnum=np.inf,
        runpath=None,
        sample=3,
        # ncols=4,
        lats=45,
        latn=50,
        func=np.mean,
        figsize=(9, 4)):
    labs = iter(string.ascii_lowercase)
    nax = len(expts)
    # fig, ax = plt.subplots(int(nax / ncols), ncols, figsize=figsize)
    fig, ax = plt.subplots(1, 3, figsize=figsize, sharey=True)
    meanTs = []
    gprimes = []
    cols = []
    #styledict = {29000: '--', 100000: '-.', 500000: '-'}
    #coldict = {0.1: 'r', 10: 'b', 1000: 'k', 10000: 'g'}
    for expt in expts:
        Trms = expt.get_tprime(lats, latn, iteration=iteration,runnum=runnum,runpath=runpath,sample=sample)
        ls = expt.lstyle  #styledict[expt.Htot]
        cs = expt.c  #coldict[expt.Q0]
        z = Trms.z / expt.Htot
        zlims = (z > -0.9) & (z < -0.1)
        ax[0].plot(
            Trms.values[zlims],
            Trms.z[zlims] / expt.Htot,
            ls=ls,
            c=cs,
            label=expt.name)
        #ax[1].plot(urms.values / expt.urot, urms.z / expt.Htot,ls=ls,c=cs, label=expt.name)
        ax[1].plot(
            Trms.values[zlims] / expt.gprime,
            Trms.z[zlims] / expt.Htot,
            ls=ls,
            c=cs,
            label=expt.name)
        ax[2].plot(
            Trms.values[zlims] / expt.gprimenonrot,
            Trms.z[zlims] / expt.Htot,
            ls=ls,
            c=cs,
            label=expt.name)
        meanTs.append(func(Trms))
        gprimes.append(expt.gprime)
        #cols.append(coldict[expt.Htot])
        cols.append(expt.c)
    return fig, meanTs, gprimes, cols


def expt_grid_powerspectrums(expts,
                             iteration=np.inf,
                             lats=35,
                             latn=45,
                             ncols=4,
                             z=10,
                             figsize=(12, 9)):
    labs = iter(string.ascii_lowercase)
    nax = len(expts)
    fig, ax = plt.subplots(int(np.ceil(nax / ncols)), ncols, figsize=figsize)
    for axc, expt in zip(ax.ravel(), expts):
        wl, ps = expt.midlatitude_scales(iteration, lats=lats, latn=latn, z=z)
        axc.semilogx(wl, np.abs(ps)**2, 'k-')
        axc.grid()
        axc.set_title('(' + next(labs) + r') ' + expt.name)
        # axc.text(
        #     0.1,
        #     -0.1,
        #     'RoC = {:1.2e}'.format(expt.rostar),
        #     transform=axc.transAxes,
        #     fontsize=12,
        #     family='monospace')
        # axc.text(
        #     0.1,
        #     -0.2,
        #     'Ek = {:1.2e}'.format(expt.E),
        #     transform=axc.transAxes,
        #     fontsize=12,
        #     family='monospace')
    return fig


def save_netcdf(expt,dic,fname):
    ds=xr.Dataset(dic)
    ds.to_netcdf(expt.base_dir + fname)
