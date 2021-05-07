# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 23:46:04 2019

@author: Karen
"""

toi       = 0.01        #toi number (number portion only - i.e. no "TOI"), set to -1 if no TOI number is available
deltamag  = 0           #delta magnitude of faintest neighbor cleared, or delta magnitude of NEB, set to 0 to leave blank
transcov  = 'Full'      #Full Ingress Egress Out of Transit (CASE SENSITIVE!!!)
notes     = ''          #public note such as "deep" etc. - do not put proprietary results here

camera    = 'SINISTRO'  #Camera name (  Spectral  SINISTRO  SBIG STX6303  SBIG STX16803  )
telsize   = 1.0         #telescope aperture in meters (  2.0  1.0  0.4  0.6  )
pixscale  = 0.389       #pixel scale in asec  (  2m0 = 0.304  1m0 = 0.389  0m4 = 0.57  ULMT = 0.39  )

skipsummaryupload = 0   #set to 1 to skip uploading observation summary, set to 0 to upload observation summary
skipfileupload = 0      #set to 1 to skip file uploads, set to 0 to upload matching files

#-----------------------------------------------------------------------------

username = 'xxxxxxxxxx'
password = 'xxxxxxxxxx'

#-----------------------------------------------------------------------------

version = '2.02'
versiondate = '2020.06.06'
print('ExoFOP-TESS Upload version '+version+' released '+versiondate)
print('')

import os ,re, requests, sys, tkinter
from astropy.table import Table
import numpy as np

homedir = os.path.expanduser('~')       # find the home directory
path = os.getcwd()       # get current directory
fileList = os.listdir(path)
for fileName in fileList:
    if os.path.isfile(os.path.join(path,fileName)) and fileName.endswith('.tbl'):
        t = Table.read(fileName, format='ascii.tab', data_start=1)
        break

fwhm = pixscale*np.median(t['FWHM_T1'])  #median FWHM in asec
psf = str(round(fwhm,2))
print("FWHM in arcsec: "+psf)

obsstart = t['JD_UTC'][0]       #start of observations (fractional JD)
obsend = t['JD_UTC'][-1]      #end of observations (fractional JD)
obsdur = str(round(abs((obsend-obsstart)*24*60)))
print("Obs length in minutes: "+obsdur)

obsnum = str(t['JD_UTC'].size)     #number of exposures
print("Number of exposures: "+obsnum)

photaprad = str(round(np.mean(t['Source_Radius']),1)) #mean aperture radius in pixels
print("Mean aperture radius in pixels: "+photaprad)

apradarcsec = str(round(pixscale*np.mean(t['Source_Radius']),1))#mean aperture radius in arcsec
print("Mean aperture radius in arcsec: "+apradarcsec)

pixscale = str(pixscale)
telsize = str(telsize)

if deltamag == 0:
    deltamag=''
else:
    deltamag = str(deltamag)
if toi==0:
    tkinter.Tk().bell()
    tkinter.Tk().bell()
    input("Enter TOI number in file")
    exit()
elif toi>0:
    toi='TOI'+str(toi)
    planet=toi.split('.')[1]
else:
    toi=''
    planet='01'

if transcov.lower() =='full':
    transcov = 'Full'
elif transcov.lower() =='ingress':
    transcov = 'Ingress'
elif transcov.lower() =='egress':
    transcov = 'Egress'
else:
    transcov = 'Out of Transit'

pieces=""

for fileName in fileList:
    if os.path.isfile(os.path.join(path,fileName)) and fileName.startswith('TIC') and not fileName.startswith('TIC '):
        pieces=fileName.split("_")
        break
print(pieces)
if len(pieces) < 4:
    input("Filenames must have at least three underscores...")
    tkinter.Tk().bell()
    tkinter.Tk().bell()
    exit()
tic=pieces[0].split('-')[0]
tic=re.search(r'([0-9]+)',tic).group(1)
date=pieces[1]
print(date[0:4]+'-'+date[4:6]+'-'+date[6:8])
observatory = pieces[2]
filterband = pieces[3].split('.')[0]
tag=date+'_'+username+'_TIC'+tic+'_'+planet
print(tic)
emailtitle='TIC '+tic+'.'+planet+' (TOI-'+toi[3:]+') on UT'+date[0:4]+'.'+date[4:6]+'.'+date[6:8]+' from '+observatory+' in '+filterband
f= open(emailtitle+".txt","w+")
f.write(emailtitle+"\r\n")
f.close()
entries = {
        'planet': toi,
        'tel': observatory,
        'telsize': telsize,
        'camera': camera,
        'filter': filterband,
        'pixscale': pixscale,
        'psf': psf,
        'photaprad': photaprad,
        'obsdate': date[0:4]+'-'+date[4:6]+'-'+date[6:8],
        'obsdur': obsdur,
        'obsnum': obsnum,
        'obstype': 'Continuous',
        'transcov': transcov,
        'deltamag': deltamag,
        'tag': tag,
        'groupname': 'tfopwg',
        'notes': notes,
        'id': tic
    }

print(entries)


if os.path.exists(path):
    credentials = {
    'username': username,
    'password': password,
    'ref': 'login_user',
    'ref_page': '/tess/'
    }
    with requests.Session() as session:
        response1 = session.post('https://exofop.ipac.caltech.edu/tess/password_check.php', data=credentials)
        if response1 and not re.search('Invalid username or password', response1.text):
            print('\nLogin OK.')
        else:
            tkinter.Tk().bell()
            tkinter.Tk().bell()
            sys.exit('\nERROR:  Login did not work; check username or password.')
            
        if not skipsummaryupload:
            response2 = session.post('https://exofop.ipac.caltech.edu/tess/insert_tseries.php', data=entries)
            if response2:
                print('\nAdded new Time Series...')
            else:
                tkinter.Tk().bell()
                tkinter.Tk().bell()
                sys.exit('\nERROR: Time Series Add failed.')
        else:
            print('Skipped observation summary upload per user request.')
        if not skipfileupload:
            fileList = os.listdir(path)
            for fileName in fileList:
                if os.path.isfile(os.path.join(path,fileName)) and fileName.startswith('TIC') and not fileName.startswith('TIC '):
                    pieces=fileName.split("_")
                    tic=pieces[0].split("-")[0]
                    tic=re.search(r'([0-9]+)',tic).group(1)
                    date=pieces[1]
                    description=''
                    if fileName.find('bjd-flux-err')>-1:
                        description='Photometry Table Subset for Joint Fitting'
                    elif fileName.endswith('field.png'):
                        description='Field Image with Apertures'
                    elif fileName.endswith('field-zoom.png'):
                        description='Zoomed Field Image with Apertures'  
                    elif fileName.endswith('field-zoom2.png'):
                        description='Zoomed Field Image with Apertures'
                    elif fileName.find('lightcurve')>-1 and fileName.endswith('.png'):
                        description='Light Curve Plot'
                    elif fileName.endswith('.apertures'):
                        description='AstroImageJ Photometry Aperture File'
                    elif fileName.endswith('.plotcfg'):
                        description='AstroImageJ Plot Configuration File'
                    elif fileName.endswith('.tbl'):
                        description='AstroImageJ Full Photometry Measurements Table'
                    elif fileName.endswith('seeing-profile.png'):
                        description='Seeing Profile'
                    elif fileName.endswith('WCS.fits'):
                        description='Plate Solved Representative FITS Image'  
                    elif fileName.endswith('WCS.fits.fz'):
                        description='Plate Solved Representative FITS Image'
                    elif fileName.endswith('notes.txt'):
                        description='Results and Observing Notes'
                    elif fileName.endswith('NEBcheck.zip'):
                        description='Light Curve Plots of Nearby Field Stars'    
                    elif fileName.endswith('NEB-table.txt'):
                        description='Summary of NEBcheck Results'    
                    elif fileName.endswith('dmagRMS-plot.png'):
                        description='Dmag vs. RMS plot'
                    elif fileName.find('merged')>-1 and fileName.endswith('.png'):
                        description='Merged Light Curve Plot'
                    if description == '':
                        print('******NOT UPLOADED: '+fileName)
                    else:
                        #print(tic, toi, date, tag, description)
                        files = {'file_name': open(fileName, 'rb')}
                        payload = {
                            'file_type': 'Light_Curve',
                            'planet': toi,
                            'file_desc': description,
                            'file_tag': tag,
                            'groupname': 'tfopwg',
                            'propflag': 'on',
                            'id': tic
                            }
                        response3 = session.post('https://exofop.ipac.caltech.edu/tess/insert_file.php', files=files, data=payload)
                        if response3:
                            print('\nUploading file: {}'.format(fileName))
                        else:
                            tkinter.Tk().bell()
                            tkinter.Tk().bell()
                            sys.exit('\nERROR: File upload failed: {}'.format(fileName))
                        print(response3.text)
                        print('UPLOADED:'+fileName)  
                else:
                    print('******NOT UPLOADED: '+fileName)
        else:
            print("Skipped file uploads per user request.")
print("FWHM in arcsec: "+psf)
print("Obs length in minutes: "+obsdur)
print("Number of exposures: "+obsnum)
print("Mean aperture radius in pixels: "+photaprad)
print("Mean aperture radius in arcsec: "+apradarcsec)
tkinter.Tk().bell()
input("Press Enter to continue...")