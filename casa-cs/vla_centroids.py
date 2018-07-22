for s in range(0,4):
	for channel in [3,7,11,15,19,23,27,31,35,19,43,47,51,55,59,63]:
		if s == 3 and channel > 48:
			break
		timerange='19:10:28.725~19:10:28.825'
		vis = '/home/jon/Documents/sun_20141101_t191020-191040.50ms.cal.ms'
		#restoringbeam='30'         
		spw = str(s)+':'+str(channel) # select the lowest four spectral windows, excluding some channels at the edge
		specmode = 'mfs'  # use multi-frequency synthesis to take into account all selected spectral channels within window!
		# mode='channel'
		imagename = '/home/jon/Documents/images/clean/centroid_'+str(s)+'_'+str(channel)
		imsize = [200, 200] # size of the out image in pixels 
		cell = ['0.75arcsec', '0.75arcsec']  
		phasecenter_sun = 'J2000 14h26m45.9 -14d31m25.0' # use RA DEC of the sun center (from, e.g., JPL Horizons)
		phasecenter_event = 'J2000 14h26m58.183 -14d30m34.521' # event center
		stokes='I' 
		interactive = False
		tclean(vis=vis,spw=spw, specmode=specmode, phasecenter=phasecenter_event,niter=500,imagename=imagename, timerange=timerange,imsize=imsize,cell=cell,stokes=stokes,interactive=False)


for s in range(0,4):
	for channel in [3,7,11,15,19,23,27,31,35,19,43,47,51,55,59,63]:
		if s == 3 and channel > 48:
			break
		timerange='19:10:28.725~19:10:28.825'
		vis = '/home/jon/Documents/sun_20141101_t191020-191040.50ms.cal.ms'
		#restoringbeam='30'         
		spw = str(s)+':'+str(channel) # select the lowest four spectral windows, excluding some channels at the edge
		specmode = 'mfs'  # use multi-frequency synthesis to take into account all selected spectral channels within window!
		# mode='channel'
		imagename = '/home/jon/Documents/images/clean1/centroid_'+str(s)+'_'+str(channel)
		imsize = [200, 200] # size of the out image in pixels 
		cell = ['1.0arcsec', '1.0arcsec']  
		phasecenter_sun = 'J2000 14h26m45.9 -14d31m25.0' # use RA DEC of the sun center (from, e.g., JPL Horizons)
		phasecenter_event = 'J2000 14h26m58.183 -14d30m34.521' # event center
		stokes='I' 
		interactive = False
		tclean(vis=vis,spw=spw, specmode=specmode, phasecenter=phasecenter_event,niter=1,imagename=imagename, timerange=timerange,imsize=imsize,cell=cell,stokes=stokes,interactive=False)



for s in range(0,4):
	for channel in [3,7,11,15,19,23,27,31,35,19,43,47,51,55,59,63]:
		if s == 3 and channel > 48:
			break
		timerange='19:10:28.725~19:10:28.825'
		vis = '/home/jon/Documents/sun_20141101_t191020-191040.50ms.cal.ms'
		#restoringbeam='30'         
		spw = str(s)+':'+str(channel) # select the lowest four spectral windows, excluding some channels at the edge
		specmode = 'mfs'  # use multi-frequency synthesis to take into account all selected spectral channels within window!
		# mode='channel'
		imagename = '/home/jon/Documents/images/dirty/centroid_'+str(s)+'_'+str(channel)
		imsize = [200, 200] # size of the out image in pixels 
		cell = ['1.0arcsec', '1.0arcsec'] 
		phasecenter_sun = 'J2000 14h26m45.9 -14d31m25.0' # use RA DEC of the sun center (from, e.g., JPL Horizons)
		phasecenter_event = 'J2000 14h26m58.183 -14d30m34.521' # event center
		stokes='I' 
		interactive = False
		tclean(vis=vis,spw=spw, specmode=specmode, phasecenter=phasecenter_event,niter=0,imagename=imagename, timerange=timerange,imsize=imsize,cell=cell,stokes=stokes,interactive=False)



for s in range(0,4):
	for channel in [3,7,11,15,19,23,27,31,35,19,43,47,51,55,59,63]:
		if s == 3 and channel > 48:
			break
		cs_alg='peak'
		timerange='19:10:28.725~19:10:28.825'
		vis = '/home/jon/Documents/sun_20141101_t191020-191040.50ms.cal.ms'
		#restoringbeam='30'         
		spw = str(s)+':'+str(channel) # select the lowest four spectral windows, excluding some channels at the edge
		specmode = 'mfs'  # use multi-frequency synthesis to take into account all selected spectral channels within window!
		# mode='channel'
		imagename = '/home/jon/Documents/images/peak/centroid_'+str(s)+'_'+str(channel)
		imsize = [100, 100] # size of the out image in pixels 
		cell = ['1.0arcsec', '1.0arcsec'] 
		phasecenter_sun = 'J2000 14h26m45.9 -14d31m25.0' # use RA DEC of the sun center (from, e.g., JPL Horizons)
		

		phasecenter_event = 'J2000 14h26m58.183 -14d30m34.521' # event center
		stokes='I' 
		interactive = False
		p7_cs(vis=vis,spw=spw, specmode=specmode, phasecenter=phasecenter_event,niter=0,imagename=imagename, timerange=timerange,imsize=imsize,cell=cell,stokes=stokes,interactive=False, cs_alg=cs_alg, psf_threshold=0.02, psf_cutoff=False)













