clean(vis='/home/jon/Documents/SNR_G55_10s.calib.ms', imagename='/home/jon/Documents/images/supernova/SNR_G55_10s.128.raw',
      imsize=128, cell='26arcsec', niter=0, interactive=False)

tclean(vis='/home/jon/Documents/SNR_G55_10s.calib.ms', imagename='/home/jon/Documents/images/supernova/SNR_G55_10s.128.raw',
      imsize=128, cell='26arcsec', niter=0, threshold='0.12mJy', stokes='I',interactive=False, wprojplanes = 1, pblimit=-1.0)
tclean(vis='/home/jon/Documents/SNR_G55_10s.calib.ms', imagename='/home/jon/Documents/images/supernova/SNR_G55_10s.128.clean',
      imsize=128, cell='26arcsec', niter=1000, threshold='0.12mJy', stokes='I',interactive=False, wprojplanes = 1, pblimit=-1.0)

cs_alg='positive_deconv'
vis = '/home/jon/Documents/SNR_G55_10s.calib.ms'
imagename='/home/jon/Documents/images/supernova/SNR_G55_10s.128.'+cs_alg
p7_cs(vis=vis,niter=0,imagename=imagename, imsize=[128,128], cell='26arcsec', interactive=False, stokes='I', cs_alg=cs_alg, psf_threshold=0.02, psf_cutoff=False, wprojplanes = 1, pblimit=-1.0)

algorithms = ['L1', 'L2','TV','haar','starlets']
for cs_alg in algorithms:
	vis = '/home/jon/Documents/SNR_G55_10s.calib.ms'
	imagename='/home/jon/Documents/images/supernova/SNR_G55_10s.128.'+cs_alg

	lamb_est_path = '/home/jon/Documents/images/supernova/SNR_G55_10s.128.positive_deconv'
	lamb_est = (lamb_est_path+'_objectiveVal.csv', lamb_est_path+'_solution.csv')

	p7_cs(vis=vis,niter=0,imagename=imagename, imsize=[128,128], cell='26arcsec', interactive=False, stokes='I', cs_alg=cs_alg, psf_threshold=0.02, psf_cutoff=False, lambda_estimate= lamb_est, wprojplanes = 1, pblimit=-1.0)


cs_alg='starlets'
vis = '/home/jon/Documents/SNR_G55_10s.calib.ms'
imagename='/home/jon/Documents/images/supernova/SNR_G55_10s.128.'+cs_alg+"2"
lamb_est_path = '/home/jon/Documents/images/supernova/SNR_G55_10s.128.positive_deconv'
lamb_est = (lamb_est_path+'_objectiveVal.csv', lamb_est_path+'_solution.csv')
p7_cs(vis=vis,niter=0,imagename=imagename, imsize=[128,128], cell='26arcsec', interactive=False, stokes='I', cs_alg=cs_alg, psf_threshold=0.02, psf_cutoff=False, lambda_estimate= lamb_est, wprojplanes = 1, pblimit=-1.0, starlet_levels=2)

cs_alg='starlets'
vis = '/home/jon/Documents/SNR_G55_10s.calib.ms'
imagename='/home/jon/Documents/images/supernova/SNR_G55_10s.128.'+cs_alg+"3"
lamb_est_path = '/home/jon/Documents/images/supernova/SNR_G55_10s.128.positive_deconv'
lamb_est = (lamb_est_path+'_objectiveVal.csv', lamb_est_path+'_solution.csv')
p7_cs(vis=vis,niter=0,imagename=imagename, imsize=[128,128], cell='26arcsec', interactive=False, stokes='I', cs_alg=cs_alg, psf_threshold=0.02, psf_cutoff=False, lambda_estimate= lamb_est, wprojplanes = 1, pblimit=-1.0, starlet_levels=3)




cs_alg='L1+L2'
vis = '/home/jon/Documents/SNR_G55_10s.calib.ms'
imagename='/home/jon/Documents/images/supernova/SNR_G55_10s.128.'+cs_alg
lamb_est_path = '/home/jon/Documents/images/supernova/SNR_G55_10s.128.positive_deconv'
lamb_est = (lamb_est_path+'_objectiveVal.csv', lamb_est_path+'_solution.csv')
p7_cs(vis=vis,niter=0,imagename=imagename, imsize=[128,128], cell='26arcsec', interactive=False, stokes='I', cs_alg=cs_alg, psf_threshold=0.02, psf_cutoff=False, lambda_estimate= lamb_est, wprojplanes = 1, pblimit=-1.0)









