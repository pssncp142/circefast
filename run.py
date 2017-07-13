import pycircefast as cf
import time 

conf = cf.init()
"""
conf.path = " "
conf.f_dark = "darks/dark_band_0341_1706_02_005_04869.fits"
conf.f_badp = "darks/badpix_band_0341_1706_02_005_04869.fits"
conf.f_flat = "darks/flat_35.fits"
conf.f_name = "data/CIRCE2017-05-23-%04d.fits"
conf.band_ndx[0] = 341
conf.band_ndx[1] = 1707
conf.seq = 1
"""
conf.st = 84
conf.seq = 2

#print "here.."

t1 = time.time()

cf.darksub(conf)
cf.doskies(conf)
cf.skysub(conf)
cf.fftcorr(conf)
cf.flatdivide(conf)
cf.badpixrem(conf)
cf.tile(conf)

print time.time()-t1
