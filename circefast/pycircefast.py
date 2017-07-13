from ctypes import *

libcf = cdll.LoadLibrary("./circefast/circefast.so")


class config(Structure):
    _fields_ = \
                 [('path'  , c_char*100),
                  ('d_flat', c_char*100),
                  ('f_flat', c_char*100),
                  ('f_name', c_char*100),
                  ('f_dark', c_char*100),
                  ('f_badp', c_char*100),
                  ('hwp'   , c_char*5),
                  ('filt'  , c_char*2),
                  ('n_dith', c_int),
                  ('f_arr' , c_int*2),
                  ('band_ndx', c_int*2),
                  ('add_back', c_int),
                  ('fft_coor', c_int),
                  ('seq'     , c_int),
                  ('ndx'     , c_int),
                  ('st'      , c_int)]

def darksub(cnf):
    return libcf.darksub_all(cnf)

def doskies(cnf):
    return libcf.doskies_all(cnf)
def skysub(cnf):
    return libcf.skysub_all(cnf)
def fftcorr(cnf):
    return libcf.fftcorr_all(cnf)
def flatdivide(cnf):
    return libcf.flatdivide_all(cnf)
def badpixrem(cnf):
    return libcf.badpixrem_all(cnf)
def tile(cnf):
    return libcf.tile_all(cnf)


def init():

    conf = config()

    conf.path = ""
    conf.f_dark = "/home/ydallilar/software/circefast/darks/dark_band_0341_1706_02_005_04869.fits"
    conf.f_badp = "darks/badpix_band_0341_1706_02_005_04869.fits"
    conf.f_flat = "darks/flat_35.fits"
    conf.f_name = "data/CIRCE2017-05-23-%04d.fits"
    conf.seq = 1
    conf.n_dith = 5
    conf.band_ndx[0] = 341
    conf.band_ndx[1] = 1707
   
    return conf

    
    
