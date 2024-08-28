"""
geodesic distances
"""
import os
import numpy as np
from scipy.signal import correlate
import pygeodesic.geodesic as geodesic
import antropy

def _clean_mesh_txt(filename):
    
    fobj = open(filename)
    mesh = fobj.read().replace(" ", "")
    mesh = mesh.replace(",", " ")
    fobj.close()
    
    os.remove(filename)
    fobj = open(filename, 'w')
    fobj.write(mesh)
    fobj.close()


def _get_5mm_patch(src, trg_mesh):
    
    points, faces = trg_mesh
    geoalg = geodesic.PyGeodesicAlgorithmExact(points, faces)

    trg = None
    distance, path = geoalg.geodesicDistances(np.array([src]), trg)
    
    patch = np.where(distance < 5)[0]
    
    return patch


def _compute_intrinsic_timescale(ts_array, tr):
    # assume ts_array is a vectorized timeseries
    # returns a single float

    nt = len(ts_array)
    # timeseries was mean centered in preprocessing script

    autocorr = correlate(ts_array, ts_array, 'full', 'direct')
    autocorr = autocorr / max(autocorr) # unit scale
    autocorr = autocorr[nt + 1:]

    # splice the autocorrelation up to first negative autocorrelation
    stop = np.where(autocorr < 0)[0][0]
    intrinsic_timescale = np.sum(autocorr[:stop])

    return intrinsic_timescale * tr


def _compute_rmssd(ts_array):
    
    nt = len(ts_array)
    
    summate = 0
    
    for i in range(nt - 1):
        diff = ts_array[i + 1] - ts_array[i]
        summate += diff
    
    rmssd = summate / (nt - 1)
    
    return rmssd

def _compute_entropy(ts_array):

    ent = antropy.sample_entropy(ts_array)

    return ent

    
class Channel:
    
    def __init__(self, subject, electrode_name, contact_pair, \
                 closest_vertex, ts_idx, 
                 patch_5k, patch_hipp, \
                 structure, coords, \
                 region=None, rr=None, frr=None, \
                 abnormal95=None, abnormal90=None, \
                 ieeg_ts_clean = None, fmri_ts_clean = None):
       
       self.subject = subject
       self.electrode_name = electrode_name
       self.contact_pair = contact_pair
       self.closest_vertex = closest_vertex
       self.coords = coords
       self.ts_idx = ts_idx
       self.patch_5k = patch_5k
       self.patch_hipp = patch_hipp
       self.structure = structure
       
       self.region = region
       self.rr = rr
       self.frr = frr
       self.abnormal95= abnormal95
       self.abnormal90= abnormal90
       
       self.ieeg_ts_clean = ieeg_ts_clean
       self.fmri_ts_clean = fmri_ts_clean
       
       
    def update_hfo(self, region, rr, frr, abnormal95, abnormal90):
        self.region = region
        self.rr = rr
        self.frr = frr
        self.abnormal95 = abnormal95
        self.abnormal90 = abnormal90
    
    def update_ieeg_ts(self, ieeg_ts_clean):
        self.ieeg_ts_clean = ieeg_ts_clean
        
    def update_fmri_ts(self, fmri_ts_clean):
        self.fmri_ts_clean = fmri_ts_clean        