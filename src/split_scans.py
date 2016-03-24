#! /usr/bin/python 

import math
import copy
import gc
#import numpy as np

from time_stream import base_single
from kiyopy import utils
from core import data_block
from core import fitsGBT

from mpi4py import MPI

class SplitScans(base_single.BaseSingle):

    prefix = 'ss_'

    params_init = {
            'scan_length' : 1812,
            }

    def mpiexecute(self, n_processes=1):
        """
        Process all data with MPI
        To start with MPI, you need to change manager.py 
        calling mpiexecute instead of execute. 
        and do,

        $ mpirun -np 9 --bynode  python  manager.py pipeline.pipe

        """

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        params = self.params
        n_files = len(params['file_middles'])
        for file_ind in range(n_files)[rank::size]:

            self.process_file(file_ind)

        comm.barrier()

    
    def process_file(self, file_ind):

        params = self.params

        scan_length = params['scan_length']

        file_middle = params['file_middles'][file_ind]
        input_fname = (params['input_root'] + file_middle +
                       params['input_end'])
        output_fname = (params['output_root']
                        + file_middle + params['output_end'])
        Writer = fitsGBT.Writer(feedback=self.feedback)
        Reader = fitsGBT.Reader(input_fname, feedback=self.feedback)
        scan_inds = params['scans']
        if len(scan_inds) == 0 or scan_inds is None :
            scan_inds = range(len(Reader.scan_set))
        # Loop over scans.
        jj = 0
        for thisscan in scan_inds :
            Blocks = Reader.read(thisscan, params['IFs'], force_tuple=True)
            for Data in Blocks:
                self.action(Data, Writer)
            del Blocks
            gc.collect()

        # Go to a new line if we are printing statistics.
        if hasattr(self, 'feedback_title') and self.feedback > 1:
            print ''
        # Finally write the data back to file.
        utils.mkparents(output_fname)
        Writer.write(output_fname)



    def action(self, Data, Writer):
        params = self.params
        Blocks = split_scans(Data, params['scan_length'] )
        for Block in Blocks:
            Writer.add_data(Block)
        del Blocks
        gc.collect()

        #return Blocks

def split_scans(Data, scan_length):

    scan_number = int(math.ceil(Data.data.shape[0]/float(scan_length)))
    out_Blocks = ()
    print "split into %d scans.."%scan_number
    for ii in range(scan_number):
        start_ind = ii*scan_length
        end_ind = (ii+1)*scan_length
        new_data = Data.data[start_ind:end_ind, ...]
        this_Block = data_block.DataBlock(new_data)
        del new_data
        this_Block.field = copy.deepcopy(Data.field)
        this_Block.field_axes = copy.deepcopy(Data.field_axes)
        this_Block.field_formats = copy.deepcopy(Data.field_formats)

        # fix the scan length
        this_Block.field['EXPOSURE'] = Data.field['EXPOSURE'][start_ind:end_ind, ...]
        this_Block.field['DATE-OBS'] = Data.field['DATE-OBS'][start_ind:end_ind, ...]
        #this_Block.field['LST']      = Data.field['LST'][start_ind:end_ind, ...]
        #this_Block.field['ELEVATION']= Data.field['ELEVATION'][start_ind:end_ind, ...]
        #this_Block.field['AZIMUTH']  = Data.field['AZIMUTH'][start_ind:end_ind, ...]
        #this_Block.field['OBSFREQ']  = Data.field['OBSFREQ'][start_ind:end_ind, ...]
        this_Block.field['RA']       = Data.field['RA'][start_ind:end_ind, ...]
        this_Block.field['DEC']      = Data.field['DEC'][start_ind:end_ind, ...]
        this_Block.field['CRVAL2']   = Data.field['CRVAL2'][start_ind:end_ind, ...]
        this_Block.field['CRVAL3']   = Data.field['CRVAL3'][start_ind:end_ind, ...]
        #print "scan number :", this_Block.field['SCAN'],
        this_Block.field['SCAN']   += ii
        #print "scan number :", this_Block.field['SCAN']
        this_Block.add_history('Split scans into different block ', 
                         ('scan_length %d'%scan_length))

        this_Block.verify()
        out_Blocks += (this_Block, )
        del this_Block
        gc.collect()

    return out_Blocks

#class SplitScans(base_single.BaseSingle):
#
#    prefix = 'ss_'
#
#    params_init = {
#            'scan_length' : 1812,
#            }
#
#    def process_file(self, file_ind):
#
#        params = self.params
#
#        scan_length = params['scan_length']
#
#        file_middle = params['file_middles'][file_ind]
#        input_fname = (params['input_root'] + file_middle +
#                       params['input_end'])
#        Reader = fitsGBT.Reader(input_fname, feedback=self.feedback)
#        scan_inds = params['scans']
#        if len(scan_inds) == 0 or scan_inds is None :
#            scan_inds = range(len(Reader.scan_set))
#        # Loop over scans.
#        jj = 0
#        for thisscan in scan_inds :
#            Blocks = Reader.read(thisscan, params['IFs'], force_tuple=True)
#            for Data in Blocks:
#                scan_number = int( math.ceil( Data.data.shape[0]/float(scan_length)))
#                for ii in range(scan_number):
#                    print Data.data.shape
#                    NewData = split_scans(Data, ii*scan_length, (ii+1)*scan_length)
#                    if not type(NewData) is tuple:
#                        NewData = (NewData,)
#                    Writer = fitsGBT.Writer(feedback=self.feedback)
#                    Writer.add_data(NewData)
#                    output_fname = (params['output_root']
#                        + file_middle + '_%02d'%jj + params['output_end'])
#                    utils.mkparents(output_fname)
#                    Writer.write(output_fname)
#                    jj += 1
#
#                    del NewData
#                    del Writer
#                    gc.collect()
#
#
#def split_scans(Data, start_ind, end_ind):
#
#    this_Block = data_block.DataBlock(Data.data[start_ind:end_ind, ...])
#    this_Block.field = copy.deepcopy(Data.field)
#    this_Block.field_axes = copy.deepcopy(Data.field_axes)
#    this_Block.field_formats = copy.deepcopy(Data.field_formats)
#
#    print this_Block.data.shape
#
#    # fix the scan length
#    this_Block.field['EXPOSURE'] = Data.field['EXPOSURE'][start_ind:end_ind, ...]
#    this_Block.field['DATE-OBS'] = Data.field['DATE-OBS'][start_ind:end_ind, ...]
#    #this_Block.field['LST']      = Data.field['LST'][start_ind:end_ind, ...]
#    #this_Block.field['ELEVATION']= Data.field['ELEVATION'][start_ind:end_ind, ...]
#    #this_Block.field['AZIMUTH']  = Data.field['AZIMUTH'][start_ind:end_ind, ...]
#    #this_Block.field['OBSFREQ']  = Data.field['OBSFREQ'][start_ind:end_ind, ...]
#    this_Block.field['RA']       = Data.field['RA'][start_ind:end_ind, ...]
#    this_Block.field['DEC']      = Data.field['DEC'][start_ind:end_ind, ...]
#    this_Block.field['CRVAL2']   = Data.field['CRVAL2'][start_ind:end_ind, ...]
#    this_Block.field['CRVAL3']   = Data.field['CRVAL3'][start_ind:end_ind, ...]
#
#    this_Block.verify()
#    return this_Block

if __name__=='__main__':
    pass
