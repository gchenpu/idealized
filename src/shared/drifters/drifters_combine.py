#!/usr/bin/env python

VERSION = "$Id$"

from Scientific.IO.NetCDF import NetCDFFile
import sys
import re
import Numeric as Num
import glob
from optparse import OptionParser


class drifters_combine:

    def __init__(self, filenames):

        self.data = {}
        self.nd = None
        self.nf = None

        self.global_atts = {}
        self.position_atts = {}
        self.field_atts = {}
        
        for f in filenames:
            nc =  NetCDFFile(f, 'r')
            index_time= nc.variables['index_time'][:]
            time      = nc.variables['time'][:]
            ids       = nc.variables['ids'][:]

            positions = nc.variables['positions'][:]
            nd = Num.shape(positions)[1]

            # get global attributes
            for a in nc.__dict__:
                self.global_atts[a] = nc.__dict__[a]
            
            # get attributes
            for a in nc.variables['positions'].__dict__:
                if not self.position_atts.has_key(a):
                    self.position_atts[a] = getattr(nc.variables['positions'], a)

            fields    = None
            nf        = 0
            if 'fields' in nc.variables:
                fields    = nc.variables['fields'][:]
                nf = Num.shape(fields)[1]
                for a in nc.variables['fields'].__dict__:
                    if not self.field_atts.has_key(a):
                        self.field_atts[a] = getattr(nc.variables['fields'], a)

            if self.nf==None: self.nf = nf
            if self.nd==None: self.nd = nd

            if nf != self.nf or nd != self.nd:
                raise 'Incompatible number of fields (nf) or space dimensions (nd) in file'% f

            for i in range(len(ids)):
                id = ids[i]
                it = index_time[i]
                tim= time[i]
                xyz= positions[i,:]
                fld = ()
                if nf>0: fld= fields[i,:]

                if not self.data.has_key(it):
                    self.data[it] = {}

                if not self.data[it].has_key(id):
                    self.data[it][id] = (tim,) + tuple(xyz) + tuple(fld)

    def save(self, outfile):

        nc = NetCDFFile(outfile,'w')
        
        nc.createDimension('it_id', None)
        nc.createDimension('nd', self.nd)
        if self.nf>0:
            nc.createDimension('nf', self.nf)
        
        nc_index_time = nc.createVariable('index_time', Num.Int, ('it_id',))
        nc_time       = nc.createVariable('time', Num.Float64, ('it_id',))
        nc_ids        = nc.createVariable('ids', Num.Int, ('it_id',))
        nc_positions  = nc.createVariable('positions', Num.Float64, \
                                          ('it_id', 'nd'))

        for a in self.global_atts:
            setattr(nc, a, self.global_atts[a])

        for a in self.position_atts:
            setattr(nc_positions, a, self.position_atts[a])

        if self.nf>0:           
            nc_fields     = nc.createVariable('fields', Num.Float64, \
                                          ('it_id', 'nf'))
            for a in self.field_atts:
                setattr(nc_fields, a, self.field_atts[a])

        k = 0
        for it in range(len(self.data)):
            print 'it=', it
            for id in self.data[it]:
                print '\tid=', id
                nc_index_time[k] = it
                nc_time[k]       = self.data[it][id][0]
                nc_ids[k]        = id
                nc_positions[k:k+1,:]= self.data[it][id][1:1+self.nd]
                if self.nf>0:
                    nc_fields[k:k+1,:]= self.data[it][id][1+self.nd:1+self.nd+self.nf]
                k += 1
               
        nc.close()
                    
                
###############################################################################



def main():

    parser = OptionParser(version=VERSION)
    parser.add_option('-p', '--pe_range', action='store', type="string",
                  dest="pe_range",
                  help='PE range (e.g. "10...20"). By default the output of all PEs will be combined.',
                  default='',
                  )
    parser.add_option('-f', '--file', action='store', type="string",
                  dest="filename",
                  help='Input file. Files FILE.XXXX to FILE.YYYY will be combined.',
                  default='drifters_out.nc',
                  )
    parser.add_option('-o', '--output', action='store', type="string",
                  dest="output",
                  help='The combined files will be saved as OUTPUT. By default OUTPUT=FILE.',
                  default='',
                  )
    
    options, args = parser.parse_args(sys.argv)
    fname = options.filename
    if not fname:
        print 'ERROR: must supply a netcdf file name ([-f|--file] FILE).'
        sys.exit(1)

    min_pe = 0
    if options.pe_range:
        p1 = re.search(r'^\s*(\d+)\s*\.', options.pe_range)
        p2 = re.search(r'\.\s*(\d+)\s*$', options.pe_range)
        if p1 and p2:
            min_pe = p1.group(1)
            max_pe = p2.group(1)
        else:
            print 'ERROR: pe range should be in format [-p|--pe_range] "min_pe..max_pe".'
            sys.exit(1)
            
        fnames = [fname + '.' + str(pe).zfill(4) for pe in range(int(min_pe), int(max_pe)+1)]
        
    else:
        fnames = glob.glob(fname + '.' + '????')

    vfc = drifters_combine(fnames)
    if options.output:
        ofile = options.output
    else:
        ofile = fname

    vfc.save(ofile)
    
if __name__=='__main__': main()
