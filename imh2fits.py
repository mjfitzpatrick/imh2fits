#!/bin/env python
#
#  IMH2FITS -- Convert IRAF OIF images in ".imh" format to FITS
#
#  Convert a 1-D or 2-D IRAF OIF image (i.e. a ".imh" file and the associated
#  ".pix" file) to FITS.  This task will read both old V1 files written in a
#  native byte-order as well as big-endian V2 files and swap bytes as needed.
#  IMH2FITS can be used to recover old IRAF image data and convert them to 
#  FITS files when 
# 
#     1) no IRAF installation is currently available
#     2) old IRAF image files cannot be read because of byte-swapping
# 
#  Output FITS files are written to the current directory using a filename
#  derived from the ".imh" image filename.
# 
#  Usage:
#         % imh2fits [<opts>] <files>
# 
#     where <opts> include:
#         -d              print debug output
#         -h              this help message
#         -l              list file rather than convert
#         -s              force byte swapping
#         -v              print verbose output
#         -w              print FITS validation warning messages
# 
#         <files>         a list of ".imh" file to convert
# 
#     The task will automatically search for a ".pix" pixel file that has the
#     same root name as the ".imh" file.  The output FITS file will be written
#     to the current directory using the root filename and a ".fits" extension.
# 
#  Examples:
#     1) Convert all ".imh" files in the current directory to FITS:
# 
#         % python imh2fits -v *.imh
# 
#        The '-v' flag will print some minimal information for each converted
#        file as they are processed.


from __future__ import print_function

import os
import sys
import ast
import struct
import fnmatch
import warnings
import numpy as np
import astropy.io.fits as fits
from astropy.io.fits.verify import VerifyWarning

__author__ = 'Mike Fitzpatrick <mjfitzpatrick@gmail.com>'
__version__ = '2020-07-15'      # yyyy-mm-dd


# Global options.
debug = False
verbose = False
force_swap = False
warn = False
list_only = False


# OIF image file types
V1HDR = 1                               # OIF v1 .imh file
V1PIX = 2                               # OIF v1 .pix file
V2HDR = 3                               # OIF v2 .imh file
V2PIX = 4                               # OIF v2 .pix file

# OIF image file magic strings
V1HDR_MAGIC_BIG_ENDIAN = b'\x00i\x00m\x00h\x00d\x00r'
V1PIX_MAGIC_BIG_ENDIAN = b'\x00i\x00m\x00p\x00i\x00x'
V1HDR_MAGIC_LITTLE_ENDIAN = b'i\x00m\x00h\x00d\x00r\x00'
V1PIX_MAGIC_LITTLE_ENDIAN = b'i\x00m\x00p\x00i\x00x\x00'
V2HDR_MAGIC = b'imhv2'
V2PIX_MAGIC = b'impv2'

# IRAF pixtype codes
TY_SHORT = 3                            # signed 16-bit integer
TY_INT = 4                              # signed 32-bit integer
TY_LONG = 5                             # signed 32-bit integer
TY_REAL = 6                             # 32-bit floating point
TY_DOUBLE = 7                           # 64-bit floating point
TY_USHORT = 11                          # unsigned 16-bit integer

pixtypes = {
    TY_SHORT: 'i2',                    # (native) short
    TY_USHORT: 'U2',                   # (native) unsigned short
    TY_INT: 'i4',                      # (native) int
    TY_LONG: 'i4',                     # (native) long
    TY_REAL: 'f4',                     # (native) real
    TY_DOUBLE: 'f8',                   # (native) double
}
oiftypes = {
    TY_SHORT: 'short',                 # (native) short
    TY_USHORT: 'ushort',               # (native) unsigned short
    TY_INT: 'int',                     # (native) int
    TY_LONG: 'long',                   # (native) long
    TY_REAL: 'real',                   # (native) real
    TY_DOUBLE: 'double',               # (native) double
}

# For converting to the IRAF 1/1/1980 epoch used in timestamps (approx).
EPOCH_1980 = (10 * (3600 * 24 * 365.25))


class OIFImage(object):
    '''A Class for OIF images.
    '''

    def __init__(self, hdrfile, pixfile):
        self.path = hdrfile
        self.base = os.path.basename(hdrfile)
        self.root = os.path.splitext(self.base)[0]

        self.hdrfile = OIFHeader(hdrfile)
        self.pixfile = OIFPixfile(pixfile)

    def fileType(self, data):
        '''Determine the OIF File type.  Do this for both header and
           pixel files.
        '''
        magic = data[:10]

        if magic == V1HDR_MAGIC_BIG_ENDIAN:
            return (V1HDR, sys.byteorder != 'big')
        elif magic == V1PIX_MAGIC_BIG_ENDIAN:
            return (V1PIX, sys.byteorder != 'big')

        elif magic == V1HDR_MAGIC_LITTLE_ENDIAN:
            return (V1HDR, sys.byteorder != 'little')
        elif magic == V1PIX_MAGIC_LITTLE_ENDIAN:
            return (V1PIX, sys.byteorder != 'little')

        elif magic[:5] == V2HDR_MAGIC:
            return (V2HDR, sys.byteorder == 'little')
        elif magic[:5] == V2PIX_MAGIC:
            return (V2PIX, sys.byteorder == 'little')

        return None, None

    def toFITS(self, path, verbose=False):
        '''Save the header/pixfile to a FITS image.
        '''
        # Create the FITS object
        hdu = fits.PrimaryHDU()

        # Set the pixel array.  Note this is transposed relative to the
        # normal FITS ordering.
        if self.pixfile.ndim == 1:
            hdu.data = self.pixfile.pixels
        else:
            hdu.data = self.pixfile.pixels.T

        # Set the PHU header keywords from the OIF user area.
        for key in self.hdrfile.keywords:
            if key.strip() == '':
                continue
            st = 7 if key[:7] in ['HISTORY', 'COMMENT'] else 8
            k = key[:st].strip().replace('/', '')
            if key[st+1:].count("'") == 2:
                v = "'" + key[st+1:].split("'")[1] + "'"
                c = key[key[st+1:].find('/')].strip()
            elif key.find('/') > 0:
                v = key[st+1:].split('/')[0].strip()
                c = key[st:].split('/')[1].strip()
            else:
                v = key[st+1:].strip()
                c = ''

            if debug:
                print('k :%s:  v :%s:  c :%s:' % (str(k), str(v), str(c)))

            # Determine the keyword value type so we don't write everything
            # as a string.
            if key[:7] not in ['HISTORY', 'COMMENT']:
                if v.lower() not in ['t','f','true','false']:
                    t = ast.literal_eval(v)
                    if type(t) == int:
                        v = int(v)
                    elif type(t) == float:
                        v = float(v)
                    else:
                        v = v.strip("'")

            hdu.header[k] = (v, c)

        # Add a HISTORY keyword to track back to the original file and
        # write to the requested path.
        hdu.header.add_history('New copy of ' + str(self.hdrfile.base))

        # Print some image statistics to see whether things look reasonable.
        if verbose or list_only:
            title = ''
            if self.hdrfile.title != '':
                title = self.hdrfile.title.decode('utf-8')
            if self.pixfile.ndim == 1:
                print('%20s[%d][%s]\t%s\n\t\tmean: %g  min: %g  max: %g' %
                      (path, self.pixfile.dims[0], 
                      oiftypes[self.pixfile.pixtype], title,
                      np.mean(self.pixfile.pixels),
                      np.min(self.pixfile.pixels),
                      np.max(self.pixfile.pixels)))
            else:
                print('%20s[%d,%d][%s]\t%s\n\t\tmean: %g  min: %g  max: %g' %
                      (path, self.pixfile.dims[0],  self.pixfile.dims[1], 
                      oiftypes[self.pixfile.pixtype], title,
                      np.mean(self.pixfile.pixels),
                      np.min(self.pixfile.pixels),
                      np.max(self.pixfile.pixels)))

        # Write out the FITS file.
        if not list_only:
            hdu.writeto(path, overwrite=True)


class OIFHeader(OIFImage):
    '''A Class for OIF image header files.
    '''

    def __init__(self, path):
        self.path = path
        self.base = os.path.basename(path)
        self.extn = os.path.splitext(self.base)[1]

        if self.extn != '.imh':
            raise Exception('Unknown image header extension')

        # Load the file contents to memory.
        with open(path, mode='rb') as fd:
            self.content = fd.read()

        # Get the file type and whether we need to swap it.
        self.type, self.swap = self.fileType(self.content)
        ifmt = '>i' if (sys.byteorder == 'little' and self.swap) else '<i'
        if force_swap:
            ifmt = '<i' if ifmt == '>i' else '>i'

        self.data = []
        self.userarea = []
        self.title = ''
        if self.type == V1HDR:
            # Parse the V1 resident disk header.  V1 images use SPP strings
            # we convert to 'normal' strings by removing the NULLs.
            self.hdr = self.content[:2052]
            self.userarea = self.content[2052:].replace(b'\x00', b'')
            self.keywords = self.userarea.split(b'\n')[:-1]
            for i, k in enumerate(self.keywords):
                self.keywords[i] = self.keywords[i].decode('utf-8')
            self.title = self.hdr[732:888].split(b'\0\0')[0].replace(b'\x00', b'')
            self.pixtype = struct.unpack(ifmt, self.hdr[16:20])[0]

            history = self.content[892:2048].replace(b'\x00', b'')
            self.history = b''.join(history.split(b'\n'))

            self.ndim = struct.unpack(ifmt, self.hdr[20:24])[0]
            self.dims = list(struct.unpack(ifmt+'iiiiii', self.hdr[24:52]))
            self.pdims = list(struct.unpack(ifmt+'iiiiii', self.hdr[52:80]))

            self.date = struct.unpack(ifmt, self.hdr[108:112])[0] + EPOCH_1980

        elif self.type == V2HDR:
            # Parse the V2 resident disk header.  V2 images use normal
            # "packed" strings.
            self.hdr = self.content[:2046]
            self.userarea = self.content[2046:].replace(b'\x00', b'')
            self.keywords = self.userarea.split(b'\n')
            for i, k in enumerate(self.keywords):
                self.keywords[i] = self.keywords[i].decode('utf-8')
            self.title = self.hdr[638:718].replace(b'\x00', b'')
            self.pixtype = struct.unpack(ifmt, self.hdr[10:14])[0]

            history = self.content[990:2048].replace(b'\x00', b'')
            self.history = b''.join(history.split(b'\n'))

            self.ndim = struct.unpack(ifmt, self.hdr[18:22])[0]
            self.dims = list(struct.unpack(ifmt+'iiiiii', self.hdr[22:50]))
            self.pdims = list(struct.unpack(ifmt+'iiiiii', self.hdr[50:78]))

            self.date = struct.unpack(ifmt, self.hdr[106:110])[0] + EPOCH_1980

        self.printHdr()

    def printHdr(self):
        if debug:
            print('==============\nhdrfile;: ' + path)
            print('type: ' + str(self.type) + '  swap: ' + str(self.swap))
            print('title: "' + str(self.title) + '"')
            print('len hdr: ' + str(len(self.hdr)))
            print('len user: ' + str(len(self.userarea)))
            print('len data: ' + str(len(self.data)))
            print('pixtype = ' + str(self.pixtype))
            print('date = ' + str(self.date))
            print('ndim = ' + str(self.ndim))
            print('dims: ' + str(self.dims) + '  phys: ' + str(self.pdims))
            print('history: "' + str(self.history) + '"')
            if verbose:
                print(self.keywords)
                print(str(type(self.keywords)))


class OIFPixfile(OIFImage):
    '''A Class for OIF image pixel files.
    '''

    def __init__(self, path):
        self.path = path
        self.base = os.path.basename(path)
        self.extn = os.path.splitext(self.base)[1]

        if self.extn != '.pix':
            raise Exception('Unknown image pixfile extension')

        # Load the file contents to memory.
        with open(path, mode='rb') as fd:
            self.content = fd.read()

        # Get the file type and whether we need to swap it.
        self.type, self.swap = self.fileType(self.content)
        ifmt = '>i' if (sys.byteorder == 'little' and self.swap) else '<i'
        if force_swap:
            ifmt = '<i' if ifmt == '>i' else '>i'

        self.data = []
        if self.type == V1PIX:
            self.hdr = self.content[:1024]
            self.swapped = self.swap
            self.pixtype = struct.unpack(ifmt, self.hdr[16:20])[0]
            self.ndim = struct.unpack(ifmt, self.hdr[20:24])[0]
            self.dims = list(struct.unpack(ifmt+'iiiiii', self.hdr[24:52]))
            self.pdims = list(struct.unpack(ifmt+'iiiiii', self.hdr[52:80]))

            self.data = self.content[1024:]

        elif self.type == V2PIX:
            self.hdr = self.content[:2048]
            self.swapped = (struct.unpack(ifmt, self.hdr[14:18])[0] == 1)
            self.pixtype = struct.unpack(ifmt, self.hdr[10:14])[0]
            self.ndim = struct.unpack(ifmt, self.hdr[18:22])[0]
            self.dims = list(struct.unpack(ifmt+'iiiiii', self.hdr[22:50]))
            self.pdims = list(struct.unpack(ifmt+'iiiiii', self.hdr[50:78]))

            self.data = self.content[2048:]

        # Convert the pixels to a Numpy array.
        self.readPixels()
        self.printHdr()

    def readPixels(self):
        '''Read the pixels from the file contents and return a Numpy array.
        '''
        if self.type == V1PIX:
            swap = '>' if (sys.byteorder == 'little' and self.swap) else '<'
        else:
            swap = '<' if (sys.byteorder == 'little' and self.swapped) else '>'
        if force_swap:
            swap = '<' if swap == '>' else '>'
        ptype = swap + pixtypes[self.pixtype]

        raw_pixels = np.frombuffer(self.data, np.dtype(ptype), -1, 0)

        # Reshape the array and extract the pixels from the physical storage.
        if self.ndim == 1:
            self.pixels = raw_pixels[:self.dims[0]]
        elif self.ndim == 2:
            phys_pixels = np.reshape(raw_pixels, (self.pdims[1], -1)).T
            self.pixels = phys_pixels[:self.dims[0], :self.dims[1]]
        else:
            raise Exception ('Error: Images >2 dimensions not supported')

        if debug:
            print('pix type: %s shape: %s' % (ptype, str(self.pixels.shape)))
            print('final min/max: %g  %g' %
                  (np.min(self.pixels), np.max(self.pixels)))

    def printHdr(self):
        if debug:
            print('==============\npixfile;: ' + path)
            print('type: ' + str(self.type) + '  swap: ' + str(self.swap))
            print('pixtype = ' + str(self.pixtype))
            print('ndim = ' + str(self.ndim))
            print('dims: ' + str(self.dims) + '  phys: ' + str(self.pdims))


def findFile(name, path):
    '''Find the first match for the named file in the specified directory path.
    '''
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)


def findAll(pattern, path):
    '''Find all files matching the given pattern in the named directory.
    '''
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


def Usage():
    '''Print a help summary.
    '''
    print('''
Usage:
       % imh2fits [<opts>] <files>

   where <opts> include:
       -d              print debug output
       -h              this help message
       -l              list file rather than convert
       -s              force byte swapping
       -v              print verbose output
       -w              print FITS validation warning messages

       <files>         a list of ".imh" file to convert

   The task will automatically search for a ".pix" pixel file that has the
   same root name as the ".imh" file.  The output FITS file will be written
   to the current directory using the root filename and a ".fits" extension.

Examples:
   1) Convert all ".imh" files in the current directory to FITS:

       % python imh2fits -v *.imh

      The '-v' flag will print some minimal information for each converted
      file as they are processed.
    ''')


############################
# Task main()
############################

if __name__ == '__main__':

    # Process commandline args.
    files = []
    for arg in sys.argv[1:]:
        if arg in ("-h", "--help"):             # Help
            Usage()
            sys.exit()
        elif arg in ("-d", "--debug"):          # debug output
            debug = True
        elif arg in ("-l", "--list"):           # list rather than convert
            list_only = True
        elif arg in ("-v", "--verbose"):        # verbose output
            verbose = True
        elif arg in ("-s", "--swap"):           # force byte-swap
            force_swap = True
        elif arg in ("-w", "--warnings"):       # print FITS warning msgs
            warn = True
        else:
            files.append(arg)

    if not warn:
        warnings.simplefilter('ignore', category=VerifyWarning)

    # If there were no files specified on the commandline, find all files in
    # and below the current directory.
    if len(files) == 0:
        files = findAll('*.imh', './')

    for f in files:
        base = os.path.basename(f)
        root = os.path.splitext(base)[0]
        extn = os.path.splitext(base)[1]

        if extn == '.imh':
            path = os.path.dirname(f)
            if path == '':
                path = './'
            pixfile = findFile(root+'.pix', path)

            if not os.path.exists(f):
                raise Exception("Header file '%s' does not exist" % f)
            if not os.path.exists(pixfile):
                raise Exception("Pixel file '%s' does not exist" % pixfile)

            # Create an OIF object from the input files.
            image = OIFImage(f, pixfile)

            # Save to a FITS file.
            image.toFITS(root+'.fits', verbose=verbose)
        else:
            continue
