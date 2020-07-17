# IMH2FITS 

## Description

Convert a 1-D or 2-D IRAF OIF image (i.e. a ".imh" file and the associated
".pix" file) to FITS.  This task will read both old V1 files written in a
native byte-order as well as big-endian V2 files and swap bytes as needed.
IMH2FITS can be used to recover old IRAF image data and convert them to 
FITS files when 

    1) no IRAF installation is currently available
    2) old IRAF image files cannot be read because of byte-swapping

Output FITS files are written to the current directory using a filename
derived from the ".imh" image filename.

## Usage:

        % imh2fits [<opts>] <files>

    where <opts> include:
        -d              print debug output
        -h              this help message
        -l              list rather than convert images
        -s              force byte swapping
        -v              print verbose output
        -w              print FITS validation warning messages

        <files>         a list of ".imh" file to convert

    The task will automatically search for a ".pix" pixel file that has the
    same root name as the ".imh" file.  The output FITS file will be written
    to the current directory using the root filename and a ".fits" extension.

## Examples:

    1) Convert all ".imh" files in the current directory to FITS:

        % python imh2fits -v *.imh

       The '-v' flag will print some minimal information for each converted
       file as they are processed.



## Disk resident OIF header layout:

IRAF OIF image files are comprised of a ".imh" header file, and a ".pix"
pixel file.  Each file contains a header record that corresponds to the 
IRAF data structure used in the IMIO kernel interface followed by either
a string containing the "user area" (i.e. header keywords in '.imh' files)
or the image pixels (in a ".pix" file).  The original V1 files were written
in native byte order and used SPP strings (basically short integers) and
are typically now byte-swapped relative to x86-based machines.  The later
V2 files were written in a fixed big-endian format using packed strings
an contain a flag to indicate whether the pixels are swapped relative to the
big-endian format. 

The IMH2FITS task will unpack the header/pixels from either V1 or V2 files
and perform the needed byte-swapping to retrieve the information needed to
construct a valid FITS image.

The OIF format is not documented specifically in terms of a byte layout 
for the file, but rather by how the IMIO interface wrote the data structure
for a particular version.  The file layout is as follows:

```
Field         OIF V1       OIF V2     Description                    Type
-----         ------       ------     -----------                    ----
MAGIC         1-12         1-6        magic string                   int4
HDRLEN        13-16        7-10       length of image header         int4
PIXTYPE       17-20        11-14      datatype of the pixels         int4
SWAPPED                    15-18      pixels are byte-swapped        int4
NDIM          21-23        18-21      number of dimensions           int4
LEN[7]        25-51        22-49      length of the dimensions       int4
PHYSLEN[7]    52-79        50-77      physical length (as stored)    int4
SSMTYPE       80-83        78-81      type of subscript mapping      int4
LUTOFF        84-87        82-85      offset to subscript map luts   int4
PIXOFF        88-91        86-89      offset of the pixels           int4
HGMOFF        92-95        90-93      offset of hgm pixels           int4
BLIST         96-99        94-97      offset of bad pixel list       int4
SZBLIST       100-103      98-101     size of bad pixel list         int4
NBPIX         104-107      102-105    number of bad pixels           int4
CTIME         108-109      106-109    time of image creation         int4
MTIME         112-115      110-113    time of last modify            int4
LIMTIME       116-119      114-117    time min,max computed          int4

MAX           120-123      118-121    max pixel value                float32
MIN           124-127      122-125    min pixel value                float32
HGM           132-208                 (not used)                     float32
CTRAN         209-412                 (not used)                     float32

PIXFILE       413-572      126-382    pixel storage file             char
HDRFILE       573-731      383-638    header storage file            char
TITLE         732-888      639-1022   image name string              char
HISTORY       889-1928     1023-2046  history comment strings        char

UserArea      2053-EOF     2047-EOF   User area (keywords)           char


Header Sizes:
  header      2052 bytes   2046 bytes
  pixfile     1024 bytes   2048 bytes

Notes:
    - field offsets may be word-aligned
    - V1 images are written in their native by order.  This may be swapped
      relative to the current machine (e.g. if on an Intel PC, old Sun Sparc
      images will need swapping)
    - V2 images are written using big-endian byte order (i.e. FITS byte
      order, swapped relative to x86)
    - SPP strings are written as short integers
    - pixfile headers are the same except for:
        - a different magic string
        - no TITLE or HISTORY fields in the header block
```
