"""
Header Class
===========


This class creates the header that will be written in the images created by
the AIS.
"""

import astropy.io.fits as fits


class Header:
    """Images header class.

    Parameters
    ----------
    dic : dictionary
        Parameter of the operation mode of the SPARC4 camera

    Returns
    -------
    None.
    """

    def __init__(self, ccd_operation_mode):
        """Initialize the Header Class."""
        self.em_mode = ccd_operation_mode['em_mode']
        self.noise_factor = 1
        self.em_gain = 1
        if self.em_mode == 1:
            self.noise_factor = 1.41
            self.em_gain = ccd_operation_mode['em_gain']
        self.preamp = ccd_operation_mode['preamp']
        self.hss = ccd_operation_mode['hss']
        self.bin = ccd_operation_mode['bin']
        self.t_exp = ccd_operation_mode['t_exp']

    def _create_header(self):
        """Create the image header.

        This functions writes a simple header with the used parameters for
        the CCD operation mode to the image cube file
        """
        hdr = fits.Header()
        hdr['NAXIS1'] = (200, 'length of data axis 1')
        hdr['NAXIS2'] = (200, 'length of data axis 2')
        hdr['EXTEND'] = ('T', 'FITS dataset may contain extensions')
        hdr['COMMENT'] = 'and Astrophysics, volume 376, page 359; bibcode:' + \
            '2001A&A...376..3'
        hdr['ACQMODE'] = ('Single  ', 'Acquisition Mode')
        hdr['READMODE'] = ('Image   ', 'Readout Mode')
        hdr['IMGRECT'] = ('1, 200,200, 1', 'Image Format')
        hdr['HBIN'] = (self.bin, 'Horizontal Binning')
        hdr['VBIN'] = (self.bin, 'Vertical Binning')
        hdr['TRIGGER'] = ('Internal', 'Trigger Mode')
        hdr['EXPOSURE'] = (self.t_exp, 'Total Exposure Time')
        hdr['TEMP'] = (self.ccd_temp, 'Temperature')
        hdr['READTIME'] = (str(1/self.hss)+'E-006', 'Pixel readout time ')
        hdr['VSHIFT'] = ('4.33E-06', 'Vertical Shift Speed')
        hdr['GAIN'] = (self.gain, 'Preamp Gain (e-/ADU)')
        em_mode = 'Conventional'
        if self.em_mode == 1:
            em_mode = 'Electron Multiplying'
        hdr['OUTPTAMP'] = (em_mode, 'Output Amplifier')
        hdr['EMGAIN'] = (self.em_gain, 'Electron Multiplying Gain')
        hdr['PREAMP'] = (str(self.preamp)+'x', 'Pre Amplifier Gain')
        hdr['SERNO'] = (self.serial_number, 'Serial Number')
        hdr['DATE'] = ('2017-07-14T00:00:58',
                       'File Creation Date (YYYY-MM-HHThh:mm:ss)')
        hdr['FRAME'] = ('2017-07-14T00:00:58.642', 'Start of Frame Exposure')
        hdr['IMAGE'] = ('hats-24_I_transito_001', 'Nome do arquivo')
        self.hdr = hdr

    def get_header(self):
        """Return the image header.

        Returns
        -------
        Image header: dictionary
            Information related to the artificial image created by the AIS
        """
        return self.hdr
