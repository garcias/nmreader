"""Parse a NMReady JCAMP-DX file to load FID of an experiment and compute FFT."

Provides a Spectrum class to represent experimental data from an NMReady spectrometer,
including the FID, FFT, related axes, and experimental parameters. When you create a Spectrum
object, it will parse a JCAMP file (usually with .dx extension to extract JCAMP headers into a
Parameters object and the FID into a complex array, and store these as its attributes. It will 
also construct and store the time axis, and then compute the FFT with frequency and chemical
shift axes. The method phase() provides a phase-corrected version of the fft.

Also provides utility functions to help with processing the spectrum. They will likely need to
work in tandem with an interactive GUI, to obtain meaningful arguments. 
    parse_intervals() takes a string denoting intervals that encompass signals
    fit_baseline() takes spectral axes and intervals to fit a polynomial to the baseline

Example:
    file = 'test.dx'
    spec = nmreader.Spectrum( file )
    phased = spec.phase( 1490, -10, -25 )  # pivot the first-order correction at index 1490
    intervals = nmreader.parse_intervals("solvent 7.3 0.1")
    phased_fitted = nmreader.fit_baseline( spec.shift, phased.real, intervals )
"""

import numpy
from scipy.fftpack import fft, fftshift
import altair
import pandas
import dataclasses

@dataclasses.dataclass
class Parameters:
    """Class to represent experimental parameters of NMReady NMR experiment.

    Most attributes are read from the JCAMP ## headers, but detalt is inferred (in accordance
    with the JCAMP spec) from firstx and lastx. In order to construct meaningful arrays to 
    represent the fid and fft, three parameters (npoints, deltat, observe_frequency) are required.

    Attributes:
        npoints (int) 
        deltat (float): inferred from (final_time - first_time) / (npoints - 1)
        observe_frequency (float): spectrometer frequency
        title (str): name of file recorded during acquisition
        sample (str): notes entered by user
        temperature (float): in Centigrade units
        scans (int): number of scans averaged
        date (str): e.g. 2016/12/09 07:54:21-0700
        nucleus (str): e.g. ^1H
        solvent (str): e.g. Chloroform-d
    """
    npoints : int
    deltat : float
    observe_frequency : float
    title : str
    temperature : float
    date : str
    sample : str
    nucleus : str
    solvent : str
    scans : int

@dataclasses.dataclass( init=False )
class Spectrum:
    """Class to represent experimental data from NMReady NMR experiment.

    Spectrum requires a path to the dx file, it will parse the file to extract the
    JCAMP headers and fid signal. It stores the header data in a Parameters object,
    and computes the fft and stores it as an array. The method phase() allows for
    linear phase correction up to 1st order, with an adjustable pivot.

    Attributes:
        parameters (Parameters): extracted or inferred from JCAMP headers
        fid (numpy.ndarray): free-induction decay signal, represented as array of complex numbers
        time (numpy.ndarray): time axis for fid, inferred from delta t and npoints
        fft (numpy.ndarray): 
        frequency (numpy.ndarray): frequency axis for fft, inferred from deltat
        shift (numpy.ndarray): frequency axis converted to chemical shift units, inferred 
            from deltat and observe_frequency

    Methods:
        __init__( filename ): Spectrum requires path to a dx file containing the data.
        phase( offset, ph0, ph1 ): Applies phase correction and returns FFT
    """

    parameters: Parameters
    fid: numpy.ndarray
    time: numpy.ndarray
    fft: numpy.ndarray
    frequency: numpy.ndarray
    shift: numpy.ndarray

    def __init__( self, filename ):
        """Parses file, sets parameters, computes and sets fid and fft and raw axes.

        Args:
            filename (str): local path to NMReady dx file. File must contain the headers -- 
                FIRST, LAST, NPOINTS, .OBSERVE FREQUENCY -- in order to construct the fid and fft.
        """

        # Don't split into lines yet, so we can use string methods to detect blocks
        with open(filename, 'r') as file:
            f = file.read()

        # split into blocks for headers, real data, and imaginary data
        header_block, data_block = f.split( '##DATA TABLE= (X++(R..R)), XYDATA', maxsplit=1 )
        real_block, imag_block = data_block.split( '##DATA TABLE= (X++(I..I)), XYDATA', maxsplit=1 )

        # parse parameters from headers block
        p = dict( [
            line.split("##")[1].split('=', maxsplit = 1) for line in header_block.split('\n')
            if line.startswith('##')
        ] )

        first_time = float( p['FIRST'].split(',')[0] )
        final_time = float( p['LAST'].split(',')[0] )
        npoints = float( p['NPOINTS'] )
        deltat = (final_time - first_time) / (npoints - 1)
        sampling_frequency = 1 / deltat
        observe_frequency = float(p['.OBSERVE FREQUENCY'])
        real_factor = float( p['FACTOR'].split(',')[1] )
        imag_factor = float( p['FACTOR'].split(',')[2] )

        parameters = {
            'npoints' : int(npoints), 'deltat' : float(deltat), 
            'observe_frequency' : float(observe_frequency),
            'title' : p['TITLE'], 'temperature' : float(p['TEMPERATURE']),
            'date' : p['LONG DATE'], 'sample' : p['SAMPLE DESCRIPTION'],
            'nucleus' : p['.OBSERVE NUCLEUS'], 'solvent' : p['.SOLVENT NAME'],
            'scans' : int(p['.AVERAGES']), 
        }

        self.parameters = Parameters( **parameters )

        # parse data blocks into arrays, then combine into complex array
        real_table = [
            line.split()[1:] for line in real_block.split('\n')
            if len( line.split()[1:] ) == 4
        ]

        imag_table = [
            line.split()[1:] for line in imag_block.split('\n')
            if len( line.split()[1:] ) == 4
        ]

        real_list = [ i for row in real_table for i in row ]
        imag_list = [ i for row in imag_table for i in row ]

        real_array = numpy.array( real_list ).astype(float) * real_factor 
        imag_array = numpy.array( imag_list ).astype(float) * imag_factor

        self.fid = real_array + 1.0j * imag_array
        self.time = numpy.linspace( first_time, final_time, self.parameters.npoints )

        # returns arrays for spec (complex), frequency (Hz), shift (ppm)
        self.fft = fftshift( fft( self.fid ) )
        self.frequency = numpy.linspace( 0, 1/deltat, self.parameters.npoints )
        self.shift = self.frequency / observe_frequency

        return None

    def phase( self, offset, ph0, ph1 ):
        """Applies phase correction and returns FFT.

        Args:
            offset (float): "pivot" position for first-order correction, as a fraction in [0,1] 
                of the entire frequency axis. The first-order correction applied will be a 
                linear function of frequency centered on this offset.
            ph0 (float): Zeroth order phase coefficient (in degrees)
            ph1 (float): First order phase coefficient (in degrees), centered on offset

        Returns:
            numpy.ndarray: complex array with same length (npoints) as fft 
        """

        factor = 1.0j * numpy.pi / 180
        index = numpy.linspace(0 - offset, 1 - offset, self.parameters.npoints)
        return self.fft * numpy.exp(factor * (ph0 + ph1*index ))

def parse_intervals( interval_string ):
    """Parses a string that specifies intervals on a spectral axis.

    Args:
        interval_string (str): Each interval specified by three attributes: label, center, width.
            Attributes separated by white space. Multiple intervals separated by newline.
    
    Returns:
        list: each item is a dictionary with keys: name, center, width
    """
    interval_list = [ interval.split() for interval in interval_string.split("\n") ]
    intervals = [
        { 'name' : interval[0], 'center' : interval[1], 'width' : interval[2] } 
        for interval in interval_list if len(interval) > 2
    ]

    df = pandas.DataFrame( columns = ['name', 'center', 'width'] )
    df = df.append(intervals)
    df['center'] = pandas.to_numeric( df['center'], errors='coerce').astype(float)
    df['width']  = pandas.to_numeric( df['width'],  errors='coerce').astype(float)
    df['left']   = df['center'] + df['width']/2
    df['right']  = df['center'] - df['width']/2

    return df.dropna().to_dict(orient='records')

def fit_baseline( shift, real, intervals, order=4 ):
    """Computes a polynomial fit to spectrum, with non-baseline intervals masked.
    
    Args:
        shift (numpy.ndarray): chemical shift axis of spectrum
        real (numpy.ndarray): real component of fft (should be phased for best results)
        intervals (list of dict): intervals specifying center frequency and width, should
            be the valid result of function `parse_intervals`
        order (int): order of the polynomial, defaults to 4
    
    Returns:
        numpy.ndarray: the best-fit baseline, in same units as shift array
    """
    df = pandas.DataFrame({'shift': shift, 'real': real})
    df['masked'] = False
    
    mask_indices = [ 
        df.query(f"{interval['right']} < shift < {interval['left']}").index
        for interval in intervals 
    ]

    for index in mask_indices:
        df.loc[index, 'masked'] = True

    ppm_array = numpy.ma.masked_array( df['shift'], df['masked'] )
    coeffs = numpy.ma.polyfit( ppm_array, df['real'], order )
    return numpy.polyval( coeffs, df['shift'] )
