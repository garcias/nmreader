import numpy
from scipy.fftpack import fft, fftshift
import altair
import pandas
import dataclasses

@dataclasses.dataclass
class Parameters(object):
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
class Spectrum(object):
    parameters: Parameters
    fid: numpy.ndarray
    time: numpy.ndarray
    fft: numpy.ndarray
    frequency: numpy.ndarray
    shift: numpy.ndarray

    def __init__( self, filename ):
        # filename is the name of an NMReady .dx file containing the FID
        # returns fid (complex array), time (array), p (dictionary)

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
        # offset represents fraction (0-1) of x-axis range 
        factor = 1.0j * numpy.pi / 180
        index = numpy.linspace(0 - offset, 1 - offset, self.parameters.npoints)
        return self.fft * numpy.exp(factor * (ph0 + ph1*index ))

def parse_intervals( interval_string ):

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
