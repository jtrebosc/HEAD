import re
import numpy as np
import sys


def loadSimpsonFile(filePath):
    """
    Loads SIMPSON file. Both ASCII and binary data are supported. As well
    as 1D and 2D data.

    Parameters
    ----------
    filePath: string
        Path to the file that should be loaded

    Returns
    -------
    SpectrumClass
        SpectrumClass object of the loaded data
    """
    with open(filePath, 'r') as f:
        Lines = f.read().split('\n')
    NI, SW, SW1, TYPE, FORMAT, REF = 1, 0, 0, '', 'Normal', 0
    DataStart = Lines.index('DATA')
    DataEnd = Lines.index('END')
    for s in range(0, DataStart):
        if Lines[s].startswith('NI='):
            NI = int(re.sub('NI=', '', Lines[s]))
        elif Lines[s].startswith('SW='):
            SW = float(re.sub('SW=', '', Lines[s]))
        elif Lines[s].startswith('SW1='):
            SW1 = float(re.sub('SW1=', '', Lines[s]))
        elif Lines[s].startswith('TYPE='):
            TYPE = re.sub('TYPE=', '', Lines[s])
        elif Lines[s].startswith('FORMAT='):
            FORMAT = re.sub('FORMAT=', '', Lines[s])
        elif Lines[s].startswith('REF='):
            REF = re.sub('REF=', '', Lines[s])
    if 'Normal' in FORMAT:
        length = DataEnd - DataStart - 1
        data = np.zeros(length, dtype=complex)
        for i in range(length):
            temp = Lines[DataStart + 1 + i].split()
            data[i] = float(temp[0]) + 1j * float(temp[1])
    elif 'BINARY' in FORMAT:
        # Binary code based on:
        # pysimpson: Python module for reading SIMPSON files
        # By: Jonathan J. Helmus (jjhelmus@gmail.com)
        # Version: 0.1 (2012-04-13)
        # License: GPL
        chardata = ''.join(Lines[DataStart + 1:DataEnd])
        nquads, mod = divmod(len(chardata), 4)
        assert mod == 0     # character should be in blocks of 4
        BASE = 33
        charst = np.frombuffer(chardata.encode('ascii'), dtype=np.uint8)
        charst = charst.reshape(nquads, 4) - BASE

        def FIRST(f, x):
            # This function takes a np.uint8 ndarray and returns it with all bits more than f places from the right set to 0.
            # For example, 102 is 01100110 in binary, and 00100110, or 38, is returned if f is 6
            return x & np.uint8((1 << f) - 1)

        def LAST(f, x):
            # This function takes a np.uint8 ndarray and returns it with all bits more than f places from the left set to 0.
            # For example, 102 is 01100110 in binary, and 01100100, or 100, is returned if f is 6
            return x & np.uint8(256 - (1 << 8 - f))

        first = FIRST(6, charst[:, 0]) | LAST(2, charst[:, 1] << 2)
        second = FIRST(4, charst[:, 1]) | LAST(4, charst[:, 2] << 2)
        third = FIRST(2, charst[:, 2]) | LAST(6, charst[:, 3] << 2)
        Bytes = np.ravel(np.transpose(np.array([first, second, third]))).astype('int64')
        # convert every 4 'bytes' to a float
        num_points, num_pad = divmod(len(Bytes), 4)
        Bytes = np.array(Bytes)
        Bytes = Bytes[:-num_pad]
        Bytes = Bytes.reshape(num_points, 4)
        mantissa = ((Bytes[:, 2] % 128) << 16) + (Bytes[:, 1] << 8) + Bytes[:, 0]
        exponent = (Bytes[:, 3] % 128) * 2 + (Bytes[:, 2] >= 128) * 1
        negative = Bytes[:, 3] >= 128
        e = exponent - 127
        m = np.abs(mantissa) / np.float64(1 << 23)
        data = np.float32((-1)**negative * np.ldexp(m, e))
        data = data.view('complex64')
    if NI != 1:  # 2D data, reshape to NI, NP
        data = data.reshape(int(NI), -1)

    return (data, NI, SW, SW1, REF, TYPE)


# read SPE file
# calculate LEFT and RIGHT
# Write file with format
# F1LEFT = 18.2412229965813 ppm. F1RIGHT = -6.751947497302158 ppm.
# F2LEFT = 255.676342688474 ppm. F2RIGHT = -244.1870671891952 ppm.
#
# NROWS = 512 ( = number of points along the F1 axis)
# NCOLS = 8192 ( = number of points along the F2 axis)
#
# followed by real values
(data, ni, sw, sw1, ref, kind) = loadSimpsonFile(sys.argv[1])
npts = data.shape[-1]
x = np.arange(npts, dtype=np.float64)/npts*sw - sw/2
y = np.arange(ni, dtype=np.float64)/ni*sw1 - sw1/2
left = x[-1]
right = x[0]
left1 = y[-1]
right1 = y[0]

print(f"# F1LEFT = {left1} Hz. F1RIGHT = {right1} Hz.")
print(f"# F2LEFT = {left} Hz. F2RIGHT = {right} Hz.")
for j in range(ni):
    for i in range(npts):
        print(data[ni-1-j][npts-1-i].real)
