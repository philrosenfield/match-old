import os

if 'Linux' in os.uname():
    import matplotlib as mpl
    mpl.use('Agg')
