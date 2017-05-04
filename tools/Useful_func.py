def update_print(phrase, appendix=' ', percent=False):
    import sys
    import time
    from time import sleep
    sys.stdout.write('\r')
    if percent==True:
        sys.stdout.write(phrase + '%d%%' % appendix )
        if appendix == 100:
            print '\n'
    if percent==False:
        sys.stdout.write(phrase + appendix )
        if appendix == 'done':
            print '\n'
    sys.stdout.flush()
    sleep(1.e-6)

def hist(values,sizebins):
	import numpy as np
        mm = min(values)
        mx = max(values)
        bins=np.arange(mm,mx,sizebins)
        hist, b_edges = np.histogram(values,bins)
        centers = (b_edges[:-1] + b_edges[1:]) / 2.0
        hist = np.array(hist, dtype=np.float64)
        return centers, hist


