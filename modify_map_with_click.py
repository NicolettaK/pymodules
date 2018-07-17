from matplotlib import pyplot as plt
import healpy as hp
import numpy as np

class LineBuilder:
    def __init__(self, line):
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        print('click', event)
        if event.inaxes!=self.line.axes: return
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()
        coords.append((event.xdata, event.ydata))

def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata
    print 'press_enter'
    coords.append((ix, iy))
    fig.canvas.mpl_disconnect(cid)
    return coords

def transform_coords(coords):
    MP = hp.projector.MollweideProj()
    vectors = []
    for p in range(len(coords)):
        vec = MP.xy2vec(coords[p][0], coords[p][1])
        vectors.append(vec)
    vectors = np.array(vectors)
    return vectors

def modify_map(hpmap, value_to_fill):
    global fig, cid, coords
    coords = []
    fig = plt.figure()
    hp.mollview(hpmap, fig=fig.number)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    raw_input()
    print coords
    line, = plt.plot(coords[0][0], coords[0][1])  # empty line
    linebuilder = LineBuilder(line)
    raw_input()
    coords = np.array(coords)
    coords = np.delete(coords, -1, 0)
    vec = transform_coords(coords)
    nside = hp.npix2nside(len(hpmap))
    to_mask = hp.query_polygon(nside, vec)
    map_out = np.copy(hpmap)
    map_out[to_mask] = value_to_fill
    return map_out
