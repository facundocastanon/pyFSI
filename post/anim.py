
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter

from pyFSI.post.plots import *
import numpy as np

class pAnimation():
    def __init__(self,
                 P,
                 tfilePath,
                 name="animation",
                 interval=0,
                 frames='all',
                 xlim=False,
                 ylim=False,
                 mode=None,
                 startFrame=0,
                 endFrame=-1):

        self.P = P # The plot object
        self.mode = mode
        self.interval = interval
        if frames == 'all':
            self.frames = np.arange(0, len(self.P.yData[:,0]), 1)
        else:
            self.time = np.loadtxt(tfilePath)
            size = len(self.time[startFrame:endFrame])
            skip = round(size / frames)
            self.frames = np.arange(startFrame, size+startFrame, skip)
        self.animate()

    def animate(self):
        self.ani = animation.FuncAnimation( self.P.fig, 
                                            self.update,
                                            frames=self.frames)

    def update(self, frame):
        self.P.update(frame)

    def save(self, name, fps=30):
        writer = animation.FFMpegWriter(fps=fps,
                                        metadata=dict(artist='Me'),
                                        bitrate=1800)
        self.ani.save(name + '.mp4', writer=writer)



