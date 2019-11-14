#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 ,-*
(_) Created on <Sa Okt 19 2019> @ 18:44:15

@authors: Boris Daszuta
@function:

"""
import time

import utils.asset_tools as uat

import numpy as np
import os

import matplotlib.pyplot as plt
import matplotlib.animation as plta

###############################################################################
ix_max = 321
###############################################################################

dir = os.path.join(uat._DIR_ATHENA, 'outputs', 'wave', 'wave_1d_cvg_Gaussian_P')
fn = lambda ix: os.path.join(dir, 'wave_1d.out1.{ix}.athdf'.format(
    ix=str(ix).zfill(5)))


def pl():
    fig = plt.figure()
    ax = plt.axes(xlim=(-4.0, 4.0), ylim=(-0.1, 1.1))
    # ax = plt.semilogy(xlim=(-4.0, 4.0), ylim=(10e-16, 10))
    line, = ax.plot([], [], lw=1)
    data = {}
    for ix in range(0, ix_max):
        tmp = uat._ar.athdf(fn(ix))
        x = tmp['x1v'].flatten()
        y = np.abs(tmp['wU'].flatten())
        data[ix] = [x, y]


    def init():
        line.set_data([], [])
        return line,

    def animate(ix):
        x, y = data[ix]
        line.set_data(x, y)
        return line,

    anim = plta.FuncAnimation(fig, animate, init_func=init,
                              frames=ix_max, interval=20, blit=True)

    # anim.save('sine_wave.gif', writer='imagemagick')
    # for ix in range(ix_max):
    #     data = uat._ar.athdf(fn(ix))
    #     plot = plt.plot(data['x1v'].flatten(), data['wU'].flatten())
    #     fig.canvas.draw()
    #     # time.sleep(0.2)



def ff(ix=0):
    leg_lev = ['ok', 'sb', 'dg', '.r']
    tmp = uat._ar.athdf(fn(ix), raw=True, return_levels=True)
    tmp_f = uat._ar.athdf(fn(ix))

    plt.plot(tmp_f['x1v'].flatten(), tmp_f['wU'].flatten(), '-b')
    plt.plot(tmp_f['x1v'].flatten(), tmp_f['wError'].flatten(), '-c')
    for lev, x in zip(tmp['Levels'], tmp['x1f']):
        plt.plot(x, 0 * np.ones_like(x), leg_lev[lev], markersize=3)

    plt.show()

#
# :D
#
