import numpy as np
import matplotlib.pyplot as plt
import argparse, time, matplotlib, weakref
import quantum_matched_filter_functions as qmffn
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

np.random.seed(int(time.time()))



def main(N, P, t, outpath, fontsize=28, ticksize=22, figsize=(10,10)):

    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'
    
    fig = plt.figure(figsize=figsize)
    ax = fig.gca()
    ax = [ax]
    
    lw=1.5
    fontsize=16
    
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['left'].set_visible(False)
    ax[0].spines['bottom'].set_visible(False)
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[0].set_xlim(-1,1)
    ax[0].set_ylim(-1,1)

    ax[0].annotate('', xy=(-1,0), xytext=(1,0), fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw))
    ax[0].annotate('', xy=(0,-1), xytext=(0,1), fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw))
    
    ax[0].text(-.07, 0.5, r'$|\rho<\rho_{th}\rangle$', fontsize=fontsize, rotation=90)
    ax[0].text(0.5, -0.07, r'$|\rho\geq\rho_{th}\rangle$', fontsize=fontsize)
    #ax[0].text(-.07, 0.5, r'$|T_{bad}\rangle$', fontsize=fontsize, rotation=90)
    #ax[0].text(0.5, -0.07, r'$|T_{good}\rangle$', fontsize=fontsize)

    theta = np.arcsin(np.sqrt(float(t)/N))
    state = [np.array([np.sqrt(float(t)/N),np.sqrt((float(N)-t)/N)])]

    ax[0].annotate(r'$|\Psi_{0} \rangle$', xy=(0,0), xytext=state[-1], fontsize=fontsize, arrowprops=dict(arrowstyle="-", lw=lw, ls=':', color='black'))
    ax[0].annotate('', xy=(0,0), xytext=-state[-1], fontsize=fontsize, arrowprops=dict(arrowstyle="-", lw=lw, ls=':', color='black'))

    lines_f = []
    lines_sf = []
    lines = []

    def init():
        lines.append(ax[0].annotate(r'$|\Psi_{0} \rangle$', xy=(0,0), xytext=state[-1], fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw, color='red')))
        lines_sf.append(ax[0].annotate('', xy=(0,0), xytext=(0,0), fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw, color='white')))
        lines_f.append(ax[0].annotate('', xy=(0,0), xytext=(0,0), fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw, color='white')))

        return fig,

    def animate(k):
        if k==0:
            return fig,
        if k%2!=0:
            lines_f[-1].remove()
            #lines_sf[-1].remove()
            lines[-1].remove()
            lines_sf.append(ax[0].annotate(r'$|\Psi_{kscore} \rangle$'.format(kscore=int(k/2)), xy=(0,0), xytext=state[-1], fontsize=fontsize, arrowprops=dict(arrowstyle="-", lw=lw, color='red', ls='--')))
            lines_f.append(ax[0].annotate(r'$|\Psi^{\ast}$'+r'$_{kscore} \rangle$'.format(kscore=int(k/2)), xy=(0,0), xytext=state[-1]* np.array([-1,1]), fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw, color='green')))
        else:
            p = k/2
            lines_sf[-1].remove()
            state.append(np.array([np.sin(((2.*p)+1)*theta),np.cos(((2.*p) +1)*theta)]))
            lines.append(ax[0].annotate(r'$|\Psi_{kscore} \rangle$'.format(kscore=int(p)), xy=(0,0), xytext=state[-1], fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw, color='red')))
            print(p)

        return fig,

    frames = (P*2)
    im_ani = matplotlib.animation.FuncAnimation(fig, animate, init_func=init, frames=frames, interval=100, blit=True)

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=.2, metadata=dict(artist='Me'), bitrate=1800)
    im_ani.save(outpath+'psi1_state_short.mp4', writer=writer)#outpath+str(P)+'_'+str(N)+'_'+str(t)+'_psi1_state_ani.mp4', writer=writer)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--N', help="", type=int, required=True)
    parser.add_argument('--P', help="", type=int, required=True)
    parser.add_argument('--t', help="", type=int, required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)

    opt = parser.parse_args()
 
    main(opt.N, opt.P, opt.t, opt.outpath)
