import numpy as np
import matplotlib.pyplot as plt
import argparse, time, matplotlib, weakref
import quantum_matched_filter_functions as qmffn
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

np.random.seed(int(time.time()))



def main(N, P, t, outpath, fontsize=28, ticksize=22, figsize=(10,15)):

    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'
    
    fig = plt.figure(figsize=figsize)
    gs = matplotlib.gridspec.GridSpec(3, 2, figure=fig)
    ax = []
    ax.append(fig.add_subplot(gs[:2,:]))
    ax.append(fig.add_subplot(gs[2,:]))

    #ax = fig.gca()
    #ax = [ax]
    
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

    theta = np.arcsin(np.sqrt(float(t)/N))
    state = [np.array([np.sqrt(float(t)/N),np.sqrt((float(N)-t)/N)])]

    #ax[0].annotate(r'$|\Psi_{0}\rangle$', xy=(0,0), xytext=state[-1], fontsize=fontsize, arrowprops=dict(arrowstyle="-", lw=lw, ls=':', color='black'))
    #ax[0].annotate('', xy=(0,0), xytext=-state[-1], fontsize=fontsize, arrowprops=dict(arrowstyle="-", lw=lw, ls=':', color='black'))

    ax[1].plot(np.arange(P), np.sin(((2.*np.arange(P)) + 1.)*theta), color='blue', label=r'$\rho\geq\rho_{th}$')
    ax[1].plot(np.arange(P), np.cos(((2.*np.arange(P)) + 1.)*theta), color='purple', label=r'$\rho<\rho_{th}$')

    ax[1].set_xlabel('Grover operations', fontsize=fontsize)
    ax[1].set_ylabel('Amplitude', fontsize=fontsize)

    ax[1].set_xticks(np.arange(P).astype(int))

    lines_f = []
    lines_sf = []
    lines = []
    lines2 = []
    lines_f2 = []
    line1 = []
    line2 = []
    linesx = []
    linesy = []

    def init():
        lines.append(ax[0].annotate('', xy=(0,0), xytext=state[-1], fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw, color='red')))
        #lines.append(ax[0].annotate(r'$|\Psi_{0}\rangle$', xy=(0,0), xytext=state[-1], fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw, color='black')))
        line1.append(ax[1].scatter(0, state[-1][0], color='blue', lw=2, marker='o'))
        line2.append(ax[1].scatter(0, state[-1][1], color='purple', lw=2, marker='o'))
        return fig,

    def animate(k):
        if k==0:
            #state.append(np.array([np.sqrt(float(t)/N),np.sqrt((float(N)-t)/N)]))
            #lines2.append(ax[0].annotate('', xy=(0,0), xytext=(0,0), fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw, color='white')))
            #lines_sf.append(ax[0].annotate('', xy=(0,0), xytext=(0,0), fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw, color='white')))
            lines_f.append(ax[0].annotate('', xy=(0,0), xytext=(0,0), fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw, color='white')))
            linesx.append(ax[0].annotate('', xy=(state[-1][0],0), xytext=state[-1], fontsize=fontsize, arrowprops=dict(arrowstyle="-", lw=lw, ls=':', color='blue')))
            linesy.append(ax[0].annotate('', xy=(0,state[-1][1]), xytext=state[-1], fontsize=fontsize, arrowprops=dict(arrowstyle="-", lw=lw, ls=':', color='purple')))
            return fig,
        if k%2!=0:
            return fig,
            #lines_f[-1].remove()
            #lines_sf[-1].remove()
            #lines[-1].remove()
            #lines_sf.append(ax[0].annotate(r'$|\Psi_{kscore} \rangle$'.format(kscore=int(k/2)), xy=(0,0), xytext=state[-1], fontsize=fontsize, arrowprops=dict(arrowstyle="-", lw=lw, color='red', ls='--')))
            #print(np.arctan(state[-1][0]/state[-1][1]))
            #state.append(state[-1] * np.array([-1,1]))
            #print(np.arctan(state[-1][0]/state[-1][1]))
            #print(state[-1],state_f)
            #lines_f.append(ax[0].annotate(r'$\Psi^{\ast}$'+r'$_{kscore} \rangle$'.format(kscore=int(k/2)), xy=(0,0), xytext=state[-1]* np.array([-1,1]), fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw, color='green')))
        else:
            p = k/2
            #p_ = (k/2)-1
            lines[-1].remove()
            linesx[-1].remove()
            linesy[-1].remove()
            #lines2[-1].remove()
            #lines_sf[-1].remove()
            #theta_ = np.arctan(state[-1][0]/state[-1][1])
            #print(theta_,theta)
            #state.append(np.matmul(np.array([[np.cos(4*theta_),-np.sin(4*theta_)],[np.sin(4*theta_),np.cos(4*theta_)]]),state[-1]))
            #print(state[-1])
            #print([np.cos(((2.*p)+1)*theta),phase*np.sin(((2.*p) +1)*theta)])
            #print([np.cos(((2.*(p-1))+1)*theta),phase*np.sin(((2.*(p-1)) +1)*theta)])
            #print([np.cos(((2.*(p+1))+1)*theta),phase*np.sin(((2.*(p+1)) +1)*theta)])
            #state.append(np.array([np.sin(((2.*p)+1)*theta),phase*np.cos(((2.*p) +1)*theta)]))#,np.cos(((2.*k)+1)*theta)]))
            #lines.append(ax[0].annotate(r'$|\Psi_{kscore} \rangle$'.format(kscore=int(p)), xy=(0,0), xytext=state[-1], fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw, color='red')))
            state.append(np.array([np.sin(((2.*p)+1)*theta),np.cos(((2.*p) +1)*theta)]))
            lines.append(ax[0].annotate(''.format(kscore=int(p)), xy=(0,0), xytext=state[-1], fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw, color='red')))
            #lines.append(ax[0].annotate(r'$|\Psi_{kscore} \rangle$'.format(kscore=int(p)), xy=(0,0), xytext=state[-1], fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw, color='red')))

            #lines2.append(ax[0].annotate(r'$|\Psi_{kscore} \rangle$'.format(kscore=int(p)), xy=(0,0), xytext=state[-1]*np.array([1,-1]), fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw, ls=':', color='red')))
            #lines_f2.append(ax[0].annotate(r'$\Psi^{\ast}$'+r'$_{kscore} \rangle$'.format(kscore=int(k/2)), xy=(0,0), xytext=state[-2]*np.array([-1,1]), fontsize=fontsize, arrowprops=dict(arrowstyle="<-", lw=lw, ls='--', color='green')))
            #lines_sf.append(ax[0].annotate(r'$|\Psi_{kscore} \rangle$'.format(kscore=int(p_)), xy=(0,0), xytext=state[-3], fontsize=fontsize, arrowprops=dict(arrowstyle="-", lw=lw, color='red', ls='--')))
            line1[-1].remove()
            line2[-1].remove()
            print(p)
            line1.append(ax[1].scatter(p, np.sin(((2.*p)+1)*theta), vmin=-1, vmax=1, color='blue', lw=2, marker='o'))
            line2.append(ax[1].scatter(p, np.cos(((2.*p)+1)*theta), vmin=-1, vmax=1, color='purple', lw=2, marker='o'))
            linesx.append(ax[0].annotate('', xy=(state[-1][0],0), xytext=state[-1], fontsize=fontsize, arrowprops=dict(arrowstyle="-", lw=lw, ls=':', color='purple')))
            linesy.append(ax[0].annotate('', xy=(0,state[-1][1]), xytext=state[-1], fontsize=fontsize, arrowprops=dict(arrowstyle="-", lw=lw, ls=':', color='blue')))

        return fig,

    leg = ax[1].legend(fontsize=fontsize)

    frames = (P*2)
    im_ani = matplotlib.animation.FuncAnimation(fig, animate, init_func=init, frames=frames, interval=100, blit=True)

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=1, metadata=dict(artist='Me'), bitrate=1800)
    im_ani.save(outpath+'psi1_state_ani_with_amp.mp4', writer=writer)#outpath+str(P)+'_'+str(N)+'_'+str(t)+'_psi1_state_ani.mp4', writer=writer)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='', description="Perform one QMF on GW150914 given a template bank")
    parser.add_argument('--N', help="", type=int, required=True)
    parser.add_argument('--P', help="", type=int, required=True)
    parser.add_argument('--t', help="", type=int, required=True)
    parser.add_argument('--outpath', help="", type=str, required=True)

    opt = parser.parse_args()
 
    main(opt.N, opt.P, opt.t, opt.outpath)
