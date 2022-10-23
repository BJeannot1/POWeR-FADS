#Chart output of NIDMOWR
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import gridspec

def make_graph(outputfile,zm,zf,zw,time,flux_m,flux_f,zbot_f_):
    #defining subplots
    nbrlignesplot=2
    nbrcolonnesplot=1
    axs = [ [ 0 for i in range(nbrcolonnesplot) ] for j in range(nbrlignesplot) ]
    fig = plt.figure(figsize=(18, 15),dpi=300)
    index=0
    gs = gridspec.GridSpec(10, 1)
    axs[0][0]=plt.subplot(gs[0:7,0])
    axs[1][0]=plt.subplot(gs[7:10,0])       
    #MAIN PLOT
    index=-1
    for nligne in range(0,nbrlignesplot,1):
        for ncolonne in range(0,nbrcolonnesplot,1):
            index=index+1
            #PLOT OF ZM, ZF and ZW
            if (index==0):                
                #zm
                y1 =zm
                x1 =time/3600
                axs[nligne][ncolonne].plot(x1, y1,"--", label = "Zm",color="black",linewidth=10)
                #zf
                y2 =zf
                x2 =time/3600
                axs[nligne][ncolonne].plot(x2, y2,"--", label = "Zf",color=(255/255,9/255,9/255),linewidth=10)
                #zw
                y3 =zw
                x3 =time/3600
                axs[nligne][ncolonne].plot(x3, y3, label = "Zw",color=(91/255,155/255,213/255),linewidth=7)
                #zbot_f
                y4 =[zbot_f_]*len(time)
                x4 =time/3600
                axs[nligne][ncolonne].plot(x4, y4,"--", label = "zbot_f",color="grey",linewidth=1.5)
                textstr =  r'$\mathrm{z_{bot}}%.2f$' +"_"+r'$\mathrm{_{f}}%.2f$'
                axs[nligne][ncolonne].text(time[-1]/3600, zbot_f_, textstr, fontsize=43,verticalalignment='bottom',horizontalalignment='right')                               
      
                #set Y ticks
                axs[nligne][ncolonne].tick_params(axis='y', which='major', labelsize=30,size=15)
                #No x ticks
                axs[nligne][ncolonne].set_xticks([])
                #Name of subplot in angle
                props = dict(boxstyle='square', facecolor="white", alpha=0.5, pad=0.0)
                axs[nligne][ncolonne].text(1.0, 0.0, "a", transform=axs[nligne][ncolonne].transAxes, fontsize=35,
                verticalalignment='bottom',horizontalalignment='right', bbox=props)
                #Title Y axis
                axs[nligne][ncolonne].set_ylabel('Z(m)', fontsize = 35, labelpad=64)
            #PLOT OF FLUXES
            elif (index==1):
                #flux well==>matrix
                y1 =flux_m
                x1 =time/3600
                axs[nligne][ncolonne].plot(x1, y1,label = "flux well==>matrix",color="black",linewidth=7)
                y2 =flux_f
                x2 =time/3600
                #flux fracture==>well
                axs[nligne][ncolonne].plot(x2, y2,label = "flux fracture==>well",color=(255/255,9/255,9/255),linewidth=7)
                axs[nligne][ncolonne].tick_params(axis='y', which='major', labelsize=30,size=15)            
                #Name of subplot in angle
                props = dict(boxstyle='square', facecolor="white", alpha=0.5, pad=0.0)
                axs[nligne][ncolonne].text(1.0, 0.0, "b", transform=axs[nligne][ncolonne].transAxes, fontsize=35,
                verticalalignment='bottom',horizontalalignment='right', bbox=props)
                #Title Y axis
                axs[nligne][ncolonne].set_ylabel('Flux (m'+chr(0x00B3)+'/h)', fontsize = 35, labelpad=35)
                #Title X axis
                fig.supxlabel('Time (h)',y=0.025,size=35)
                #X ticks        
                axs[nligne][ncolonne].tick_params(axis='x', which='major', labelsize=30,size=15)
                axs[nligne][ncolonne].xaxis.set_major_locator(ticker.MultipleLocator(15))
                axs[nligne][ncolonne].tick_params(axis='x', which='minor',size=5)

    #No space between subplots
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.0, hspace=0)
    #save subplots
    plt.savefig(outputfile,bbox_inches='tight')
    plt.close()
