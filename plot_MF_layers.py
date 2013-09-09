# functions to visualize layer structure in a MODFLOW model
# 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from collections import defaultdict
import pdb

# input file
infile="plot_MF_layers.in"

# Get input files locations
infile=open(infile,'r').readlines()
inputs=defaultdict()
for line in infile:
    if "#" not in line.split()[0]:
        varname=line.split("=")[0]
        var=line.split("=")[1].split()[0]
        inputs[varname]=var.replace("'","") # strip extra quotes

modelbottoms=inputs["modelbottoms"]
L1top=inputs["L1top"]
xsname=inputs["xsname"]
colors=inputs["colors"].strip('[]').split(',')
interval=int(inputs["interval"])


# get dimmensions
temp=open(L1top).readlines()
cols=len(temp[0].strip().split())
rows=len(temp)

def set_rowscols(interval,nrows_cols,inputsname):
    try:
        start,end=map(int,inputs[inputsname].strip('[]').split(','))
        start,end=start-1,end
    except:
        start,end=0,nrows_cols-1
    if interval>(end-start):
        interval=end-start    
    return start,end,interval
 
startrow,endrow,interval=set_rowscols(interval,rows,"rows")
startcol,endcol,interval=set_rowscols(interval,cols,"columns")


temp=np.fromfile(L1top,sep=" ")
L1top=np.reshape(temp,(rows,cols))
bots=np.fromfile(modelbottoms,sep=" ")
nlayers=np.size(bots)/np.size(L1top)
bots=np.reshape(bots,(nlayers,rows,cols))
all=np.append(L1top,bots)
all=np.reshape(all,(nlayers+1,rows,cols))

def make_imshow(layer,array,title):
    plt.figure()
    plt.imshow(array[layer,:,:])
    plt.colorbar()
    plt.title('%s %s' %(title,layer))


def make_planimages(fnames,arrays,titles):
    # figs: list of variable names, one for each multi-page PDF
    # arrays: list of numpy ndarrays, one for each multi-page PDF
    # titles: list of titles, one for each multi-page PDF
    
    nlay=len(arrays[0]) # number of model layers
    
    for c in range(len(fnames)):
        print fnames[c]
        fig=PdfPages(fnames[c])
        for i in range(nlay):
            make_imshow(i,arrays[c],titles[c])
            fig.savefig()
        fig.close()


def make_xsectplot(nlayers,slice1,slice2,interval,slice1_ind,X,label,xlabel,colors):
    # colors: list of html colors, one for each layer going from top to bottom
    
    
    plt.figure()
    ax1=plt.subplot(211)
    ax2=plt.subplot(212)
    
    for l in range(nlayers)[:-1]:
        ax1.fill_between(X,slice1[:,l-1],slice1[:,l],facecolor=colors[l-1])
        ax1.set_ylabel('Elevation (ft.)')
        
        ax1.text(0.99,0.05,label %(slice1_ind),verticalalignment='bottom',horizontalalignment='right',transform=ax1.transAxes)        
        
        if slice2<>None:
            ax2.fill_between(X,slice2[:,l-1],slice2[:,l],facecolor=colors[l-1])
            
            ax2.set_ylabel('Elevation (ft.)')
            ax2.text(0.99,0.05,label %(slice1_ind+interval),verticalalignment='bottom',horizontalalignment='right',transform=ax2.transAxes)
            #ax2.set_xlabel('Distance along section (miles)')
                
            ax2.set_xlabel(xlabel)            


def make_xsections(arrays,rc_range,interval,fnames,colors):
    # arrays: list of numpy ndarrays to plot
    # interval: index interval at which to make slices
    # fnames: list of filenames for PDFs to make (one per array)
    # colors: list of html colors, one per layer
    
    nlayers,nrows,ncols=np.shape(arrays[0])
    startrow,endrow,startcol,endcol=rc_range   
    
    for i in range(len(fnames)):
        
        fig=PdfPages(fnames[i])
        
        print "generating cross section plots at:"
        print "Rows: ",
        for c in range(nrows)[startrow:endrow:interval][::2]:
            print "%s,%s" %(c,c+interval),
            X=np.arange(ncols)
            #X=X*float(cellsize)/5280.0
            cslice1=np.transpose(arrays[i][:,c,:])
            try:
                cslice2=np.transpose(arrays[i][:,c+interval,:])
            except IndexError:
                cslice2=None
            
            label='Model row %s'
            xlabel='Column'
            
            make_xsectplot(nlayers,cslice1,cslice2,interval,c,X,label,xlabel,colors)
            
            fig.savefig()
        print "\n\nColumns: ",
        for c in range(ncols)[startcol:endcol:interval][::2]:
            print "%s,%s" %(c,c+interval),
            X=np.arange(nrows)
            #X=X*float(cellsize)/5280.0
            cslice1=np.transpose(arrays[i][:,:,c])
            try:
                cslice2=np.transpose(arrays[i][:,:,c+interval])
            except IndexError:
                cslice2=None
                
            label='Model column %s'
            xlabel='Row'
            
            make_xsectplot(nlayers,cslice1,cslice2,interval,c,X,label,xlabel,colors)
            fig.savefig()
        print "\n\nsaving to %s" %(fnames[i])
        fig.close()

# MAIN program        
make_xsections([all],[startrow,endrow,startcol,endcol],interval,[xsname],colors)
print "\nDone!"