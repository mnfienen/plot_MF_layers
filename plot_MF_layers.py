# functions to visualize layer structure in a MODFLOW model
# authored by Andrew Leaf, USGS - Wisconsin Water Science Center

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from collections import defaultdict
import discomb_utilities
import binaryfile_recarray as bf
import pdb

def set_rowscols(interval,nrows_cols,inputsname):
    try:
        start,end=map(int,inputs[inputsname].strip('[]').split(','))
        if start>0:
            start,end=start-1,end
        else:
            start,end=0,end
    except:
        start,end=0,nrows_cols-1
    if interval>(end-start):
        interval=end-start    
    return start,end,interval
 
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


def make_xsectplot(nlayers,slice1,slice2,interval,slice1_ind,X,label,xlabel,colors,**kwargs):
    # colors: list of html colors, one for each layer going from top to bottom
    try:
        Xunits=kwargs["Xunits"]
    except:
        Xunits=False

    try:
        linewidth=kwargs["linewidth"]
    except:
        linewidth=0.0

    try:
        WT1=kwargs["WT1"]
        watertable='On'
    except:
        watertable='Off'

    try:
        WTcolor = kwargs["WTcolor"]
    except:
        WTcolor = 'SkyBlue'
    try:
        WTthick = kwargs["WTthick"]
    except:
        WTthick = 1.0

    fig=plt.figure()
    ax1=plt.subplot(211)
    ax2=plt.subplot(212) 
    
    for c in range(nlayers)[:-1]:
        
        # in case there aren't enough colors
        try:
            colors[c-1]
            facecolor=colors[c-1]
        except IndexError:
            facecolor=colors[-1]
            
        ax1.fill_between(X,slice1[:,c-1],slice1[:,c],facecolor=facecolor,linewidth=linewidth)

        if c==0:
            # make room for label
            ax1.set_ylim(ax1.get_ylim()[0],ax1.get_ylim()[1]*1.1)             
            
        # figure out vertical exaggeration
        def vert_exaggeration(subplot,Xunits):
            
            defaultaspect=subplot.bbox.bounds[3]/subplot.bbox.bounds[2]
            ylen=subplot.get_ylim()[1]-subplot.get_ylim()[0]
            xlen=subplot.get_xlim()[1]-subplot.get_xlim()[0]
            
            if Xunits=='miles':
                VE=defaultaspect*(xlen*5280.0/ylen)
            else:
                VE=defaultaspect*(xlen/ylen)
            return VE
        
        VE1=vert_exaggeration(ax1,Xunits)

        ax1.set_ylabel('Elevation (ft.)')
        ax1.text(0.99,0.95,label %(slice1_ind+1,round(VE1,1)),verticalalignment='top',horizontalalignment='right',transform=ax1.transAxes)
        
        
        if slice2<>None:
            ax2.fill_between(X,slice2[:,c-1],slice2[:,c],facecolor=facecolor,linewidth=linewidth)
            
            if c==0:
                # make room for label
                ax2.set_ylim(ax2.get_ylim()[0],ax2.get_ylim()[1]*1.1)
                
            VE2=vert_exaggeration(ax2,Xunits)
            ax2.set_ylabel('Elevation (ft.)')
            ax2.text(0.99,0.95,label %(slice1_ind+1+interval,round(VE2,1)),verticalalignment='top',horizontalalignment='right',transform=ax2.transAxes)
            #ax2.set_xlabel('Distance along section (miles)')
                
            ax2.set_xlabel('%s' %(xlabel))       
        
        if Xunits=='cells':
            tickz=ax1.get_xticks()
            # convert back to cells
            tickz_cells=[]
            for t in tickz:
                idx = np.argmin(np.abs(X - t)) # find cell closest to the tick distance
                tickz_cells.append(idx)
            
            ax1.set_xticklabels(tickz_cells)
            ax2.set_xticklabels(tickz_cells)

    if watertable == "On":
        ax1.plot(X, WT1, linewidth=WTthick, color=WTcolor)
        if slice2<>None:
            ax2.plot(X, kwargs["WT2"], linewidth=WTthick, color=WTcolor)


def make_xsections(arrays,watertable,DX,DY,rc_range,row_interval,col_interval,fnames,colors,**kwargs):
    # arrays: list of numpy ndarrays to plot
    # DX,DY: column/row coordinates (in ft. or miles)
    # interval: index interval at which to make slices
    # fnames: list of filenames for PDFs to make (one per array)
    # colors: list of html colors, one per layer
    # Xunits
    try:
        Xunits=kwargs["Xunits"]
    except:
        Xunits=False
    try:
        linewidth=kwargs["linewidth"]
    except:
        linewidth=0.0

    try:
        WTcolor = kwargs["WTcolor"]
    except:
        WTcolor = 'SkyBlue'
    try:
        WTthick = kwargs["WTthick"]
    except:
        WTthick = 1.0

    nlayers,nrows,ncols=np.shape(arrays[0])
    startrow,endrow,startcol,endcol=rc_range   
    
    for i in range(len(fnames)):
        
        fig=PdfPages(fnames[i])
        
        print "generating cross section plots at:"
        print "Rows: ",
        for c in range(nrows)[startrow:endrow:row_interval][::2]:
            print "%s,%s" %(c+1,c+row_interval+1),

            cslice1=np.transpose(arrays[i][:,c,:])
            try:
                cslice2=np.transpose(arrays[i][:,c+row_interval,:])
            except IndexError:
                cslice2=None
                WT2=None
            if Xunits=='cells':    
                xlabel='Column'
            else:
                xlabel='Distance (%s) ' %(Xunits)            
                
            label='Model row %s; VE: %sx'
            try:
                WT1 = watertable[c,:]
                if cslice2 != None:
                    WT2 = watertable[c+row_interval,:]
                make_xsectplot(nlayers,cslice1,cslice2,row_interval,c,DX,label,xlabel,colors,WT1=WT1,WT2=WT2,
                               Xunits=Xunits,linewidth=linewidth,WTcolor=WTcolor,WTthick=WTthick)
            except:
                make_xsectplot(nlayers,cslice1,cslice2,row_interval,c,DX,label,xlabel,colors,Xunits=Xunits,linewidth=linewidth)
            
            fig.savefig()
        print "\n\nColumns: ",
        for c in range(ncols)[startcol:endcol:col_interval][::2]:
            print "%s,%s" %(c+1,c+col_interval+1),

            cslice1=np.transpose(arrays[i][:,:,c])
            try:
                cslice2=np.transpose(arrays[i][:,:,c+col_interval])
            except IndexError:
                cslice2=None
                WT2=None
            if Xunits=='cells':
                xlabel='Row'  
            else:
                xlabel='Distance (%s) ' %(Xunits)            
                
            label='Model column %s; VE: %sx'
            
            try:
                WT1 = watertable[:,c]
                if cslice2 != None:
                    WT2 = watertable[:,c+row_interval]
                make_xsectplot(nlayers,cslice1,cslice2,col_interval,c,DY,label,xlabel,colors,WT1=WT1,WT2=WT2,
                               Xunits=Xunits,linewidth=linewidth,WTcolor=WTcolor,WTthick=WTthick)
            except:
                make_xsectplot(nlayers,cslice1,cslice2,col_interval,c,DY,label,xlabel,colors,Xunits=Xunits,linewidth=linewidth)
            fig.savefig()
        print "\n\nsaving to %s" %(fnames[i])
        fig.close()

# MAIN program
# input file
infile="plot_MF_layers.in"

# Get input files locations
infile=open(infile,'r').readlines()
inputs=defaultdict()
for line in infile:
    try:
        if "#" not in line.split()[0]:
            varname=line.split("=")[0]
            var=line.split("=")[1].split()[0]
            inputs[varname]=var.replace("'","") # strip extra quotes
    except IndexError: # means that line is empty
        continue

modelbottoms=inputs["modelbottoms"]
L1top=inputs["L1top"]
xsname=inputs["xsname"]
colors=inputs["colors"].strip('[]').split(',')
interval=int(inputs["interval"])
DISfile=inputs["DIS_file"]
Xunits=inputs["Xunits"]
linewidth=float(inputs["linewidth"])
try:
    watertable=inputs["watertable"]
    HDSfile=inputs["HDS_file"]
    WTcolor=inputs["WTcolor"]
    WTthick=float(inputs["WTthickness"])
except:
    watertable='Off'



# get dimmensions
DX, DY, NLAY, NROW, NCOL, i = discomb_utilities.read_meta_data(DISfile)
if Xunits=='miles':
    DX,DY = DX/5280.0, DY/5280.0

startrow, endrow, row_interval = set_rowscols(interval, NROW, "rows")
startcol, endcol, col_interval = set_rowscols(interval, NCOL, "columns")


#temp=np.fromfile(L1top,sep=" ")
#L1top=np.reshape(temp,(NROW,NCOL))
#bots=np.fromfile(modelbottoms,sep=" ")

# get layer tops/bottoms from DIS file
layer_elevs=np.zeros((NLAY+1, NROW, NCOL))
for c in range(NLAY+1):
    tmp, i = discomb_utilities.read_nrow_ncol_vals(DISfile, NROW, NCOL, 'float', i)
    layer_elevs[c, :, :] = tmp
'''
#nlayers=np.size(bots)/np.size(L1top)
bots=np.reshape(bots, (NLAY, NROW, NCOL))
all=np.append(L1top, bots)
all=np.reshape(all, (NLAY+1, NROW, NCOL))
'''
if watertable.lower()=='on':
    hds = bf.HeadFile(HDSfile)
    heads = hds.get_data(kstp=1, kper=1)
    watertable = heads[0,:,:]


make_xsections([layer_elevs],watertable,DX,DY,[startrow,endrow,startcol,endcol],row_interval,col_interval,[xsname],colors,
               Xunits=Xunits,linewidth=linewidth,WTcolor=WTcolor,WTthick=WTthick)
print "\nDone!"