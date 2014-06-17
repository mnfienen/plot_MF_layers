# functions to visualize layer structure in a MODFLOW model
# authored by Andrew Leaf, USGS - Wisconsin Water Science Center

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import discomb_utilities
import binaryfile_recarray as bf
import sys
import xml.etree.ElementTree as ET



newparams = {'legend.fontsize':12, 'axes.labelsize':12, 'xtick.labelsize':12, 'ytick.labelsize':12,
             'font.sans-serif':'Univers 57 Condensed',
             'pdf.fonttype':42, 'pdf.compression':0}

mpl.rcParams.update(newparams)


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

def tf2flag(intxt):
    # converts text written in XML file to True or False flag when parsing XML input
    if intxt.lower()=='true':
        return True
    else:
        return False

def make_imshow(layer,array,title):
    plt.figure()
    plt.imshow(array[layer,:,:])
    plt.colorbar()
    plt.title('%s %s' %(title,layer))

def make_xsectplot(nlayers,slice1,slice1_ind,X,label,xlabel,colors,VE,**kwargs):
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
    ax = plt.subplot(111)
    if VE < 0:
        VE=vert_exaggeration(ax,Xunits)
    else:
        VE_x_conversion = 1.0
        if Xunits.lower() == 'cells':
            VE_x_conversion = 1/np.mean(DX)
        elif Xunits.lower() == 'miles':
            VE_x_conversion = 1/5280
        #VE *= VE_x_conversion

        
    for c in range(nlayers)[:-1]:
        
        # in case there aren't enough colors
        try:
            colors[c-1]
            facecolor=colors[c-1]
        except IndexError:
            facecolor=colors[-1]
            
        ax.fill_between(X,slice1[:,c-1],slice1[:,c],facecolor=facecolor,linewidth=linewidth)

        if c==0:
            # make room for label
            ax.set_ylim(ax.get_ylim()[0],ax.get_ylim()[1]*1.1)             


        plt.ylabel('Elevation (ft.)')
        plt.xlabel(xlabel)
        ax.text(0.99,0.95,label %(slice1_ind+1,VE),verticalalignment='top',
                 horizontalalignment='right',transform=ax.transAxes)
        
        


        if Xunits=='cells':
            tickz=ax.get_xticks()
            # convert back to cells
            tickz_cells=[]
            for t in tickz:
                idx = np.argmin(np.abs(X - t)) # find cell closest to the tick distance
                tickz_cells.append(idx)

            ax.set_xticklabels(tickz_cells)

        if watertable == "On":
            ax.plot(X, WT1, linewidth=WTthick, color=WTcolor)
        ax.set_aspect(VE, adjustable='box')

def make_xsections(arrays,watertable,DX,DY,rc_range,row_interval,col_interval,fnames,colors,VE,**kwargs):
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
        for c in range(startrow,endrow,row_interval):
            print "%s,%s" %(c+1,c+row_interval+1),

            cslice1=np.transpose(arrays[i][:,c,:])

            if Xunits=='cells':    
                xlabel='Column'
            else:
                xlabel='Distance (%s) ' %(Xunits)            
                
            label='Model row %s; VE: %sx'
            try:
                WT1 = watertable[c,:]

                make_xsectplot(nlayers,cslice1,row_interval,c,DX,label,xlabel,colors,VE,WT1=WT1,
                               Xunits=Xunits,linewidth=linewidth,WTcolor=WTcolor,WTthick=WTthick)
            except:
                make_xsectplot(nlayers,cslice1,c,DX,label,xlabel,colors,VE,Xunits=Xunits,
                               linewidth=linewidth)
            
            fig.savefig()
            plt.close('all')
        print "\n\nColumns: ",
        for c in range(startcol,endcol,col_interval):
            print "%s,%s" %(c,c+col_interval),

            cslice1=np.transpose(arrays[i][:,:,c])

            if Xunits=='cells':
                xlabel='Row'  
            else:
                xlabel='Distance (%s) ' %(Xunits)            
                
            label='Model column %s; VE: %sx'
            
            try:
                WT1 = watertable[:,c]

                make_xsectplot(nlayers,cslice1,c,DY,label,xlabel,colors,VE,WT1=WT1,
                               Xunits=Xunits,linewidth=linewidth,WTcolor=WTcolor,WTthick=WTthick)
            except:
                make_xsectplot(nlayers,cslice1,c,DY,label,xlabel,colors,VE,
                               Xunits=Xunits,linewidth=linewidth)
            fig.savefig()

        print "\n\nsaving to {0}".format(fnames[i])
        fig.close()

# MAIN program
# input file
try:
    infile = sys.argv[1]
except:
    infile = "plot_MF_layers.in"

# Read in the XML input file
##
inpardat = ET.parse(infile)

inpars = inpardat.getroot()
UseDis = tf2flag(inpars.findall('.//layers/UseDis')[0].text)
if UseDis == False:
    modelbottoms = inpars.findall('.//layers/modelbottoms')[0].text
    L1top = inpars.findall('.//layers/L1top')[0].text

xsname=inpars.findall('.//files/xsname')[0].text
colors = []
for ccol in inpars.findall('.//colors/color'):
    colors.append(ccol.text)
interval=int(inpars.findall('.//options/interval')[0].text)
DISfile=inpars.findall('.//files/DIS_file')[0].text
Xunits=inpars.findall('.//options/Xunits')[0].text
linewidth = float(inpars.findall('.//options/linewidth')[0].text)

startrow = int(inpars.findall('.//rows/startrow')[0].text)
endrow = int(inpars.findall('.//rows/endrow')[0].text)
startcol = int(inpars.findall('.//columns/startcolumn')[0].text)
endcol = int(inpars.findall('.//columns/endcolumn')[0].text)

if inpars.findall('.//options/watertable')[0].text.lower() == "on":
    watertable = "On"
    HDSfile=inpars.findall('.//WaterTableSettings/HDS_file')[0].text
    WTcolor=iinpars.findall('.//WaterTableSettings/WTColor')[0].text
    WTthick=float(inpars.findall('.//WaterTableSettings/WTthickness')[0].text)
    # read in water table values
else:
    watertable = "Off"
try:
    VE = float(inpars.findall('.//options/VE')[0].text)
except:
    VE = -999

# get dimmensions
DX, DY, NLAY, NROW, NCOL, i = discomb_utilities.read_meta_data(DISfile)
if Xunits=='miles':
    DX,DY = DX/5280.0, DY/5280.0

row_interval = interval
col_interval = interval

# get layer tops/bottoms from DIS file
layer_elevs=np.zeros((NLAY+1, NROW, NCOL))
for c in range(NLAY+1):
    tmp, i = discomb_utilities.read_nrow_ncol_vals(DISfile, NROW, NCOL, 'float', i)
    layer_elevs[c, :, :] = tmp


if watertable.lower()=='on':
    hds = bf.HeadFile(HDSfile)
    heads = hds.get_data(kstp=1, kper=1)
    watertable = heads[0,:,:]
else:
    WTcolor = None
    WTthick = None

make_xsections([layer_elevs],watertable,DX,DY,[startrow,endrow,startcol,endcol],
               row_interval,col_interval,[xsname],colors,VE,
               Xunits=Xunits,linewidth=linewidth,WTcolor=WTcolor,WTthick=WTthick)
print "\nDone!"