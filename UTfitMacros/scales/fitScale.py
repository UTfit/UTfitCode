import matplotlib.pyplot as plt
import matplotlib.patches as patches
import argparse
import numpy as np
import pandas as pd
import math
import json
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, TextArea, VPacker
from matplotlib.patches import Patch

gray = ["lightgray", "darkgray", "gray"]
colors = ["#8dd3c7", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f"]

#define style
style = {"pgf.rcfonts":False,
         "pgf.texsystem": "pdflatex",   
         "text.usetex": True,                
         "font.family": "sans-serif",
         "font.size": 20
        }
#set
plt.rcParams.update(style)

# set font in latex
#plt.rcParams['text.latex.preamble'] = (
#    r'\usepackage{siunitx}'
#    r'\sisetup{detect-all}'
#    r'\usepackage{helvet}'
#    r'\usepackage{sansmath}'
#    r'\sansmath'
#)

plt.rcParams['text.latex.preamble'] = (
    r'\usepackage{siunitx}'
    r'\sisetup{detect-all}'
    r'\usepackage{mathastext}'
    r'\MTfamily{lmss}'   # Latin Modern Sans
)

def read_input_txt(inputFileName):
    inputFile = open(inputFileName, "r")
    scenario = ""
    operators = []
    observables = []

    for line in inputFile:
        if line.find("scenario") != -1:
            scenario = line[:-1].replace("scenario: ", "")
        if line.find("operators") != -1:
            operator_names = line[:-1].replace("operators: ", "")
            operators = operator_names.split(", ")

    n_operators = len(operators)
    inputFile.close()

    inputFile = open(inputFileName, "r")
    observables = []
    values = []
    for line in inputFile:
        if line.find("scenario") != -1: continue
        if line.find("operators") != -1: continue
        observables.append(line.split(", ")[0])
        print(observables)
        these_values = []
        for value in line.split(", ")[1:]:
            these_values.append(float(value))
        print(these_values)
        values.append(these_values)
    inputFile.close()
    
    return scenario, observables, operators, np.array(values)
    
def drawBars(ax, value, scenarios, observables, operators, plot_name, alpha, fill=True, inLegend=True):

    width = 0.15
    nObs = len(observables)
    nOps = len(operators)
    x = np.arange(nOps)
    for i in range(nObs):
        if inLegend==True:
            ax.bar(x + (i-int(nOps/2.)+0.5)*width, np.ndarray.flatten(value[i]), width,
                   label=observables[i], alpha=alpha, color=colors[i], edgecolor=colors[i], linewidth=2, fill=fill)
        else:
            ax.bar(x + (i-int(nOps/2.)+0.5)*width, np.ndarray.flatten(value[i]), width,
                   label='_nolegend_', alpha=alpha, color=colors[i], edgecolor=colors[i], linewidth=2, fill=fill)
            

    #ax.set_title(' ')
    ax.set_xticks(x)
    ax.set_xticklabels(operators)

    return ax

def draw_MultipleQuantities(value, plot_name, utfit):
    fig, ax = plt.subplots()

    plt.savefig(plot_name)

    
def drawHEPfitLogo(ax):
    x_limits= ax.get_xlim()
    y_limits= ax.get_ylim()

    xLogo = x_limits[0] + (x_limits[1]-x_limits[0])*0.905
    yLogo = y_limits[0] + math.pow(10., math.log10(y_limits[1])*0.93)
        
    logo = plt.imread('../common/hepfit_logo.png')
    imagebox = OffsetImage(logo, zoom = 0.08)
    ab = AnnotationBbox(imagebox, (xLogo, yLogo), frameon = False)
    ax.add_artist(ab)
    
    return ax

def drawUTfitLogo(ax):
    x_limits= ax.get_xlim()
    y_limits= ax.get_ylim()

    xLogo = x_limits[0] + (x_limits[1]-x_limits[0])*0.89
    yLogo = y_limits[0] + math.pow(10., math.log10(y_limits[1])*0.86)

    logo = plt.imread('../common/logo.png')
    imagebox = OffsetImage(logo, zoom = 0.1)
    #ab = AnnotationBbox(imagebox, (xLogo, yLogo), frameon = False)

    text='summer25'
    text_size=18
    
    if text is not None:
        # Create a TextArea for the caption
        text_area = TextArea(text, textprops=dict(size=text_size, color="#1434A4", fontweight='bold'))
        # Stack image and text vertically
        vpack = VPacker(children=[imagebox, text_area], align="center", pad=0, sep=2)
        ab = AnnotationBbox(vpack, (xLogo, yLogo), frameon=False)
    else:
        # Only logo
        ab = AnnotationBbox(imagebox, (xLogo, yLogo), frameon=False)
    
    ax.add_artist(ab)
    
    return ax

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i1','--input1', type=str, help="First input file", required=True)
    parser.add_argument('-i2','--input2', type=str, help="Second input file", default="none")
    parser.add_argument('-i3','--input3', type=str, help="Third input file", default="none")
    parser.add_argument('-o','--output', type=str, help="output file", default="test")
    parser.add_argument('-t','--title', type=str, help="figure title", default=" ")
    parser.add_argument('--utfit', action='store_true',  help="Use utfit style")
    args = parser.parse_args()

    if args.utfit == True:
        green = '#322B8F'
        red = '#D7322C'
    
    scenario1, observables, operators, values1 = read_input_txt(args.input1)
    scenario2 = ""
    values2 = np.array([])
    scenario3 = ""
    values3 = np.array([])
    scenarios = [scenario1]
    values = [values1]
    if args.input2 != "none":
        scenario2, obs_null, ops_null, values2 = read_input_txt(args.input2)
        scenarios.append(scenario2)
        values.append(values2)
    if args.input3 != "none":
        scenario3, obs_null, ops_null, values3 = read_input_txt(args.input3)
        scenarios.append(scenario3)
        values.append(values3)

    fig, ax = plt.subplots(figsize=(8, 6), dpi=120)

    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['mathtext.default'] = 'it'
#    ax.xaxis.get_major_formatter()._usetex = False
#    ax.yaxis.get_major_formatter()._usetex = True
    ax.set_ylim(0.1, np.max(values)*100.)
    ax.set_yscale('log')
    ax.set_ylabel(r'$\Lambda$ (TeV)',fontsize="26")

    alpha = [0.03, 0.05, 1.0]

    #re-draw lines
    if len(values)>1:
        ax = drawBars(ax, values[0], scenarios, observables, operators, args.output, alpha[3-len(values)], True, False)
        ax = drawBars(ax, values[0], scenarios, observables, operators, args.output, 1.0, False, False)
        ax = drawBars(ax, values[1], scenarios, observables, operators, args.output,  alpha[4-len(values)])
        if len(values)>2: ax = drawBars(ax, values[2], scenarios, observables, operators, args.output, alpha[2], True, False)
    else:
        ax = drawBars(ax, values[0], scenarios, observables, operators, args.output, alpha[3-len(values)])
        
    x_limits= ax.get_xlim()
    y_limits= ax.get_ylim()

    # draw grid lines
    max_y = int(math.log10(np.max(values)*100.))-1
    for i in range(1,max_y):
        ax.plot((-2,2.*len(operators)),(math.pow(10.,i), math.pow(10.,i)), linestyle="--", color="lightgray", linewidth = 1, zorder=0)
    ax.set_xlim(x_limits)
    plt.tick_params(left = False)
    plt.tick_params(bottom  = False)
            
    leg1 = plt.legend(loc="upper left", ncol=3, fontsize="16", frameon=False)
    plt.gca().add_artist(leg1)
    
    if len(values) == 1:
        ax.set_title(scenario1)
    elif len(values) == 2:
        legend_elements = [Patch(facecolor="white", edgecolor='dimgray', linewidth=2, label= scenario1),
                           Patch(facecolor=gray[2], edgecolor='dimgray', linewidth=2, label= scenario2)]
        leg2 = plt.legend(handles=legend_elements, loc="upper left", bbox_to_anchor=(0.0, 1.1), ncol=len(values), fontsize="16", frameon=False)
        plt.gca().add_artist(leg2)
    else:
        legend_elements = [Patch(facecolor="white", edgecolor='dimgray', linewidth=2, label= scenario1),
                           Patch(facecolor=gray[1], edgecolor='dimgray', linewidth=2, label= scenario2),
                           Patch(facecolor=gray[2], edgecolor='dimgray', linewidth=2, label= scenario3)]
        leg2 = plt.legend(handles=legend_elements, loc="upper left", bbox_to_anchor=(0.0, 1.1), ncol=len(values), fontsize="16", frameon=False)
        plt.gca().add_artist(leg2)

    if  args.utfit == True:
        drawUTfitLogo(ax)
    else:
        drawHEPfitLogo(ax)


    #plt.savefig(args.output, bbox_inches='tight')
    plt.savefig(args.output+".png", dpi=300, bbox_inches='tight')
    plt.savefig(args.output+".pdf", bbox_inches='tight')


