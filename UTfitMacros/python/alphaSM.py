import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from matplotlib.patches import Patch
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, TextArea, VPacker
from matplotlib.patches import Rectangle
import matplotlib.image as mpimg

# Parameters
pipi_alphaSM = 81.9
spipi_alphaSM = 3.9

rhorho_alphaSM = 94.7
srhorho_alphaSM = 3.8

#average: 88.4662 +/- 2.72167
#norm chi2: 5.5258 so S: 2.3507
#corrected average: 88.4662 +/- 6.39784

mupdg_alphaSM = 88.5 
spdg_alphaSM = 6.4

all_alphaSM1 = 88.0 
sall_alphaSM1 = 3.6
all_alphaSM2 = 102.9
sall_alphaSM2 = 5.0

pred_alphaSM = 91.7
spred_alphaSM = 1.4

num_samples = 10000000

# Define blue gradient colours
#colours = ['#6baed6', '#2171b5', '#08306b']  # dark to light blue
#colours = ['#3399FF', '#66CCFF', '#99DDFF']  # bright blue gradient
#colours = ['#0066CC', '#3399FF', '#66CCFF']
#colours = ['#005FCC', '#008CFF', '#00BFFF'] 

blues = ['#1E26C8',  # darker blue
          '#3137FD',  # base vivid blue
          '#7F84FF']  # lighter blue

#bzyellow = '#F7B500'
yellows = '#FFD700'
oranges = 'orangered'

colours = ['gray', blues[0], blues[2], oranges, yellows]
alphas = [0.5, 0.8, 0.5, 0.7, 0.6] 

# X range covering both distributions
x = np.linspace(60., 130., 500)

# Compute PDF values
pipi_pdf = norm.pdf(x, pipi_alphaSM, spipi_alphaSM)
rhorho_pdf = norm.pdf(x, rhorho_alphaSM, srhorho_alphaSM)
pdg_pdf = norm.pdf(x, mupdg_alphaSM, spdg_alphaSM)

# Compute the individual PDFs
all1_pdf = norm.pdf(x, all_alphaSM1, sall_alphaSM1)
all2_pdf = norm.pdf(x, all_alphaSM2, sall_alphaSM2)

# Sum of the two distributions
all_pdf = 0.4*all1_pdf + 0.6*all2_pdf 

# 
pred_pdf = norm.pdf(x, pred_alphaSM, spred_alphaSM)

# Plotting
fig, ax = plt.subplots(figsize=(12,12))

#plt.plot(x, pred_pdf, color=colours[0])
#plt.fill_between(x, pred_pdf, color=colours[0], alpha=alphas[0])
plt.plot(x, pipi_pdf, color=colours[1])
plt.fill_between(x, pipi_pdf, color=colours[1], alpha=alphas[1])
plt.plot(x, rhorho_pdf, color=colours[2])
plt.fill_between(x, rhorho_pdf, color=colours[2], alpha=alphas[2])
plt.plot(x, pdg_pdf, color=colours[3])
plt.fill_between(x, pdg_pdf, color=colours[3], alpha=alphas[3])
#plt.plot(x, all_pdf, color=colours[4])
#plt.fill_between(x, all_pdf, color=colours[4], alpha=alphas[4])
plt.plot(x, all1_pdf, color=colours[4])
plt.fill_between(x, all1_pdf, color=colours[4], alpha=alphas[4])
plt.xlabel("x")
plt.ylabel("Probability Density")
plt.ylim(bottom=0)  # Set minimum of y-axis to zero

ax.set_xlabel(r"$\alpha_{SM}\;[^{\circ}]$", fontsize=30)
ax.set_ylabel(r"Probability Density", fontsize=28)
# Move x label to right end and y label to top end
ax.xaxis.set_label_coords(0.90, -0.05)  # (x, y) in axis-relative coordinates                                                                                                          
ax.yaxis.set_label_coords(-0.07, 0.76)
ax.tick_params(axis='both', which='major', labelsize=16)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontweight('bold')

handles = [Rectangle((0,0),1,1, color=c, alpha=a) for c, a in zip(colours, alphas)]
labels = [r'UTfit Pred.', r'$\mathbf{\pi\pi}$', r'$\mathbf{\rho\rho}$', 'PDG-style Average', 'UTfit SM Comb.']

#plt.legend(handles, labels, prop={'size':20,'weight':'bold'})

#selected_handles = [handles[0], handles[-1]]
#selected_labels = [labels[0], labels[-1]]

selected_handles = [handles[1], handles[2], handles[3], handles[4]]
selected_labels = [labels[1], labels[2], labels[3], labels[4]]
plt.legend(selected_handles, selected_labels, prop={'size':20, 'weight':'bold'})

# Load logo
logo = mpimg.imread('../common/logo.png')
imagebox = OffsetImage(logo, zoom=0.17)
# Create TextArea for the label below the logo
#text = TextArea("summer25", textprops=dict(color="#cc0000", fontsize=20, fontweight='bold'))
text = TextArea("summer25", textprops=dict(color="#1434A4", fontsize=20, fontweight='bold'))

# Combine logo and text vertically
logo_and_text = VPacker(children=[imagebox, text],
                        align="center", pad=0, sep=5)
#ab = AnnotationBbox(imagebox, (0.22, 0.98), xycoords='axes fraction',
#                    frameon=False, box_alignment=(1,1))

# Create AnnotationBbox to place it in axes fraction coordinates
ab = AnnotationBbox(logo_and_text, (0.22, 0.98),
                    xycoords='axes fraction',
                    frameon=False,
                    box_alignment=(1,1))

ax.add_artist(ab)

plt.savefig("alphaSM_2025.png", dpi=300, bbox_inches='tight')
plt.savefig("alphaSM_2025.pdf", bbox_inches='tight')

plt.show()
