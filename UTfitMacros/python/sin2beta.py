import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from matplotlib.patches import Patch
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, TextArea, VPacker
from matplotlib.patches import Rectangle
import matplotlib.image as mpimg

# Parameters
sin2beta = 0.700
ssin2beta = 0.015

pred_sin2beta = 0.768
spred_sin2beta = 0.029

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
x = np.linspace(0.6, 0.85, 500)
dx = x[1] - x[0]

# Compute PDF values
sin2beta_pdf = norm.pdf(x, sin2beta, ssin2beta)
integral = np.sum(sin2beta_pdf * dx)
sin2beta_pdf /= integral

# 
pred_pdf = norm.pdf(x, pred_sin2beta, spred_sin2beta)
integral = np.sum(pred_pdf * dx)
pred_pdf /= integral

# Plotting
fig, ax = plt.subplots(figsize=(12,12))

plt.plot(x, sin2beta_pdf, color=colours[0])
plt.fill_between(x, sin2beta_pdf, color=colours[0], alpha=alphas[0])
plt.plot(x, pred_pdf, color=colours[1])
plt.fill_between(x, pred_pdf, color=colours[1], alpha=alphas[1])
plt.xlabel("x")
plt.ylabel("Probability Density")
plt.ylim(bottom=0)  # Set minimum of y-axis to zero

ax.set_xlabel(r"$\sin 2\beta$", fontsize=30)
ax.set_ylabel(r"Probability Density", fontsize=28)
# Move x label to right end and y label to top end
ax.xaxis.set_label_coords(0.90, -0.05)  # (x, y) in axis-relative coordinates                                                                                                          
ax.yaxis.set_label_coords(-0.07, 0.76)
ax.tick_params(axis='both', which='major', labelsize=16)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontweight('bold')

handles = [Rectangle((0,0),1,1, color=c, alpha=a) for c, a in zip(colours, alphas)]
labels = ['UTfit Pred.', 'UTfit SM Comb.']

#plt.legend(handles, labels, prop={'size':20,'weight':'bold'})

#selected_handles = [handles[0], handles[-1]]
#selected_labels = [labels[0], labels[-1]]

selected_handles = [handles[0], handles[1]]
selected_labels = [labels[0], labels[1]]
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

plt.savefig("sin2beta_2025.png", dpi=300, bbox_inches='tight')
plt.savefig("sin2beta_2025.pdf", bbox_inches='tight')

plt.show()
