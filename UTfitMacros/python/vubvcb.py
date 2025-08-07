import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import multivariate_normal
from scipy.optimize import brentq
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, TextArea, VPacker
import matplotlib.image as mpimg

# UTfit-lattice_202411.txt
# VcbExcl = 0.04012 +- 0.00055 (Guido, Ludovico & Silvano 2105.08674, 2204.05925, 2310.03680)
# VcbIncl = 0.04197 +- 0.00048 (Finauri&Gambino arXiv:2310.20324 [hep-ph])
# VubExcl = 0.00363 +- 0.00026 (Guido, Ludovico & Silvano update of 2202.10285)
# VubIncl = 0.00413 +- 0.00026 (PDG2025)
# VuboVcb = 0.087 +- 0.009 (Bs to Klnu & Bs to Dslnu: arXiv:2310.03680) 
# 
# Skeptic 2D average:
# UTfit-lattice_202411.txt
# from Silvano: without B to taunu
# Vcb = 0.04118 +- 0.00076
# Vub = 0.00382 +- 0.00034
# correlation 0.04
#
# measurements from full fit
# Vcb = 0.04187 +- 0.00037
# Vub = 0.00374 +- 0.00008

vcb    = 0.04118 * 1000. 
svcb   = 0.00076 * 1000. 
vub    = 0.00382 * 1000. 
svub   = 0.00034 * 1000. 
rho    = 0.043 

vcbff  = 0.04187 * 1000. 
svcbff = 0.00037 * 1000. 
vubff  = 0.00374 * 1000. 
svubff = 0.00008 * 1000.

vcbmin = 0.03 * 1000. 
vcbmax = 0.048 * 1000. 
vubmin = 0.0020 * 1000.
vubmax = 0.0060 * 1000.

ratio = (vubmax-vubmin)/(vcbmax-vcbmin)

VcbExcl = 0.04012 * 1000
sVcbExcl = 0.00055  * 1000
VcbIncl = 0.04197 * 1000
sVcbIncl = 0.00048 * 1000
VubExcl = 0.00363 * 1000
sVubExcl = 0.00026 * 1000
VubIncl = 0.00413 * 1000
sVubIncl = 0.00026 * 1000
VuboVcb = 0.087 
sVuboVcb = 0.009 

colour1 = '#D4B062'
colour2 = '#DB3318'
colour3 = 'red'
colour4 = 'darkred'

##############################
# Define grid for evaluation
x = np.linspace(vcbmin, vcbmax, 500)
y = np.linspace(vubmin, vubmax, 500)
X, Y = np.meshgrid(x, y)

##############################
# Define the inputs with correlation
mean = [vcb, vub]
sigma_x = svcb
sigma_y = svub
cov = [[sigma_x**2, rho * sigma_x * sigma_y],
       [rho * sigma_x * sigma_y, sigma_y**2]]

rv = multivariate_normal(mean, cov)
pos = np.dstack((X, Y))
pdf = rv.pdf(pos)

# Function to compute the cumulative probability above a given threshold
def integral_above_threshold(thresh):
    return pdf[pdf >= thresh].sum() * (x[1]-x[0]) * (y[1]-y[0])

# Function to find the contour level for a given target probability mass
def find_level(target_prob):
    def func(thresh):
        return integral_above_threshold(thresh) - target_prob
    return brentq(func, 1e-10, np.max(pdf))

# Compute contour levels for 68% and 95% confidence
level_68 = find_level(0.68)
level_95 = find_level(0.95)

levels = [0, level_95, level_68, np.max(pdf)]
colors = [colour4, colour3]


##############################
# Define the full fit contours
mean2 = [vcbff, vubff]
cov2 = [[svcbff**2, 0], [0, svubff**2]]

# Create multivariate normal distribution
rv2 = multivariate_normal(mean2, cov2)
pos2 = np.dstack((X, Y))
pdf2 = rv2.pdf(pos2)

# Flatten and sort the pdf to integrate from top
pdf_flat = pdf2.flatten()
pdf_sort = np.sort(pdf_flat)[::-1]
pdf_cumsum = np.cumsum(pdf_sort)
pdf_cumsum /= pdf_cumsum[-1]

# Find levels corresponding to 68% and 95% enclosed probability
level_68ff = pdf_sort[np.searchsorted(pdf_cumsum, 0.68)]
level_95ff = pdf_sort[np.searchsorted(pdf_cumsum, 0.95)]

levels2 = [0, level_95ff, level_68ff, np.max(pdf2)]
colors2 = ['darkgrey', 'grey']

fig, ax = plt.subplots(figsize=(12,12))

# Vertical red band (note: x is multiplied by 1000)
x1 = VcbExcl
xw1 = sVcbExcl 
ax.axvspan(x1 - xw1, x1 + xw1, color=colour2, alpha=0.6)
# Add a vertical line at x_center
ax.axvline(x1, color=colour2, linestyle='-', linewidth=2)

ax.text(x1-xw1, 5.5, 'Exclusive', rotation=90,
        verticalalignment='center', fontsize=20, fontweight='bold', color='black')

# Vertical red band (note: x is multiplied by 1000)
x2 = VcbIncl
xw2 = sVcbIncl 
ax.axvspan(x2 - xw2, x2 + xw2, color=colour1, alpha=0.6)
# Add a vertical line at x_center
ax.axvline(x2, color=colour1, linestyle='-', linewidth=2)

ax.text(x2 - xw2, 5.5, 'Inclusive', rotation=90,
        verticalalignment='center', fontsize=20, fontweight='bold', color='black')

# Horizontal red band at y = -2 with width 1 (y multiplied by 1000)
y1 = VubExcl
yw1 = sVubExcl 
ax.axhspan(y1 - yw1, y1 + yw1, color=colour2, alpha=0.6)
ax.axhline(y1, color=colour2, linestyle='-', linewidth=2)

ax.text(32, y1 + yw1/5, 'Exclusive',
        verticalalignment='center', fontsize=20, fontweight='bold', color='black')

# Horizontal red band at y = -2 with width 1 (y multiplied by 1000)
y2 = VubIncl
yw2 = sVubIncl 
ax.axhspan(y2 - yw2, y2 + yw2, color=colour1, alpha=0.6)
ax.axhline(y2, color=colour1, linestyle='-', linewidth=2)

ax.text(32, y2 + yw2/4, 'Inclusive',
        verticalalignment='center', fontsize=20, fontweight='bold', color='black')

# Sort x to ensure fill_between works correctly
x_sorted = np.sort(x)

# Define the generalized band boundaries
y_central = VuboVcb * x_sorted
y_lower_bound = (VuboVcb - sVuboVcb) * x_sorted
y_upper_bound = (VuboVcb + sVuboVcb) * x_sorted

# Fill the band
plt.fill_between(x_sorted, y_lower_bound, y_upper_bound, color='green', alpha=0.3)
# Plot the band boundaries (optional, but good for clarity)
plt.plot(x_sorted, y_lower_bound, color='green', linestyle='--')
plt.plot(x_sorted, y_upper_bound, color='green', linestyle='--')
plt.plot(x_sorted, y_central, color='darkgreen', linestyle='-', linewidth=2)

ax.text(31, 3.0, r"$\mathbf{V_{ub}/V_{cb}}$ ($\mathbf{B_{s}}$ to $\mathbf{K}$)", rotation=16.2,
        verticalalignment='center', fontsize=20, fontweight='bold', color='black')

# Plotting inputs
contourf = ax.contourf(X, Y, pdf, levels=levels, colors=colors, alpha=0.4)
contourf.collections[0].set_alpha(0.0)  # Make lowest region transparent
contours = ax.contour(X, Y, pdf, levels=[level_95, level_68], colors=[colour3,colour4], linewidths=4)
ax.set_aspect('auto')
ax.clabel(contours, fmt={level_95: "95%", level_68: "68%"}, inline=True, fontsize=14)

# Plotting full fit
contourf2 = ax.contourf(X, Y, pdf2, levels=levels2, colors=colors2, alpha=0.4)
contourf2.collections[0].set_alpha(0.0)  # Make lowest region transparent
contours2 = ax.contour(X, Y, pdf2, levels=[level_95ff, level_68ff], colors=['grey','black'], linewidths=4)
ax.set_aspect('auto')
ax.clabel(contours2, fmt={level_95ff: "95%", level_68ff: "68%"}, inline=True, fontsize=14)

ax.set_xlabel(r"$|\mathbf{V_{cb}}|[10^{-3}]$", fontsize=28)
ax.set_ylabel(r"$|\mathbf{V_{ub}}|[10^{-3}]$", fontsize=28)
# Move x label to right end and y label to top end
ax.xaxis.set_label_coords(0.88, -0.05)  # (x, y) in axis-relative coordinates
ax.yaxis.set_label_coords(-0.07, 0.86)
ax.tick_params(axis='both', which='major', labelsize=16)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontweight('bold')

for spine in ax.spines.values():
    spine.set_linewidth(2)
    
#ax.grid(True)

# Set exact plot limits AFTER plotting
ax.set_xlim(vcbmin, vcbmax)
ax.set_ylim(vubmin, vubmax)

# Remove margins around data
ax.margins(0)

# Add legend for filled areas
legend_elements = [
    Patch(facecolor='none', edgecolor='none', label='Measurements'),  # first section title
    Patch(facecolor=colour4, edgecolor=colour4, label='68% region'),
    Patch(facecolor=colour3, edgecolor=colour3, label='95% region'),
    Patch(facecolor='none', edgecolor='none', label='Full Fit'),  # second section title
    Patch(facecolor='darkgrey', edgecolor='black', label='68% region'),
    Patch(facecolor='lightgrey', edgecolor='darkgrey', label='95% region')
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=18,
          handleheight=0.5, handlelength=1.0)

# Adjust layout to minimize whitespace
#plt.subplots_adjust(left=0.05, right=0.95, top=0.85, bottom=0.15)

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

plt.savefig("vubvcb_2025.png", dpi=300, bbox_inches='tight')
plt.savefig("vubvcb_2025.pdf", bbox_inches='tight')

plt.show()
