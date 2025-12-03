#!/usr/bin/env python3

import uproot
import argparse
import sys
import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    from scipy.interpolate import CubicSpline
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# Setup logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s:%(lineno)d - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def print_object(obj):
    for member in obj.all_members.keys():
        print(member, obj.all_members[member])

def get_global_bin_index(bin_indices, thn_axis_nbins):
    """
    Calculates the global bin index from a list of 1-based bin indices for each axis.
    The calculation method must match the one used in SplineCalculator.cxx
    global_bin_id = sum_{i=0 to N-1} (FindBin(v_i)-1)*stride_i
    where stride_i = product_{j=0 to i-1} nbins_j
    Note: uproot bin indices are 0-based for underflow, 1-based for first bin.
    We use 0-based indices internally for calculation.
    """
    strides = [1] * len(thn_axis_nbins)
    for i in range(1, len(thn_axis_nbins)):
        strides[i] = strides[i-1] * thn_axis_nbins[i-1]

    # The bin_indices are 1-based, so subtract 1 for 0-based calculation
    global_bin = sum(b * s for b, s in zip(bin_indices, strides))
    return global_bin


def main():
    parser = argparse.ArgumentParser(description="Plot splines from a ROOT file in a 2D grid using uproot and matplotlib.")
    parser.add_argument("--file", required=True, help="Input ROOT file with splines.")
    parser.add_argument("--syst", required=True, help="Systematic name (directory in the ROOT file).")
    parser.add_argument("--var1", required=True, help="Name of the first variable for the 2D plot (X-axis of grid).")
    parser.add_argument("--var2", required=True, help="Name of the second variable for the 2D plot (Y-axis of grid).")
    parser.add_argument("--output", default="my_spline_plot.pdf", help="Output file name.")
    parser.add_argument("--slice-bins", nargs='*', help="Values for slicing other dimensions, in order. e.g. --slice-bins E_true=1.5 Gen_vtx_z=100. If not provided, the first bin is taken for other dimensions.")

    args = parser.parse_args()

    if not HAS_SCIPY:
        logger.warning("scipy is not installed. Splines will be linearly interpolated. For cubic splines, please 'pip install scipy'")

    with uproot.open(args.file) as root_file:
        if "BinningTemplate" not in root_file:
            logger.error("Could not find 'BinningTemplate' in the file.")
            sys.exit(1)
        
        binning_template = root_file["BinningTemplate"]
        # print(root_file["BinningTemplate"].all_members.keys())
        
        axis_info = {}
        thn_axis_nbins = []
        for i, axis in enumerate(binning_template.all_members["fAxes"]):
            title = axis.all_members['fTitle']
            axis_info[title] = {
                "index": i,
                "nbins": len(axis.edges()) + 2 - 1, # +2 for underflow and overflow
                "edges": [-np.inf] + list(axis.edges()) + [np.inf]
            }
            thn_axis_nbins.append(len(axis.edges()) + 2 - 1)

        if args.var1 not in axis_info:
            logger.error(f"Variable '{args.var1}' not found in BinningTemplate axes.")
            logger.info("Available variables: %s", ", ".join(axis_info.keys()))
            sys.exit(1)
        if args.var2 not in axis_info:
            logger.error(f"Variable '{args.var2}' not found in BinningTemplate axes.")
            sys.exit(1)

        ax1_idx = axis_info[args.var1]["index"]
        ax2_idx = axis_info[args.var2]["index"]
        # ax1_idx = 0
        # ax2_idx = 1
        nbins1 = axis_info[args.var1]["nbins"]
        bins1 = axis_info[args.var1]["edges"]
        nbins2 = axis_info[args.var2]["nbins"]
        bins2 = axis_info[args.var2]["edges"]

        if args.syst not in root_file:
            logger.error(f"Systematic '{args.syst}' not found in file.")
            logger.info("Available systematics: %s", ", ".join(root_file.keys()))
            sys.exit(1)

        syst_dir = root_file[args.syst]
        
        splines = {}
        g_ymin, g_ymax = 1.0, 1.0
        g_xmin, g_xmax = 0.0, 0.0

        for key in syst_dir.keys():
            obj = syst_dir[key]
            # A bit of a hack to check if it's a spline

            if "fNp" in obj.all_members and "Spline" in obj.all_members["fName"]:
                spline_name = obj.all_members["fTitle"]
                splines[spline_name] = obj
                
                x_vals = [pt.member('fX') for pt in obj.member("fPoly")]
                y_vals = [pt.member('fY') for pt in obj.member("fPoly")]

                if len(splines) == 1:
                    g_xmin = x_vals[0]
                    g_xmax = x_vals[-1]
                
                g_ymin = min(g_ymin, np.min(y_vals))
                g_ymax = max(g_ymax, np.max(y_vals))

        if not splines:
            logger.warning(f"No splines found for systematic '{args.syst}'.")
            sys.exit(0)

        logger.info(f"Found {len(splines)} splines for systematic '{args.syst}'.")
        for spline in splines:
            logger.info(f"Spline found: {spline}")
            
        margin = (g_ymax - g_ymin) * 0.15
        g_ymin -= margin
        g_ymax += margin

        fig, axes = plt.subplots(nrows=nbins2, ncols=nbins1, figsize=(12, 10), sharex=True, sharey=True)
        fig.suptitle(f"Splines for {args.syst}", fontsize=16)

        slice_bin_indices = {}
        if args.slice_bins:
            for s in args.slice_bins:
                var, val_str = s.split('=')
                if var not in axis_info:
                    logger.warning(f"Slicing variable '{var}' not found in axes. Ignoring.")
                    continue
                val = float(val_str)
                axis_idx = axis_info[var]['index']
                axis_edges = axis_info[var]['edges']
                # np.digitize is 1-based, perfect for our use case
                slice_bin_indices[axis_idx] = np.digitize(val, axis_edges)

        for i in range(0, nbins1):
            for j in range(0, nbins2):
                ax = axes[nbins2 - j - 1, i] # matplotlib indexing vs ROOT's canvas
                
                bin_indices = [1] * len(thn_axis_nbins)
                for axis_idx, bin_val in slice_bin_indices.items():
                    bin_indices[axis_idx] = bin_val
                
                bin_indices[ax1_idx] = i
                bin_indices[ax2_idx] = j

                print(i, j, bin_indices, thn_axis_nbins)

                global_bin = get_global_bin_index(bin_indices, thn_axis_nbins)
                spline_name = f"{args.syst}_bin{global_bin}"

                if spline_name in splines:
                    spline_obj = splines[spline_name]
                    x_pts = [pt.member('fX') for pt in spline_obj.member("fPoly")]
                    y_pts = [pt.member('fY') for pt in spline_obj.member("fPoly")]

                    if HAS_SCIPY:
                        spline_func = CubicSpline(x_pts, y_pts)
                        x_fine = np.linspace(g_xmin, g_xmax, 100)
                        y_fine = spline_func(x_fine)
                    else: # Fallback to linear interpolation
                        x_fine = x_pts
                        y_fine = y_pts

                    ax.plot(x_fine, y_fine, 'r-')
                    ax.axhline(1.0, color='gray', linestyle='--', linewidth=1)
                else:
                    logger.warning(f"Spline '{spline_name}' not found. Skipping plot for bin ({i}, {j}).")

                ax.set_ylim(g_ymin, g_ymax)
                ax.set_xlim(g_xmin, g_xmax)

                if j == 0:
                    ax.set_xlabel(f"{bins1[i]:.2f}:{bins1[i+1]:.2f}", fontsize=14)
                if i == 0:
                    ax.set_ylabel(f"{bins2[j]:.2f}:{bins2[j+1]:.2f}", fontsize=14)

                ax.tick_params(axis='x', which='both', bottom=False, top=True, labelbottom=False)
                ax.tick_params(axis='y', which='both', left=False, right=True, labelleft=False)

                if i == nbins2 - 1:
                    ax.tick_params(axis='y', which='both', labelright=True, labelsize=14)
                if j == nbins1 - 1:
                    ax.tick_params(axis='x', which='both', labeltop=True, labelsize=14, rotation=90)
                ax.grid(True, linestyle=':', linewidth=0.5)


    # Common axis labels
    fig.text(0.5, 0.04, args.var1, ha='center', va='center', fontsize=16)
    fig.text(0.06, 0.5, args.var2, ha='center', va='center', rotation='vertical', fontsize=16)

    fig.subplots_adjust(hspace=0, wspace=0)

    # plt.tight_layout(rect=[0.08, 0.05, 0.98, 0.95]) # Adjust layout to make room for titles
    
    plt.savefig(args.output)

    logger.info(f"Plot saved to {args.output}")

if __name__ == "__main__":
    main()