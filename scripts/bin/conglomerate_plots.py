import copy
import sys
import shutil
import json

import ROOT
import numpy

import utilities.utilities as utilities
import utilities.root_style as root_style
import conglomerate
from conglomerate.compare_config import CompareConfig
from conglomerate.conglomerate_merge import ConglomerateMerge
from conglomerate.conglomerate_one import ConglomerateOne
from conglomerate.conglomerate_one import ConglomerateContainer
from conglomerate.merge_cuts_summary_tex import MergeCutsSummaryTex

def mc_mod(file_name, axis_range, right_labels, top_labels, verticals = None):
    modifiers = {
        "merge_options":{
            "right_labels":right_labels,
            "top_labels":top_labels
        },
        "redraw":{
            "draw_option":[""],
            "fill_color":[ROOT.kOrange-2],
            "x_range":axis_range,
        },
        "normalise_hist":True,
        "defit":True,
        "write_plots":{"file_name":file_name},
    }
    if verticals != None:
        modifiers = vertical(verticals, modifiers)
    return modifiers

EMIT_COLORS = [ROOT.TColor(10000, 255, 255, 255),
               ROOT.TColor(10001, 245, 245, 255),
               ROOT.TColor(10002, 235, 235, 255)]
def emit_colors():
    # must be the same as EMIT_COLORS declaration
    # risk of collision with existing colors in the palette (need ROOT>6 to
    # avoid) so set index to high number and pray
    return [0, 0, 0]

class CompareMCConfig(CompareConfig):
    def __init__(self, beam, target_dir, top_labels, right_labels):
        dir_list = [
            #target_dir+"plots_"+beam+"/",
            target_dir+"plots_Simulated_"+beam
        ]
        self.setup(beam, target_dir, "mc_plots/", "compare_mc/", dir_list)

        self.conglomerate_list = [
            self.get_conglomerate_3("mc_residual_tof01_t", "t", "Res(t_{01}) [ns]",
                                    None, modifiers = mc_mod("tof01_t", [-1., 1.], 
                                    right_labels, top_labels, [-0.2, 0, 0.2])),
            self.get_conglomerate_3("mc_residual_tku_tp_x", "x", "TKU Res(x) [mm]",
                                    None, modifiers = mc_mod("tku_x", [-2., 2.], 
                                    right_labels, top_labels, [-1, 0, 1])),
            self.get_conglomerate_3("mc_residual_tku_tp_px", "px", "TKU Res(p_{x}) [MeV/c]", 
                                    None, modifiers = mc_mod("tku_px", [-5., 5.], 
                                    right_labels, top_labels, [-2, 0, 2])),
            self.get_conglomerate_3("mc_residual_tku_tp_y", "y", "TKU Res(y) [mm]", 
                                    None, modifiers = mc_mod("tku_y", [-5., 5.], 
                                    right_labels, top_labels, [-1, 0, 1])),
            self.get_conglomerate_3("mc_residual_tku_tp_py", "py", "TKU Res(p_{y}) [MeV/c]", 
                                    None, modifiers = mc_mod("tku_py", [-5., 5.], 
                                    right_labels, top_labels, [-2, 0, 2])),
            self.get_conglomerate_3("mc_residual_tku_tp_pz", "pz", "TKU Res(p_{z}) [MeV/c]", 
                                    None, modifiers = mc_mod("tku_pz", [-10., 10.], 
                                    right_labels, top_labels, [-2, 0, 2])),
            self.get_conglomerate_3("mc_residual_tkd_tp_x", "x", "TKD Res(x) [mm]", 
                                    None, modifiers = mc_mod("tkd_x", [-2., 2.], 
                                    right_labels, top_labels, [-1, 0, 1])),
            self.get_conglomerate_3("mc_residual_tkd_tp_px", "px", "TKD Res(p_{x}) [MeV/c]", 
                                    None, modifiers = mc_mod("tkd_px", [-5., 5.], 
                                    right_labels, top_labels, [-2, 0, 2])),
            self.get_conglomerate_3("mc_residual_tkd_tp_y", "y", "TKD Res(y) [mm]", 
                                    None, modifiers = mc_mod("tkd_y", [-5., 5.], 
                                    right_labels, top_labels, [-1, 0, 1])),
            self.get_conglomerate_3("mc_residual_tkd_tp_py", "py", "TKD Res(p_{y}) [MeV/c]", 
                                    None, modifiers = mc_mod("tkd_py", [-5., 5.], 
                                    right_labels, top_labels, [-2, 0, 2])),
            self.get_conglomerate_3("mc_residual_tkd_tp_pz", "pz", "TKD Res(p_{z}) [MeV/c]", 
                                    None, modifiers = mc_mod("tkd_pz", [-10., 10.], 
                                    right_labels, top_labels, [-2, 0, 2])),
        ]

def fractional_emittance_mod(fraction, y_range, beam, target_dir, top_labels, right_labels):
    output_dir = "compare_fractional_emittance/"
    modifiers = {
        "merge_options":{
            "right_labels":top_labels, # intentionally reversed
            "top_labels":right_labels,
            "row_fill":emit_colors(),
        },
        "hist_title":"",
        "file_name":"fractional_emittance",
        "canvas_name":"fractional_emittance",
        "normalise_hist":False,
        "histogram_names":[],
        "graph_names":["mg_feps", "reco", "all_mc"], #"Reco "+str(fraction), "MC "+str(fraction)],
        "mice_logo":False,
        "legend":False,
        "calculate_errors":[],
        "redraw":{
            "draw_option":["P"]*2, #E1 PLC
            "fill_color":[1]*2,
            "transparency":None,
            "line_color":[1]*2,
            "marker_style":[1]*2,
            "marker_color":[1]*2,
            "draw_order":[1]*2,
            "x_range":[15010, 18990],
            "y_range":y_range,
            "graph":{
                    "draw_option":["P", "P", "P", "P", "P"], # reco, extrap reco, reco mc, mc truth
                    "marker_style":[1, 20, 1, 21, 25],
                    "marker_color":[1, 1, ROOT.kRed, ROOT.kRed, ROOT.kRed],
                    "line_color":[1, 1, ROOT.kRed, ROOT.kRed, ROOT.kRed],
                    "fill_color":[1, 1, ROOT.kRed, ROOT.kRed, ROOT.kRed],
                    "fill_style":None,
                    "transparency":None,
                    "draw_order":[0, 2, 4, 3, 1],
            },
            "ignore_more_histograms":False,
        },
        "write_plots":{
            "dir":target_dir+"/"+output_dir+"/"+beam,
            "file_name":"fractional_emittance_"+str(int(fraction)).rjust(3, '0')
        },
        "axis_title":{
            "x":"z [mm]",
            "y":"Amplitude [mm]",
        },

    }
    return modifiers

def frac_y_range(fraction, beam):
    print "Finding y_range for", fraction, beam
    if abs(fraction-50.0) < 0.1:
        return [12, 36]
    elif abs(fraction-9.0) < 0.1:
        if "4-140" in beam:
            return [3.1, 5.9]
        elif "6-140" in beam:
            return [5.1, 7.9]
        elif "10-140" in beam:
            return [9.1, 11.9]
    print "Failed to find a y range, applying default", fraction, beam
    return [0., 100.]

class CompareFractionalEmittanceConfig(CompareConfig):
    def __init__(self, beam, target_dir, top_labels, right_labels):
        dir_list = [
            target_dir+"plots_"+beam+"/",
            target_dir+"plots_Simulated_"+beam
        ]
        output_dir = "compare_fractional_emittance/"
        self.setup(beam, target_dir, "fractional_emittance/", output_dir, dir_list)
        self.conglomerate_list = []
        for frac in [9.0]:
            y_range = frac_y_range(frac, beam)
            mods = fractional_emittance_mod(frac, y_range, beam, target_dir, top_labels, right_labels)
            self.conglomerate_list.append(
                self.get_conglomerate_0(modifiers = mods),
            )
        self.data_caption = [[],]


def density_y_range(beam):
    y_range = [1, 109]
    if "4-140" in beam:
        y_range = [1, 109]
    elif "6-140" in beam:
        y_range = [1, 59]
    elif "10-140" in beam:
        y_range = [0, 24]
    return y_range


def density_recon_mod(name, beam, target_dir, top_labels, right_labels):
    output_dir = "compare_density/"
    x_range = [0.01, 1.1]
    modifiers = {
        "merge_options":{
            "right_labels":right_labels,
            "top_labels":top_labels,
            "col_fill":emit_colors(),
            "outer_box":{"fill":True},
            "legend":[
                {"fill_color":ROOT.kOrange+4, "transparency":0.8,
                    "marker_color":None, "marker_style":None,
                    "text":"Upstream"},
                {"fill_color":ROOT.kGreen+3, "transparency":0.8,
                    "marker_color":None, "marker_style":None,
                    "text":"Downstream"},
            ]
        },
        "hist_title":"",
        "file_name":name,
        "canvas_name":name,
        "normalise_hist":False,
        "histogram_names":[name],
        "graph_names":[name, "Graph"],
        "mice_logo":False,
        "legend":False,
        "calculate_errors":[],
        "rescale_x":x_range,
        "redraw":{ 
            "draw_option":["p"],
            "fill_color":[0],
            "transparency":None,
            "line_color":[0],
            "marker_style":None,
            "marker_color":[0],
            "draw_order":[0],
            "x_range":x_range,
            "y_range":[0.01, 1.4], #density_y_range(beam),
            "graph":{
                    "draw_option":["L 3", "L 3", "L 3", "L 3", "L 3"], # reco, extrap reco, reco mc, mc truth
                    "marker_style":None,
                    "marker_color":None,
                    "fill_color":[ROOT.kOrange+4, ROOT.kOrange+4, ROOT.kOrange+4, ROOT.kGreen+3, ROOT.kGreen+3],
                    "line_color":[ROOT.kOrange+4, ROOT.kOrange+4, ROOT.kOrange+4, ROOT.kGreen+3, ROOT.kGreen+3],
                    "line_width":[2, 2, 2, 2, 2],
                    "fill_style":[1001, 1001, 1001, 1001, 1001],
                    "transparency":[0, 0.5, 0.5, 0.5, 0.5],
                    "draw_order":[0, 1, 2, 3, 4],
            },
            "ignore_more_histograms":False,
        },
        "write_plots":{
            "dir":target_dir+"/"+output_dir+"/"+beam,
            "file_name":name,
        },
        "defit":False,
        "axis_title":{
            "x":"Fraction of upstream sample",
            "y":"Normalised Density",
        },
    }
    return modifiers

def density_ratio_mod(name, beam, target_dir, top_labels, right_labels):
    output_dir = "compare_density/"
    x_range = [0.01, 1.1]
    modifiers = {
        "merge_options":{
            "right_labels":right_labels,
            "top_labels":top_labels,
        },
        "hist_title":"",
        "file_name":name,
        "canvas_name":name,
        "normalise_hist":False,
        "histogram_names":[""],
        "graph_names":["", "sys", "stats"],
        "mice_logo":False,
        "legend":False,
        "calculate_errors":[],
        "extra_lines":{
            "verticals":[],
            "horizontals":[{
                "y_value":1.0,
                "line_color":1,
                "line_style":2,
                "line_width":2,
            }]
        },
        "rescale_x":x_range,
        "redraw":{ 
            "draw_option":["p", "p"],
            "fill_color":[0, 0],
            "transparency":None,
            "line_color":[0, 0],
            "marker_style":None,
            "marker_color":[0, 0],
            "draw_order":[0, 1],
            "x_range":x_range,
            "y_range":[0.01, 1.4],
            "graph":{
                    "draw_option":["L 3", "L 3", "L 3", "L 3", "L 3", "L 3"],
                    "marker_style":None,
                    "marker_color":None,
                    "fill_color":[ROOT.kBlue, ROOT.kBlue, ROOT.kBlue, ROOT.kRed, ROOT.kRed, ROOT.kRed],
                    "line_color":[ROOT.kBlue, ROOT.kBlue, ROOT.kBlue, ROOT.kRed, ROOT.kRed, ROOT.kRed],
                    "line_width":[2, 2, 2, 2, 2, 2],
                    "fill_style":[1001, 1001, 1001, 1001, 1001, 1001],
                    "transparency":[0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
                    "draw_order":[3, 4, 5, 0, 1, 2, ],
            },
            "ignore_more_histograms":False,
        },
        "write_plots":{
            "dir":target_dir+"/"+output_dir+"/"+beam,
            "file_name":name,
        },
        "defit":False,
        "axis_title":{
            "x":"Fraction of upstream sample",
            "y":"Ratio of TKD to TKU Density",
        },
    }
    return modifiers

class CompareDensityConfig(CompareConfig):
    def __init__(self, beam, target_dir, top_labels, right_labels):
        dir_list = [
            target_dir+"plots_"+beam+"/",
            #target_dir+"plots_Simulated_"+beam
        ]
        output_dir = "compare_density/"
        self.setup(beam, target_dir, "density/", output_dir, dir_list)
        self.conglomerate_list = []
        mods = density_recon_mod("normalised_density_profile_reco", beam, target_dir, top_labels, right_labels)
        self.conglomerate_list.append(
            self.get_conglomerate_0(modifiers = mods),
        )

class CompareDensityRatioConfig(CompareConfig):
    def __init__(self, beam, target_dir, top_labels, right_labels):
        dir_list = [
            target_dir+"plots_"+beam+"/",
            target_dir+"plots_Simulated_"+beam
        ]
        output_dir = "compare_density/"
        self.setup(beam, target_dir, "density/", output_dir, dir_list)
        self.conglomerate_list = []
        mods = density_ratio_mod("density_ratio_reco", beam, target_dir, top_labels, right_labels)
        #mods = self.density_ratio_lines(beam, mods)
        self.conglomerate_list.append(
            self.get_conglomerate_0(modifiers = mods),
        )
        self.data_caption = [[],]


    def density_ratio_lines(self, beam, mods):
        values = {
            "2017-2.7_4-140_None":0.9857,#   0.0103  0.04609  0.04722
            "2017-2.7_6-140_None":0.9809,#  0.01253  0.04618  0.04786
            "2017-2.7_10-140_None":1.063,#  0.03582  0.05633  0.06676
            "2017-2.7_4-140_lH2_empty":0.9831,#  0.01253  0.04518  0.04689
            "2017-2.7_6-140_lH2_empty":1.015,#  0.01847  0.04866  0.05204
            "2017-2.7_10-140_lH2_empty":1.053,#  0.06238  0.05356  0.08222
            "2017-2.7_4-140_lH2_full":1.089,#  0.02459   0.0487  0.05456
            "2017-2.7_6-140_lH2_full":1.174,#  0.01624  0.05405  0.05644
            "2017-2.7_10-140_lH2_full":1.229,#  0.04927  0.06481  0.08141
            "2017-2.7_4-140_LiH":1.076,#  0.01489  0.05012  0.05228
            "2017-2.7_6-140_LiH":1.173,#  0.01685  0.05696   0.0594
            "2017-2.7_10-140_LiH":1.266,#  0.04694  0.06698  0.08179
        }
        mods["extra_lines"]["horizontals"].append({
            "y_value":values[beam],
            "line_color":ROOT.kGray+2,
            "line_style":1,
            "line_width":2,
        })
        return mods



def vertical(x_value_list, modifier):
    verticals = [
        {"x_value":x_value, "line_color":1, "line_style":2, "line_width":2} for x_value in x_value_list
    ]
    modifier = copy.deepcopy(modifier)
    modifier["extra_lines"] = {
            "verticals":verticals,
            "horizontals":[],
        }
    return modifier

class CompareCutsConfig(CompareConfig):
    def __init__(self, beam, target_dir, top_labels, right_labels):
        dir_list = [
            target_dir+"plots_"+beam+"/",
            target_dir+"plots_Simulated_"+beam+"/",
        ]
        self.setup(beam, target_dir, "cut_plots/", "compare_cuts/", dir_list)
        mod = {
            "merge_options":{
                "right_labels":right_labels,
                "top_labels":top_labels
            },
            "mice_logo":False,
            "legend":False,
            "redraw":{
                    "draw_option":["P"],
                    "fill_color":[1],
                    "line_color":[1],
            }
        }
        nu = 9
        nd = 3+nu
        nu1 = str(nu+1)
        nd1 = str(nd+1)
        ndm1 = str(nd-1)
        nu = str(nu)
        nd = str(nd)
        self.conglomerate_list = [
            self.get_conglomerate_2("global_through_virtual_diffuser_us_r_"+nu+"_0", None, "Radius at diffuser (upstream) [mm]", None, True, [0.5, 0.5, 0.9, 0.9], vertical([100], mod)),
            self.get_conglomerate_2("global_through_virtual_diffuser_ds_r_"+nu+"_0", None, "Radius at diffuser (downstream) [mm]", None, True, [0.5, 0.5, 0.9, 0.9], vertical([100], mod)),
            self.get_conglomerate_2("tkd_n_tracks_"+nd+"_0", None, "Number of tracks in TKD", None, True, [0.5, 0.5, 0.9, 0.9], vertical([0.5, 1.5], mod)),
            self.get_conglomerate_2("tkd_chi2_"+nd+"_0", None, "#chi^{2}/D.o.F. in TKD", None, True, [0.5, 0.5, 0.9, 0.9], vertical([8], mod)),
            self.get_conglomerate_2("tkd_max_r_"+nd+"_0", None, "Maximum radius in TKD stations [mm]", [0., 300.], True, [0.5, 0.5, 0.9, 0.9], vertical([150], mod)),
            self.get_conglomerate_2("tkd_p_"+nd1+"_0", None, "Momentum at TKD Reference Plane [MeV/c]", [50., 250.], True, [0.5, 0.5, 0.9, 0.9], vertical([90, 170], mod)),
            self.get_conglomerate_2("tku_chi2_"+nu+"_0", None, "#chi^{2}/D.o.F. in TKU", None, True, [0.5, 0.5, 0.9, 0.9], vertical([8], mod)),
            self.get_conglomerate_2("tku_max_r_"+nu+"_0", None, "Maximum radius in TKU stations [mm]", None, True, [0.5, 0.5, 0.9, 0.9], vertical([150], mod)),
            self.get_conglomerate_2("tku_n_tracks_"+nu+"_0", None, "Number of tracks in TKU", None, True, [0.5, 0.5, 0.9, 0.9], vertical([0.5, 1.5], mod)), # disable two cuts
            self.get_conglomerate_2("tku_p_"+nu+"_0", None, "Momentum at TKU Reference Plane [MeV/c]", [50., 250.], [135, 145], [0.5, 0.5, 0.9, 0.9], vertical([135, 145], mod)),
            self.get_conglomerate_2("tof_tof0_n_sp_"+nu+"_0", None, "Number of space points in ToF0", None, True, [0.5, 0.5, 0.9, 0.9], vertical([0.5, 1.5], mod)),
            self.get_conglomerate_2("tof_tof1_n_sp_"+nu+"_0", None, "Number of space points in ToF1", None, True, [0.5, 0.5, 0.9, 0.9], vertical([0.5, 1.5], mod)),
            self.get_conglomerate_2("tof_tof01_"+nu+"_0", None, "Time between ToF0 and ToF1 [ns]", [2.5, 5.], True, [0.5, 0.5, 0.9, 0.9], vertical([3.0, 3.5, 4.0, 4.5], mod)),
            #self.get_conglomerate_2("tof_delta_tof01_"+nu1+"_0", None, "t(ToF01) - extrapolated t(ToF01) [ns]", None, True, [0.1, 0.5, 0.5, 0.9], mod),

            self.get_conglomerate_2("tku_scifi_n_planes_with_clusters_"+nu1+"_0", None, "Number of planes with clusters in TKU", None, True, [0.1, 0.5, 0.5, 0.9], modifiers = mod),
            self.get_conglomerate_2("tku_scifi_n_planes_with_clusters_"+nu+"_0",  None, "Number of planes with clusters in TKU", None, True, [0.1, 0.5, 0.5, 0.9], modifiers = mod),
            self.get_conglomerate_2("tku_scifi_n_planes_with_clusters_"+nu+"_1",  None, "Number of planes with clusters in TKU", None, True, [0.1, 0.5, 0.5, 0.9], modifiers = mod),
            self.get_conglomerate_2("tkd_scifi_n_planes_with_clusters_"+nu1+"_1", None, "Number of planes with clusters in TKD", None, True, [0.1, 0.5, 0.5, 0.9], modifiers = mod),
            self.get_conglomerate_2("tkd_scifi_n_planes_with_clusters_"+ndm1+"_1", None, "Number of planes with clusters in TKD", None, True, [0.1, 0.5, 0.5, 0.9], modifiers = mod),
            self.get_conglomerate_2("tkd_scifi_n_planes_with_clusters_"+nd1+"_0", None, "Number of planes with clusters in TKD", None, True, [0.1, 0.5, 0.5, 0.9], modifiers = mod),
        ]
        self.data_caption = [[],]
        self.data_caption[0] = [
            " Samples are listed for None and lH2 empty datasets.",
            " Samples are listed for lH2 full and LiH datasets.",
        ]
        self.data_caption.append(copy.deepcopy(self.data_caption[0]))
        self.data_caption.append(copy.deepcopy(self.data_caption[0]))
        self.mc_caption = copy.deepcopy(self.data_caption)
        self.mc_caption.append(copy.deepcopy(self.mc_caption[0]))
        self.mc_caption.append(copy.deepcopy(self.mc_caption[0]))
        data_prefix = ["upstream", "downstream", "extrapolated"]
        mc_prefix = ["upstream reconstructed", "downstream reconstructed", "extrapolated reconstructed",
                     "upstream truth", "downstream truth"]
        for i, prefix in enumerate(data_prefix):
            for j, item in enumerate(self.data_caption[i]):
                self.data_caption[i][j] = "The "+prefix+" reconstructed data sample is listed. "+item
        for i, prefix in enumerate(mc_prefix):
            for j, item in enumerate(self.mc_caption[i]):
                self.mc_caption[i][j] = "The "+prefix+" simulated sample is listed. "+item


class CompareData1DConfig(CompareConfig):
    def __init__(self, beam, target_dir, top_labels, right_labels):
        dir_list = [
            target_dir+"plots_"+beam+"/",
            target_dir+"plots_Simulated_"+beam
        ]
        self.setup(beam, target_dir, "data_plots/", "compare_data/", dir_list)
        modifiers = {
            "merge_options":{
                "right_labels":right_labels,
                "top_labels":top_labels
            }
        }
      
        self.conglomerate_list = [
            self.get_conglomerate_1("tof0_slab_dt", "tof0_slab_dt", "us_cut", "Slab dt for ToF0 [ns]", [-2, 2], False, [0.55, 0.7, 0.9, 0.9], modifiers = modifiers),
            self.get_conglomerate_1("tof1_slab_dt", "tof1_slab_dt", "us_cut", "Slab dt for ToF1 [ns]", [-2, 2], False, [0.55, 0.7, 0.9, 0.9], modifiers = modifiers),
            self.get_conglomerate_1("tof2_slab_dt", "tof2_slab_dt", "ds_cut", "Slab dt for ToF2 [ns]", [-2, 2], False, [0.55, 0.7, 0.9, 0.9], modifiers = modifiers),
            self.get_conglomerate_1("tku_x",  "tku_x",  "us_cut", "x at TKU Reference Plane [mm]",        [-175, 300], False, [0.55, 0.7, 0.9, 0.9], modifiers = modifiers),
            self.get_conglomerate_1("tku_y",  "tku_y",  "us_cut", "y at TKU Reference Plane [mm]",        [-175, 300], False, [0.55, 0.7, 0.9, 0.9], modifiers = modifiers),
            self.get_conglomerate_1("tku_px", "tku_px", "us_cut", "p_{x} at TKU Reference Plane [MeV/c]", [-120, 200], False, [0.55, 0.7, 0.9, 0.9], rescale_y = [0., 0.15], modifiers = modifiers),
            self.get_conglomerate_1("tku_py", "tku_py", "us_cut", "p_{y} at TKU Reference Plane [MeV/c]", [-120, 200], False, [0.55, 0.7, 0.9, 0.9], modifiers = modifiers),
            self.get_conglomerate_1("tku_p", "tku_p", "us_cut", "p at TKU Reference Plane [MeV/c]",       [89, 200], False, [0.55, 0.7, 0.9, 0.9], modifiers = modifiers),
            self.get_conglomerate_1("tkd_x",  "tkd_x",  "ds_cut", "x at TKD Reference Plane [mm]",        [-175, 300], False, [0.55, 0.7, 0.9, 0.9], modifiers = modifiers),
            self.get_conglomerate_1("tkd_y",  "tkd_y",  "ds_cut", "y at TKD Reference Plane [mm]",        [-175, 300], False, [0.55, 0.7, 0.9, 0.9], modifiers = modifiers),
            self.get_conglomerate_1("tkd_px", "tkd_px", "ds_cut", "p_{x} at TKD Reference Plane [MeV/c]", [-120, 200], False, [0.55, 0.7, 0.9, 0.9], modifiers = modifiers),
            self.get_conglomerate_1("tkd_py", "tkd_py", "ds_cut", "p_{y} at TKD Reference Plane [MeV/c]", [-120, 200], False, [0.55, 0.7, 0.9, 0.9], modifiers = modifiers),
            self.get_conglomerate_1("tkd_p", "tkd_p", "ds_cut", "p at TKD Reference Plane [MeV/c]",       [89, 200], False, [0.55, 0.7, 0.9, 0.9], modifiers = modifiers),
            self.get_conglomerate_1("p_res", "p_res", "ds_cut", "p(TKU) - p(TKD) [MeV/c]",       [-20., 40.], True, [0.55, 0.7, 0.9, 0.9], modifiers = vertical([0., 5., 10., 15., 20.], modifiers)),
        ]


class CompareData2DConfig(CompareConfig):
    def __init__(self, beam, target_dir, top_labels, right_labels):
        dir_list = [
            target_dir+"plots_"+beam+"/",
            #target_dir+"plots_Simulated_"+beam
        ]
        self.setup(beam, target_dir, "data_plots/", "compare_data_2d/", dir_list)
        mod = {
            "merge_options":{
                "right_labels":right_labels,
                "top_labels":top_labels
            },
            "redraw":{
                "draw_option":["COL"],
            },
            "canvas_fill_color":root_style.get_frame_fill(),
            "extra_lines":{
                "verticals":[
                  {"x_value":x_value, "line_color":ROOT.kRed+2, "line_style":2, "line_width":2} for x_value in [0, 20]
                ],
                "horizontals":[],
            }
        }

        self.conglomerate_list = [
            self.get_conglomerate_3("p_res_vs_global_through_virtual_absorber_centre_r_ds_cut", "p_res_vs_global_through_virtual_absorber_centre_r_ds_cut", "P(TKU) - P(TKD) [MeV/c]", "radius [mm]", modifiers = mod),
            self.get_conglomerate_3("p_res_vs_global_through_virtual_absorber_centre_y_ds_cut", "p_res_vs_global_through_virtual_absorber_centre_y_ds_cut", "P(TKU) - P(TKD) [MeV/c]", "y [mm]", modifiers = mod),
        ]
        mod = {
            "merge_options":{
                "right_labels":right_labels,
                "top_labels":top_labels
            },
            "redraw":{
                "draw_option":["COL"],
                "graph":{
                    "marker_style":None,
                    "marker_color":None,
                    "draw_option":["l"]*3,
                    "transparency":None,
                    "draw_order":None,
                    "fill_color":None,
                    "fill_style":None,
                }
            },
            "rescale_x":[-1.5, 9.5],
            "rescale_y":[81., 181.],
            "graph_names":["tramlines_upper", "tramlines_lower", "cuts_graph"],
            "canvas_fill_color":root_style.get_frame_fill(),
        }
        self.conglomerate_list += [
            self.get_conglomerate_3("p_tku_vs_tof01_all", "p_tot_vs_tof", "ToF_{01} [ns]", "P(TKU) [MeV/c]", modifiers = mod),
            self.get_conglomerate_3("p_tku_vs_tof01_us_cut", "p_tot_vs_tof", "ToF_{01} [ns]", "P(TKU) [MeV/c]", modifiers = mod),
        ]



class CompareData2DMCConfig(CompareConfig):
    def __init__(self, beam, target_dir, top_labels, right_labels):
        dir_list = [
            #target_dir+"plots_"+beam+"/",
            target_dir+"plots_Simulated_"+beam
        ]
        self.setup(beam, target_dir, "data_plots/", "compare_data_2d_mc/", dir_list)
        mod = {
            "merge_options":{
                "right_labels":right_labels,
                "top_labels":top_labels
            },
            "redraw":{
                "draw_option":["COL"],
            },
            "canvas_fill_color":root_style.get_frame_fill(),
            "extra_lines":{
                "verticals":[
                  {"x_value":x_value, "line_color":ROOT.kRed+2, "line_style":2, "line_width":2} for x_value in [0, 20]
                ],
                "horizontals":[],
            }
        }

        self.conglomerate_list = [
            self.get_conglomerate_3("p_res_vs_global_through_virtual_absorber_centre_r_ds_cut", "p_res_vs_global_through_virtual_absorber_centre_r_ds_cut", "P(TKU) - P(TKD) [MeV/c]", "radius [mm]", modifiers = mod),
            self.get_conglomerate_3("p_res_vs_global_through_virtual_absorber_centre_y_ds_cut", "p_res_vs_global_through_virtual_absorber_centre_y_ds_cut", "P(TKU) - P(TKD) [MeV/c]", "y [mm]", modifiers = mod),
        ]
        mod = {
            "merge_options":{
                "right_labels":right_labels,
                "top_labels":top_labels
            },
            "redraw":{
                "draw_option":["COL"],
                "graph":{
                    "marker_style":None,
                    "marker_color":None,
                    "draw_option":["l"]*3,
                    "transparency":None,
                    "draw_order":None,
                    "fill_color":None,
                    "fill_style":None,
                }
            },
            "rescale_x":[-1.5, 9.5],
            "rescale_y":[81., 181.],
            "graph_names":["tramlines_upper", "tramlines_lower", "cuts_graph"],
            "canvas_fill_color":root_style.get_frame_fill(),
        }
        self.conglomerate_list += [
            self.get_conglomerate_3("p_tku_vs_tof01_all", "p_tot_vs_tof", "ToF_{01} [ns]", "P(TKU) [MeV/c]", modifiers = mod),
            self.get_conglomerate_3("p_tku_vs_tof01_us_cut", "p_tot_vs_tof", "ToF_{01} [ns]", "P(TKU) [MeV/c]", modifiers = mod),
        ]


class CompareOpticsConfig(CompareConfig):
    def __init__(self, beam, target_dir, top_labels, right_labels):
        dir_list = [
            target_dir+"plots_"+beam+"/",
            #target_dir+"plots_Simulated_"+beam
        ]
        self.setup(beam, target_dir, "optics_plots/", "compare_optics/", dir_list)
        graphs = ["_source_tku", "_source_tkd", "_virtual_absorber_centre",]
        for tof in ["_tof0", "_tof1", "_tof2"]:
            for element in "_us", "", "_ds":
                graphs.append(tof+element)
        for tracker in ["_tku", "_tkd"]:
            for station in range(2, 6)+["tp"]:
                graphs.append(tracker+"_"+str(station))
        modifiers = {
            "merge_options":{
                "right_labels":right_labels,
                "top_labels":top_labels
            },
            "redraw":{
                "graph":{
                    "marker_style":None,
                    "marker_color":None,
                    "draw_option":["same lp"]*len(graphs),
                    "draw_order":None,
                    "fill_color":None,
                }
            },
        }

        self.conglomerate_list = [
            self.get_conglomerate_graph("beta_4d_ds", "z [m]", "#beta_{4D} [mm]", graph_list = ["beta_4d_ds"+name for name in graphs], x_range = [12.9, 21.2], y_range = [0., 1200.0], modifiers = modifiers),
            self.get_conglomerate_graph("beta_x_ds",  "z [m]", "#beta_{x} [mm]", graph_list = ["beta_x_ds"+name for name in graphs], x_range = [12.9, 21.2], y_range = [0., 1200.0], modifiers = modifiers),
            self.get_conglomerate_graph("beta_y_ds",  "z [m]", "#beta_{y} [mm]", graph_list = ["beta_y_ds"+name for name in graphs], x_range = [12.9, 21.2], y_range = [0., 1200.0], modifiers = modifiers),
            self.get_conglomerate_graph("sigma_0_ds", "z [m]", "#sigma_{x} [mm]", graph_list = ["sigma_0_ds"+name for name in graphs], x_range = [12.9, 21.2], y_range = [20., 80.0], modifiers = modifiers),
            self.get_conglomerate_graph("sigma_2_ds", "z [m]", "#sigma_{y} [mm]", graph_list = ["sigma_2_ds"+name for name in graphs], x_range = [12.9, 21.2], y_range = [20., 80.0],  modifiers = modifiers),
        ]
 
class CompareOpticsMCConfig(CompareConfig):
    def __init__(self, beam, target_dir, top_labels, right_labels):
        dir_list = [
            #target_dir+"plots_"+beam+"/",
            target_dir+"plots_Simulated_"+beam
        ]
        self.setup(beam, target_dir, "optics_plots/", "compare_optics_mc/", dir_list)
        graphs = ["_source_tku", "_source_tkd", "_virtual_absorber_centre",]
        for tof in ["_tof0", "_tof1", "_tof2"]:
            for element in "_us", "", "_ds":
                graphs.append(tof+element)
        for tracker in ["_tku", "_tkd"]:
            for station in range(2, 6)+["tp"]:
                graphs.append(tracker+"_"+str(station))
        modifiers = {
            "merge_options":{
                "right_labels":right_labels,
                "top_labels":top_labels
            },
            "redraw":{
                "graph":{
                    "marker_style":None,
                    "marker_color":None,
                    "draw_option":["same lp"]*len(graphs),
                    "draw_order":None,
                    "fill_color":None,
                }
            },
        }

        self.conglomerate_list = [
            self.get_conglomerate_graph("beta_4d_ds", "z [m]", "#beta_{4D} [mm]", graph_list = ["beta_4d_ds"+name for name in graphs], x_range = [12.9, 21.2], y_range = [0., 1200.0], modifiers = modifiers),
            self.get_conglomerate_graph("beta_x_ds",  "z [m]", "#beta_{x} [mm]", graph_list = ["beta_x_ds"+name for name in graphs], x_range = [12.9, 21.2], y_range = [0., 1200.0], modifiers = modifiers),
            self.get_conglomerate_graph("beta_y_ds",  "z [m]", "#beta_{y} [mm]", graph_list = ["beta_y_ds"+name for name in graphs], x_range = [12.9, 21.2], y_range = [0., 1200.0], modifiers = modifiers),
            self.get_conglomerate_graph("sigma_0_ds", "z [m]", "#sigma_{x} [mm]", graph_list = ["sigma_0_ds"+name for name in graphs], x_range = [12.9, 21.2], y_range = [20., 80.0], modifiers = modifiers),
            self.get_conglomerate_graph("sigma_2_ds", "z [m]", "#sigma_{y} [mm]", graph_list = ["sigma_2_ds"+name for name in graphs], x_range = [12.9, 21.2], y_range = [20., 80.0],  modifiers = modifiers),
        ]


class CompareGlobalsConfig(CompareConfig):
    def __init__(self, beam, target_dir, top_labels, right_labels):
        dir_list = [
            target_dir+"plots_"+beam+"/",
            target_dir+"plots_Simulated_"+beam
        ]
        self.setup(beam, target_dir, "global_plots/", "compare_globals/", dir_list)
        mod1 = {
            "legend":False,
            "rebin":10,
            "write_fit":False,
            "merge_options":{
                "right_labels":right_labels,
                "top_labels":top_labels
            },
            "extra_lines":{
                    "verticals":[{"x_value":0., "line_color":1, "line_style":2, "line_width":2}],
                    "horizontals":[],
            },
        }
        mod2 = copy.deepcopy(mod1)
        mod2["rebin"] = 10
        self.conglomerate_list = [
            self.get_conglomerate_2("global_through_residual_tof0_t", "res_ex_cut", "Residual t in ToF0 [ns]", [-2, 2], False, [0.5, 0.5, 0.9, 0.9], modifiers = mod1),
            self.get_conglomerate_2("global_through_residual_tkd_tp_p", "res_ex_cut", "Residual Momentum in TKD [MeV/c]", [-20, 20], False, [0.5, 0.5, 0.9, 0.9], modifiers = mod2),
            self.get_conglomerate_2("global_through_residual_tkd_tp_px", "res_ex_cut", "Residual P_{x} in TKD [MeV/c]", None, False, [0.5, 0.5, 0.9, 0.9], modifiers = mod2),
            self.get_conglomerate_2("global_through_residual_tkd_tp_py", "res_ex_cut", "Residual P_{y} in TKD [MeV/c]", None, False, [0.5, 0.5, 0.9, 0.9], modifiers = mod2),
            self.get_conglomerate_2("global_through_residual_tkd_tp_pz", "res_ex_cut", "Residual P_{z} in TKD [MeV/c]", [-20, 20], False, [0.5, 0.5, 0.9, 0.9], modifiers = mod2),
            self.get_conglomerate_2("global_through_residual_tkd_tp_x", "res_ex_cut", "Residual x in TKD [mm]", None, False, [0.5, 0.5, 0.9, 0.9], modifiers = mod2),
            self.get_conglomerate_2("global_through_residual_tkd_tp_y", "res_ex_cut", "Residual y in TKD [mm]", None, False, [0.5, 0.5, 0.9, 0.9], modifiers = mod2),
            self.get_conglomerate_2("global_through_residual_tof1_x", "res_ex_cut", "Residual x in ToF1 [mm]", None, False, [0.5, 0.5, 0.9, 0.9], modifiers = mod2),
            self.get_conglomerate_2("global_through_residual_tof1_y", "res_ex_cut", "Residual y in ToF1 [mm]", None, False, [0.5, 0.5, 0.9, 0.9], modifiers = mod2),
            self.get_conglomerate_2("global_ds_residual_tof2_x", "res_ex_cut", "Residual x in ToF2 [mm]", None, False, [0.5, 0.5, 0.9, 0.9], modifiers = mod2),
            self.get_conglomerate_2("global_ds_residual_tof2_y", "res_ex_cut", "Residual y in ToF2 [mm]", None, False, [0.5, 0.5, 0.9, 0.9], modifiers = mod2),
            self.get_conglomerate_2("global_through_residual_tof2_t", "res_ex_cut", "Residual t in ToF2 [ns]", [-2, 2], False, [0.5, 0.5, 0.9, 0.9], modifiers = mod2),
        ]

def amplitude_x_range(beam):
    x_range = [1., 99.9]
    if "4-140" in beam:
        x_range = [1., 49.9]
    elif "6-140" in beam:
        x_range = [1., 59.9]
    elif "10-140" in beam:
        x_range = [1., 79.9]
    return x_range


class CompareAmplitudeConfigMC(CompareConfig): # MC corrections
    def __init__(self, beam, target_dir, top_labels, right_labels):
        dir_list = [
            target_dir+"plots_Simulated_"+beam
        ]
        self.setup(beam, target_dir, "amplitude/", "compare_amplitude_mc/", dir_list)
        ratio_modifiers = {
            "extra_lines":{
                "horizontals":[
                    {"y_value":1., "line_color":1, "line_style":2, "line_width":2},
                ],
                "verticals":[],
            },
            "merge_options":{
                "right_labels":right_labels,
                "top_labels":top_labels,
                "col_fill":emit_colors(),
            },
            "rescale_x":amplitude_x_range(beam),
        }
        absolute_modifiers = {
            "merge_options":{
                "right_labels":right_labels,
                "top_labels":top_labels,
                "col_fill":emit_colors(),
            },
            "normalise_graph":"Upstream stats",
            "redraw":{
                    "graph":{
                        "draw_option":["p", "3", "p", "3"],
                        "marker_style":[20, 20, 22, 22],
                        "marker_color":[ROOT.kOrange+4, ROOT.kOrange+4, ROOT.kGreen+3, ROOT.kGreen+3],
                        "fill_color":[ROOT.kOrange+4, ROOT.kOrange+4, ROOT.kGreen+3, ROOT.kGreen+3],
                        "fill_style":[1001, 1001, 1001, 1001],
                        "transparency":[0.5, 0.5, 0.5, 0.5],
                        "draw_order":[0, 1, 2, 3],
                    }
            },
            "rescale_x":amplitude_x_range(beam),
        }

        self.conglomerate_list = [
            self.get_conglomerate_graph("amplitude_pdf_all_mc", "Amplitude [mm]",
                                        "Number",
                                        "amplitude_pdf_all_mc", ["Upstream sys hist"],
                                        ["Upstream stats", "Upstream sys", "Downstream stats", "Downstream sys"],
                                        y_range = [0.001, 1.24], x_range = [0.01, 99.9], replace_hist = True,
                                        modifiers = absolute_modifiers),
            self.get_conglomerate_graph("amplitude_pdf_reco", "Amplitude [mm]",
                                        "Normalised Number of Events",
                                        "amplitude_pdf_reco", ["Upstream sys hist"],
                                        ["Upstream stats", "Upstream sys", "Downstream stats", "Downstream sys"],
                                        y_range = [0.001, 1.24], x_range = [0.01, 99.9], replace_hist = True,
                                        modifiers = absolute_modifiers),
            self.get_conglomerate_graph("amplitude_cdf_reco", "Amplitude [mm]",
                                        "Cumulative density",
                                        "amplitude_cdf_reco", ["Upstream CDF stats hist"],
                                        ["Upstream CDF stats", "Upstream CDF sys", "Downstream CDF stats", "Downstream CDF sys"],
                                        y_range = [0.001, 1.099], x_range = [0.01, 99.9], replace_hist = True,
                                        modifiers = absolute_modifiers),

            self.get_conglomerate_graph("pdf_ratio*", "Amplitude [mm]",
                                        "P_{Amp}",
                                        "pdf_ratio", ["PDF Ratio stats hist"],
                                        ["PDF Ratio stats", "PDF Ratio sys"], x_range = [0.01, 99.9], y_range = [0.501, 2.499], replace_hist = True,
                                        graph_draw_option = ["p", "2"], graph_marker_style=[20, 20], graph_marker_color=[ROOT.kRed, ROOT.kRed], graph_draw_order=[1,0], modifiers=ratio_modifiers),
            self.get_conglomerate_graph("cdf_ratio*", "amplitude [mm]",
                                        "R_{Amp}",
                                        "cdf_ratio", ["CDF_Ratio stats hist"],
                                        ["CDF_Ratio stats", "CDF_Ratio sys"], x_range = [0.01, 99.9], y_range = [0.801, 1.399], replace_hist = True,
                                        graph_draw_option = ["p", "2"], graph_marker_style=[20, 20], graph_marker_color=[ROOT.kRed, ROOT.kRed], graph_draw_order=[1,0], modifiers=ratio_modifiers),
        ]

class CompareAmplitudeConfigData(CompareConfig): # data plots
    def __init__(self, beam, target_dir, top_labels, right_labels):
        dir_list = [
            target_dir+"plots_"+beam+"/",
            #target_dir+"plots_Simulated_"+beam
        ]
        self.setup(beam, target_dir, "amplitude/", "compare_amplitude_data/", dir_list)
        ratio_modifiers = {
            "extra_lines":{
                "horizontals":[
                    {"y_value":1., "line_color":1, "line_style":2, "line_width":2},
                ],
                "verticals":[
                    {"x_value":30., "line_color":4, "line_style":2, "line_width":2},
                ],
            },
            "merge_options":{
                "right_labels":right_labels,
                "top_labels":top_labels,
                "col_fill":emit_colors(),
                "outer_box":{"fill":True},
                "legend":[
                    {"fill_color":ROOT.kOrange+4, "transparency":0.5,
                     "marker_color":ROOT.kOrange+4, "marker_style":20,
                     "text":"Upstream"},
                    {"fill_color":ROOT.kGreen+3, "transparency":0.5,
                     "marker_color":ROOT.kGreen+3, "marker_style":22,
                     "text":"Downstream"},
                ]
            },
            "rescale_x":amplitude_x_range(beam),
        }
        absolute_modifiers = {
            "merge_options":{
                "right_labels":right_labels,
                "top_labels":top_labels,
                "outer_box":{"fill":True},
                "legend":[
                    {"fill_color":ROOT.kOrange+4, "transparency":0.5,
                     "marker_color":ROOT.kOrange+4, "marker_style":20,
                     "text":"Upstream"},
                    {"fill_color":ROOT.kGreen+3, "transparency":0.5,
                     "marker_color":ROOT.kGreen+3, "marker_style":22,
                     "text":"Downstream"},
                ]
            },
            "normalise_graph":"Upstream stats",
            "redraw":{
                    "graph":{
                        "draw_option":["p", "3", "p", "3"],
                        "marker_style":[20, 20, 22, 22],
                        "marker_color":[ROOT.kOrange+4, ROOT.kOrange+4, ROOT.kGreen+3, ROOT.kGreen+3],
                        "fill_color":[ROOT.kOrange+4, ROOT.kOrange+4, ROOT.kGreen+3, ROOT.kGreen+3],
                        "fill_style":[1001, 1001, 1001, 1001],
                        "transparency":[0.5, 0.5, 0.5, 0.5],
                        "draw_order":[0, 1, 2, 3],
                    }
            },
            "rescale_x":amplitude_x_range(beam),
            "extra_lines":{
                "horizontals":[],
                "verticals":[
                    {"x_value":30., "line_color":4, "line_style":2, "line_width":2},
                ],
            }
        }

        absolute_modifiers_err = copy.deepcopy(absolute_modifiers)
        absolute_modifiers_err["write_plots"] = {}
        absolute_modifiers_err = {
            "merge_options":{
                "right_labels":right_labels,
                "top_labels":top_labels
            },
            "normalise_graph":True,
            "redraw":{
                    "graph":{
                        "draw_option":["p", "p", "3", "p", "p", "3"],
                        "marker_style":[24, 20, 20, 26, 22, 22],
                        "marker_color":[ROOT.kOrange+4, ROOT.kOrange+4, ROOT.kOrange+4, ROOT.kGreen+3, ROOT.kGreen+3, ROOT.kGreen+3],
                        "fill_color":[ROOT.kOrange+4, ROOT.kOrange+4, ROOT.kOrange+4, ROOT.kGreen+3, ROOT.kGreen+3, ROOT.kGreen+3],
                        "fill_style":[1001, 1001, 1001, 1001, 1001, 1001],
                        "transparency":[0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
                        "draw_order":[1, 2, 4, 5, 0, 3],
                    }
            },
            "rescale_x":amplitude_x_range(beam),
            "write_plots":{
                "file_name":"amplitude_pdf_reco_correction",
                "dir":self.beam_plot_dir,
                "formats":["png", "root", "pdf"],
            }
        }

        self.conglomerate_list = [
            self.get_conglomerate_graph("amplitude_pdf_reco", "Amplitude [mm]",
                                        "Normalised Number of Events",
                                        "amplitude_pdf_reco", ["Upstream sys hist"],
                                        ["Upstream stats", "Upstream sys", "Downstream stats", "Downstream sys"],
                                        x_range = [0.01, 99.9], y_range = [0.001, 1.24], replace_hist = True,
                                        modifiers = absolute_modifiers),
            self.get_conglomerate_graph("amplitude_cdf_reco", "Amplitude [mm]",
                                        "Cumulative density",
                                        "amplitude_cdf_reco", ["Upstream CDF stats hist"],
                                        ["Upstream CDF stats", "Upstream CDF sys", "Downstream CDF stats", "Downstream CDF sys"],
                                        y_range = [0.001, 1.399], x_range = [0.01, 99.9], replace_hist = True,
                                        modifiers = absolute_modifiers),
            self.get_conglomerate_graph("pdf_ratio*", "Amplitude [mm]",
                                        "P_{Amp}",
                                        "pdf_ratio", ["PDF Ratio stats hist"],
                                        ["PDF Ratio stats", "PDF Ratio sys"], x_range = [0.01, 99.9], y_range = [0.501, 2.499], replace_hist = True,
                                        graph_draw_option = ["p", "2"], graph_marker_style=[20, 20], graph_marker_color=[1, 1], graph_draw_order=[1,0], modifiers=ratio_modifiers),
            self.get_conglomerate_graph("cdf_ratio*", "Amplitude [mm]",
                                        "R_{Amp}",
                                        "cdf_ratio", ["CDF_Ratio_stats hist"],
                                        ["CDF_Ratio stats", "CDF_Ratio sys"], x_range = [0.01, 99.9], y_range = [0.601, 1.399], replace_hist = True,
                                        graph_draw_option = ["p", "2"], graph_marker_style=[20, 20], graph_marker_color=[1, 1], graph_draw_order=[1,0], modifiers=ratio_modifiers),
            self.get_conglomerate_graph("amplitude_pdf_reco", "Amplitude [mm]",
                                        "Normalised Number of Events",
                                        "amplitude_pdf_reco", ["Upstream sys hist"],
                                        ["Raw upstream", "Upstream stats", "Upstream sys", "Raw downstream", "Downstream stats", "Downstream sys"],
                                        x_range = [0.01, 99.9], y_range = [0.001, 1.24], replace_hist = True,
                                        modifiers = absolute_modifiers_err),
        ]

class CompareAmplitudeConfigBoth(CompareConfig): # comparisons
    def __init__(self, beam, target_dir, top_labels, right_labels):
        dir_list = [
            target_dir+"plots_"+beam+"/",
            target_dir+"plots_Simulated_"+beam
        ]
        self.setup(beam, target_dir, "amplitude/", "compare_amplitude_both/", dir_list)
        ratio_modifiers = {
            "extra_lines":{
                "horizontals":[
                    {"y_value":1., "line_color":1, "line_style":2, "line_width":2},
                ],
                "verticals":[
                    {"x_value":30., "line_color":4, "line_style":2, "line_width":2},
                ],
            },
            "merge_options":{
                "right_labels":right_labels,
                "top_labels":top_labels,
                #"col_fill":emit_colors(),
                "outer_box":{"fill":True},
                "legend":[
                    {"fill_color":ROOT.kBlue, "transparency":0.5,
                     "marker_color":1, "marker_style":20,
                     "text":"Measured"},
                    {"fill_color":ROOT.kRed, "transparency":0.2,
                     "marker_color":ROOT.kRed, "marker_style":22,
                     "text":"Simulated"},
                ]
            },
            "redraw":{
                "graph":{
                    "fill_style":[1001, 1001, 1001, 1001],
                    "fill_color":[ROOT.kBlue, ROOT.kBlue, ROOT.kRed, ROOT.kRed],
                    "transparency":[0.5, 0.5, 0.2, 0.2],
                }
            },
            "rescale_x":amplitude_x_range(beam),
        }

        self.conglomerate_list = [
            self.get_conglomerate_graph("pdf_ratio*", "Amplitude [mm]",
                                        "#frac{Number out}{Number in}",
                                        "pdf_ratio", ["PDF Ratio stats_hist"],
                                        ["PDF_Ratio stats", "PDF_Ratio sys"], x_range = [0.01, 99.9], y_range = [0.01, 1.49], replace_hist = True,
                                        graph_draw_option = ["p", "3", "p", "3"], graph_marker_style=[20, 20, 22, 22], 
                                        graph_marker_color=[1, 1, ROOT.kRed-7, ROOT.kRed-7], graph_draw_order=[3, 1, 2, 0,], 
                                        modifiers=ratio_modifiers),
            self.get_conglomerate_graph("cdf_ratio*", "Amplitude [mm]",
                                        "R_{Amp}",
                                        "cdf_ratio", ["CDF_Ratio_stats hist"],
                                        ["CDF_Ratio stats", "CDF_Ratio sys"], x_range = [0.01, 99.9], y_range = [0.601, 1.399], replace_hist = True,
                                        graph_draw_option = ["p", "3", "p", "3"], graph_marker_style=[20, 20, 22, 22], 
                                        graph_marker_color=[1, 1, ROOT.kRed-7, ROOT.kRed-7], graph_draw_order=[3, 1, 2, 0,], 
                                        modifiers=ratio_modifiers),
        ]

def cuts_summary(dir_list, target_dir):
    data_cuts_summary = MergeCutsSummaryTex()
    mc_cuts_summary = MergeCutsSummaryTex()
    for beam in dir_list:
        config = CompareCutsConfig(beam, target_dir, [], [])
        data_cuts_summary.append_summary(config, [0])
        mc_cuts_summary.append_summary(config, [1])
    data_cuts_summary.caption = config.data_caption
    mc_cuts_summary.caption = config.mc_caption
    data_cuts_summary.merge_summaries(target_dir+"/cuts_summary/data/", "data_cuts_summary")
    mc_cuts_summary.merge_summaries(target_dir+"/cuts_summary/mc/", "mc_cuts_summary")



def run_conglomerate(batch_level, config_list, dir_lists, do_cuts_summary, target_dir, top_labels, right_labels):
    rows = len(dir_lists)
    cols = min([len(sub_list) for sub_list in dir_lists])
    dir_list = []
    for sub_list in dir_lists:
        dir_list += sub_list
    if do_cuts_summary:
        cuts_summary(dir_list, target_dir)
    fail_dict = {}
    for ConfigClass in config_list:
        fail_list = []
        conglomerate_list = []
        print "Doing", ConfigClass.__name__
        for beam in dir_list:
            #if batch_level > 5:
                #ROOT.gROOT.SetBatch(False)
            try:
                config = ConfigClass(beam, target_dir, top_labels, right_labels)
                cong = ConglomerateContainer(config)
                cong.conglomerate()
                conglomerate_list.append(cong)
            except Exception:
                sys.excepthook(*sys.exc_info())
                fail_list.append(beam)
            ROOT.gROOT.SetBatch(True)
        if batch_level > 1:
            ROOT.gROOT.SetBatch(False)
        try:
            if len(dir_list) > 1:
                merge = ConglomerateMerge(conglomerate_list)
                merge.merge_all(rows, cols)
        except Exception:
            sys.excepthook(*sys.exc_info())
        fail_dict[ConfigClass.__name__] = fail_list
    return fail_dict

def print_fail_dict(fail_dict):
    print "Failed:"
    for key in fail_dict:
        print "    ", key, "had", len(fail_dict[key]), "fails"
        for fail in fail_dict[key]:
            print "        ", fail

def main_paper(batch_level = 0):
    """
    Main program; 
    - batch_level tells how much output for ROOT: batch_level 0 is silent, 10 is most verbose
    """
    fd_1, fd_2 = {}, {}
    root_style.setup_gstyle()
    ROOT.gROOT.SetBatch(True)
    target_dir = "output/2017-02-7-v14/"
    batch_level = 0
    hide_root_errors = True
    do_cuts_summary = True
    if batch_level < 10 and hide_root_errors:
        ROOT.gErrorIgnoreLevel = 6000
    my_dir_list = [
        ["2017-2.7_4-140_lH2_empty", "2017-2.7_6-140_lH2_empty", "2017-2.7_10-140_lH2_empty",],
        ["2017-2.7_4-140_lH2_full",  "2017-2.7_6-140_lH2_full",  "2017-2.7_10-140_lH2_full",],
        ["2017-2.7_4-140_None",      "2017-2.7_6-140_None",      "2017-2.7_10-140_None",],
        ["2017-2.7_4-140_LiH",       "2017-2.7_6-140_LiH",       "2017-2.7_10-140_LiH"],
    ]
    top_labels = ["4-140", "6-140", "10-140"]
    right_labels = ["Empty\nLH_{2}", "Full\nLH_{2}", "No\nabsorber", "LiH"]
    config_list = [CompareCutsConfig, CompareData1DConfig,
                   CompareOpticsConfig, CompareOpticsMCConfig,
                   CompareGlobalsConfig, CompareMCConfig,
                   CompareData2DConfig, CompareData2DMCConfig,
                  ] 
    config_list = [
                   CompareAmplitudeConfigBoth,
                   CompareAmplitudeConfigData,
                   CompareDensityConfig,
    ]
                   #CompareAmplitudeConfigMC,
    fd_1 = run_conglomerate(batch_level, config_list, my_dir_list, do_cuts_summary, target_dir, top_labels, right_labels)
    config_list = [CompareDensityConfig]#, CompareDensityRatioConfig,  CompareFractionalEmittanceConfig, 
    #fd_2 = run_conglomerate(batch_level, config_list, my_dir_list, False, target_dir, top_labels, right_labels)
    print_fail_dict(fd_1)
    print_fail_dict(fd_2)


if __name__ == "__main__":
    main_paper()
    if not ROOT.gROOT.IsBatch():
        raw_input("Finished - press <CR> to end")
