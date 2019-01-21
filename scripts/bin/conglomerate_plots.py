import copy
import sys
import shutil

import ROOT

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
        "extra_labels":{
            "right":right_labels,
            "top":top_labels
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
            target_dir+"plots_Simulated_"+beam
        ]
        self.setup(beam, target_dir, "cut_plots/", "compare_cuts/", dir_list)
        mod = {
            "extra_labels":{
                "right":right_labels,
                "top":top_labels
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
            self.get_conglomerate_2("tkd_chi2_"+nd+"_0", None, "#chi^{2}/D.o.F. in TKD", None, True, [0.5, 0.5, 0.9, 0.9], vertical([4], mod)),
            self.get_conglomerate_2("tkd_max_r_"+nd+"_0", None, "Maximum radius in TKD stations [mm]", [0., 300.], True, [0.5, 0.5, 0.9, 0.9], vertical([150], mod)),
            self.get_conglomerate_2("tkd_p_"+nd1+"_0", None, "Momentum at TKD Reference Plane [MeV/c]", [50., 250.], True, [0.5, 0.5, 0.9, 0.9], vertical([90, 170], mod)),
            self.get_conglomerate_2("tku_chi2_"+nu+"_0", None, "#chi^{2}/D.o.F. in TKU", None, True, [0.5, 0.5, 0.9, 0.9], vertical([4], mod)),
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
            " Samples are listed for 3-140 and 4-140 datasets.",
            " Samples are listed for 6-140 and 10-140 datasets.",
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
            "extra_labels":{
                "right":right_labels,
                "top":top_labels
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
            "extra_labels":{
                "right":right_labels,
                "top":top_labels
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

class CompareData2DMCConfig(CompareConfig):
    def __init__(self, beam, target_dir, top_labels, right_labels):
        dir_list = [
            #target_dir+"plots_"+beam+"/",
            target_dir+"plots_Simulated_"+beam
        ]
        self.setup(beam, target_dir, "data_plots/", "compare_data_2d_mc/", dir_list)
        mod = {
            "extra_labels":{
                "right":right_labels,
                "top":top_labels
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
            "extra_labels":{
                "right":right_labels,
                "top":top_labels
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
            "extra_labels":{
                "right":right_labels,
                "top":top_labels
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
            "extra_labels":{
                "right":right_labels,
                "top":top_labels
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

class CompareAmplitudeConfigMC(CompareConfig): # MC corrections
    def __init__(self, beam, target_dir, top_labels, right_labels):
        dir_list = [
            #target_dir+"plots_"+beam+"/",
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
            "extra_labels":{
                "right":right_labels,
                "top":top_labels
            }
        }
        absolute_modifiers = {
            "extra_labels":{
                "right":right_labels,
                "top":top_labels
            },
            "normalise_graph":True,
            "redraw":{
                    "graph":{
                        "draw_option":["p", "2", "p", "2"],
                        "marker_style":[20, 20, 22, 22],
                        "marker_color":[ROOT.kOrange+4, ROOT.kOrange+4, ROOT.kGreen+3, ROOT.kGreen+3],
                        "fill_color":[ROOT.kOrange+4, ROOT.kOrange+4, ROOT.kGreen+3, ROOT.kGreen+3],
                        "draw_order":[0, 1, 2, 3],
                    }
            },

        }

        self.conglomerate_list = [
            #self.get_conglomerate_3("crossing_probability_downstream", "crossing_probability_downstream", None, None),
            #self.get_conglomerate_3("crossing_probability_upstream", "crossing_probability_upstream", None, None),
            #self.get_conglomerate_graph("amplitude_inefficiency_all_upstream", "US True Amplitude [mm]",
            #                            "#frac{Number in MC True (all)}{Number in MC True (reco)}",
            #                            "inefficiency_all_upstream", ["inefficiency"], ["inefficiency"]), 
            #self.get_conglomerate_graph("amplitude_inefficiency_all_downstream", "DS True Amplitude [mm]",
            #                            "#frac{MC True (all) Number}{MC True (reco) Number}", 
            #                            "inefficiency_all_downstream", ["inefficiency"], ["inefficiency"]), 
            self.get_conglomerate_graph("amplitude_pdf_reco", "Reconstructed Amplitude [mm]",
                                        "Number",
                                        "amplitude_pdf_reco", ["Upstream sys hist"],
                                        ["Upstream stats", "Upstream sys", "Downstream stats", "Downstream sys"],
                                        y_range = [0.001, 1.399], x_range = [0.01, 99.9], replace_hist = True,
                                        modifiers = absolute_modifiers),
            self.get_conglomerate_graph("amplitude_cdf_reco", "Reconstructed Amplitude [mm]",
                                        "Cumulative density",
                                        "amplitude_cdf_reco", ["Upstream CDF stats hist"],
                                        ["Upstream CDF stats", "Upstream CDF sys", "Downstream CDF stats", "Downstream CDF sys"],
                                        y_range = [0.001, 1.399], x_range = [0.01, 99.9], replace_hist = True,
                                        modifiers = absolute_modifiers),

            self.get_conglomerate_graph("pdf_ratio*", "Reconstructed amplitude [mm]",
                                        "P_{Amp}",
                                        "pdf_ratio", ["PDF Ratio stats hist"],
                                        ["PDF Ratio stats", "PDF Ratio sys"], x_range = [0.01, 99.9], y_range = [0.501, 2.499], replace_hist = True,
                                        graph_draw_option = ["p", "2"], graph_marker_style=[20, 20], graph_marker_color=[ROOT.kRed, ROOT.kRed], graph_draw_order=[1,0], modifiers=ratio_modifiers),
            self.get_conglomerate_graph("cdf_ratio*", "Reconstructed amplitude [mm]",
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
                "verticals":[],
            },
            "extra_labels":{
                "right":right_labels,
                "top":top_labels
            },
        }
        absolute_modifiers = {
            "extra_labels":{
                "right":right_labels,
                "top":top_labels
            },
            "normalise_graph":True,
            "redraw":{
                    "graph":{
                        "draw_option":["p", "2", "p", "2"],
                        "marker_style":[20, 20, 22, 22],
                        "marker_color":[ROOT.kOrange+4, ROOT.kOrange+4, ROOT.kGreen+3, ROOT.kGreen+3],
                        "fill_color":[ROOT.kOrange+4, ROOT.kOrange+4, ROOT.kGreen+3, ROOT.kGreen+3],
                        "draw_order":[0, 1, 2, 3],
                    }
            },

        }

        self.conglomerate_list = [
            self.get_conglomerate_graph("amplitude_pdf_reco", "Reconstructed Amplitude [mm]",
                                        "Number",
                                        "amplitude_pdf_reco", ["Upstream sys hist"],
                                        ["Upstream stats", "Upstream sys", "Downstream stats", "Downstream sys"],
                                        y_range = [0.001, 1.399], x_range = [0.01, 99.9], replace_hist = True,
                                        modifiers = absolute_modifiers),
            self.get_conglomerate_graph("amplitude_cdf_reco", "Reconstructed Amplitude [mm]",
                                        "Cumulative density",
                                        "amplitude_cdf_reco", ["Upstream CDF stats hist"],
                                        ["Upstream CDF stats", "Upstream CDF sys", "Downstream CDF stats", "Downstream CDF sys"],
                                        y_range = [0.001, 1.399], x_range = [0.01, 99.9], replace_hist = True,
                                        modifiers = absolute_modifiers),
            self.get_conglomerate_graph("pdf_ratio*", "Reconstructed amplitude [mm]",
                                        "P_{Amp}",
                                        "pdf_ratio", ["PDF Ratio stats hist"],
                                        ["PDF Ratio stats", "PDF Ratio sys"], x_range = [0.01, 99.9], y_range = [0.501, 2.499], replace_hist = True,
                                        graph_draw_option = ["p", "2"], graph_marker_style=[20, 20], graph_marker_color=[1, 1], graph_draw_order=[1,0], modifiers=ratio_modifiers),
            self.get_conglomerate_graph("cdf_ratio*", "Reconstructed amplitude [mm]",
                                        "R_{Amp}",
                                        "cdf_ratio", ["CDF_Ratio_stats hist"],
                                        ["CDF_Ratio stats", "CDF_Ratio sys"], x_range = [0.01, 99.9], y_range = [0.601, 1.399], replace_hist = True,
                                        graph_draw_option = ["p", "2"], graph_marker_style=[20, 20], graph_marker_color=[1, 1], graph_draw_order=[1,0], modifiers=ratio_modifiers),
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
                "verticals":[],
            },
            "extra_labels":{
                "right":right_labels,
                "top":top_labels
            },
            "redraw":{
                "transparency":[0.5, 0.5, 0.5, 0.5],
                "graph":{
                    #"fill_style":[4050, 4050, 4050, 4050],
                    "fill_color":[ROOT.kBlack, ROOT.kBlack, ROOT.kRed, ROOT.kRed],
                    "transparency":[0.5, 0.5, 0.5, 0.5],
                }
            }
        }

        self.conglomerate_list = [
            self.get_conglomerate_graph("pdf_ratio*", "Reconstructed Amplitude [mm]",
                                        "#frac{Number out}{Number in}",
                                        "pdf_ratio", ["PDF Ratio stats_hist"],
                                        ["PDF_Ratio stats", "PDF_Ratio sys"], x_range = [0.01, 99.9], y_range = [0.01, 2.001], replace_hist = True,
                                        graph_draw_option = ["p 3", "3", "p 3", "3"], graph_marker_style=[20, 20, 22, 22], 
                                        graph_marker_color=[1, 1, ROOT.kRed, ROOT.kRed], graph_draw_order=[3, 2, 1, 0,], 
                                        modifiers=ratio_modifiers),
            self.get_conglomerate_graph("cdf_ratio*", "Reconstructed amplitude [mm]",
                                        "R_{Amp}",
                                        "cdf_ratio", ["CDF_Ratio_stats hist"],
                                        ["CDF_Ratio stats", "CDF_Ratio sys"], x_range = [0.01, 99.9], y_range = [0.601, 1.399], replace_hist = True,
                                        graph_draw_option = ["p", "3", "p", "3"], graph_marker_style=[20, 20, 22, 22], 
                                        graph_marker_color=[1, 1, ROOT.kRed, ROOT.kRed], graph_draw_order=[3, 2, 1, 0,], 
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
            if len(dir_lists) > 1:
                merge = ConglomerateMerge(conglomerate_list)
                merge.merge_all(rows, cols)
        except Exception:
            sys.excepthook(*sys.exc_info())
        fail_dict[ConfigClass.__name__] = fail_list
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
    root_style.setup_gstyle()
    ROOT.gROOT.SetBatch(True)
    target_dir = "output/2017-02-7-test/"
    batch_level = 0
    hide_root_errors = True
    do_cuts_summary = True
    if batch_level < 10 and hide_root_errors:
        ROOT.gErrorIgnoreLevel = 6000
    my_dir_list = [
        ["2017-2.7_4-140_None",      "2017-2.7_6-140_None",      "2017-2.7_10-140_None",],
        ["2017-2.7_4-140_lH2_empty", "2017-2.7_6-140_lH2_empty", "2017-2.7_10-140_lH2_empty",],
        ["2017-2.7_4-140_lH2_full",  "2017-2.7_6-140_lH2_full",  "2017-2.7_10-140_lH2_full",],
        ["2017-2.7_4-140_LiH",       "2017-2.7_6-140_LiH",       "2017-2.7_10-140_LiH"],
    ]
    top_labels = ["4-140", "6-140", "10-140"]
    right_labels = ["Empty\nLH2", "Full\nLH2", "LiH"] #"No\nabsorber", 
    config_list = [CompareData1DConfig, CompareCutsConfig, 
                   CompareOpticsConfig, CompareOpticsMCConfig,
                   CompareGlobalsConfig, CompareMCConfig,
                   CompareData2DConfig, CompareData2DMCConfig,]
    config_list += [CompareAmplitudeConfigBoth,
                   CompareAmplitudeConfigMC,
                   CompareAmplitudeConfigData]
    run_conglomerate(batch_level, config_list, my_dir_list, do_cuts_summary, target_dir, top_labels, right_labels)

if __name__ == "__main__":
    main_paper()
    if not ROOT.gROOT.IsBatch():
        raw_input("Finished - press <CR> to end")
