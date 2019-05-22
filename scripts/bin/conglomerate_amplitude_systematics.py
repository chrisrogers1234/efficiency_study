import copy
import sys
import os
import shutil
import json

import ROOT

import utilities.utilities as utilities
import utilities.root_style as root_style
import conglomerate
from conglomerate.compare_config import CompareConfig
from conglomerate.conglomerate_merge import ConglomerateMerge
from conglomerate.conglomerate_one import ConglomerateOne
from conglomerate.conglomerate_one import ConglomerateContainer
from conglomerate.merge_cuts_summary_tex import MergeCutsSummaryTex

class CompareCutsSystematicConfig(CompareConfig):
    def __init__(self, beam, absorber, sys_abs, systematic, ref):
        self.sys_abs = sys_abs
        self.ref = ref
        self.systematic = systematic
        if type(self.systematic) == type(""):
            self.systematic = [self.systematic]
        self.beam = beam
        self.absorber = absorber

    def set_dirs(self, data_src, mc_src, sys_src, output_target, dir_name, anal_dir):
        self.data_src = data_src
        self.mc_src = mc_src
        self.sys_src = sys_src
        self.analysis_dir = anal_dir
        self.output_target = output_target
        self.dir_name = dir_name
        self.output_dir = self.output_target+"/"+self.dir_name+"/"+self.beam+"_"+self.absorber+"/"

    def set_hist_data(self, top_labels, right_labels, plot, hist, x_range, normalise):
        self.top_labels = top_labels
        self.right_labels = right_labels
        self.plot = plot
        self.hist = hist
        self.x_range = x_range
        self.normalise = normalise

    def finalise_setup(self):
        sys_root = self.sys_src+"/plots_Simulated_"+self.beam+"_"+self.sys_abs+"_Systematics_"
        dir_list = [
            self.data_src+"plots_"+self.beam+"_"+self.absorber+"/",
            self.mc_src+"plots_Simulated_"+self.beam+"_"+self.absorber+"/",
            sys_root+self.ref+"/",
        ]
        for sys in self.systematic:
            dir_list.append(sys_root+sys+"/")

        n_sys = 1+len(self.systematic)
        sys_colors = [46, 38, 8, 9][:n_sys+1]
        sys_styles = [21, 34, 28, 33][:n_sys+1]
        for i in range(n_sys-1):
            if self.systematic[i] == self.ref:
                sys_colors[i+1] = sys_colors[0]
                sys_styles[i+1] = 1
        title = (self.beam+" "+self.absorber+"    "+self.systematic[0]).replace("_", " ")
        modifiers = {
            "merge_options":{
                "top_labels":self.top_labels,
                "right_labels":self.right_labels,
            },
            "hist_title":title,
            "file_name":self.plot,
            "canvas_name":self.plot,
            "normalise_hist":self.normalise,
            "histogram_names":[self.hist],
            "mice_logo":False,
            "legend":False,
            "calculate_errors":[],
            "redraw":{
                "draw_option":["P", "H"]+["P"]*n_sys,
                "fill_color":[1, ROOT.kOrange-2]+[1]*n_sys,
                "transparency":None,
                "line_color":[1, 1]+[1]*n_sys,
                "marker_style":[20, 20]+sys_styles,
                "marker_color":[1, 1]+sys_colors, # base is light red; sys is light blue
                "draw_order":[1]+range(2, 2+n_sys)+[0],
                "x_range":self.x_range,
                "y_range":None,
                "graph_draw_option":None,
                "ignore_more_histograms":False,
            },
            "rescale_x":self.x_range,
            "rescale_y":True,
            "write_plots":{
                "dir":self.output_dir,
                "file_name":self.plot+"_vs_"+self.systematic[0]
            },
            "axis_title":{"x":self.labels[self.plot]},
        }
        self.conglomerate_list = [
            self.get_conglomerate_0(modifiers = modifiers),
        ]
        self.data_caption = [[],]
        self.setup(self.beam, self.output_target, self.analysis_dir, self.dir_name, dir_list)

    labels = {
        "tku_p":"p at TKU Reference Plane [MeV/c]",
        "tku_x":"x at TKU Reference Plane [mm]",
        "tku_y":"y at TKU Reference Plane [mm]",
        "tku_px":"px at TKU Reference Plane [MeV/c]",
        "tku_py":"py at TKU Reference Plane [MeV/c]",
        "chi2_tku":"#chi^{2}/D.o.F. in TKU",
        "chi2_tkd":"#chi^{2}/D.o.F. in TKD",
        "p_res":"p(TKU) - p(TKD) [MeV/c]",
        "tkd_p":"p at TKD Reference Plane [MeV/c]",
        "tkd_x":"x at TKD Reference Plane [mm]",
        "tkd_y":"y at TKD Reference Plane [mm]",
        "tkd_px":"px at TKD Reference Plane [MeV/c]",
        "tkd_py":"py at TKD Reference Plane [MeV/c]",
        "tof01":"ToF_{01} [ns]",
        "tkd_max_r*":"Max R in TKD [mm]",
    }

def systematics_cut_summary(target_dir):
    dir_list = [
        "2017-2.7_4-140_lH2_empty_Systematics_tku_base",
        "2017-2.7_6-140_lH2_empty_Systematics_tku_base",
        "2017-2.7_10-140_lH2_empty_Systematics_tku_base", 
    ]
    mc_cuts_summary = MergeCutsSummaryTex()
    mc_cuts_summary.table_ref_pre = "systematics_"
    for beam in dir_list:
        config = CompareConfig()
        dir_list = [
            #target_dir+"plots_"+beam+"/",
            target_dir+"plots_Simulated_"+beam
        ]
        config.setup(beam, target_dir, "cut_plots/", "compare_cuts/", dir_list)
        config.mc_caption = [["" for i in range(5)] for i in range(5)]
        mc_prefix = ["upstream reconstructed", "downstream reconstructed", "extrapolated reconstructed",
                     "upstream truth", "downstream truth"]
        for i, prefix in enumerate(mc_prefix):
            for j, item in enumerate(config.mc_caption[i]):
                config.mc_caption[i][j] = "The "+prefix+" simulated sample is listed. "+item
        mc_cuts_summary.append_summary(config, [0])

    mc_cuts_summary.caption = config.mc_caption
    mc_cuts_summary.merge_summaries(target_dir+"/cuts_summary/mc/", "mc_cuts_summary")
    #raw_input("Done systematics - press <CR> to continue to the rest")

class SystematicsConglomerate(object):
    def __init__(self, top_labels, right_labels, target_dir, sys_run, dir_name):
        self.top_labels = top_labels
        self.right_labels = right_labels
        self.target_dir = target_dir
        self.sys_run = sys_run
        self.dir_name = dir_name

  
    def lh2_empty(self, dir_lists):
        sys_abs = "lH2_empty"
        ref = "tku_base"
        sys_list =  []
        for systematic in ["tku_base", "tku_scale_SSUE2_plus",
                           ["tku_scale_SSUC_plus", "tku_scale_SSUC_neg"],
                           "tku_scale_SSUE1_plus", "tku_pos_plus", "tku_rot_plus",
                           "tku_density_plus", "tku_full-p",
                          ]:
            sys_list += [
                (ref, systematic, sys_abs, "tof01", "tof01 us cut",   [2., 6.],      True, "data_plots/"),
                (ref, systematic, sys_abs, "tku_p", "tku_p us cut",   [130., 150.],  False,"data_plots/"),
                (ref, systematic, sys_abs, "tku_x", "tku_x us cut",   [-150., 150.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tku_y", "tku_y us cut",   [-150., 150.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tku_px", "tku_px us cut", [-100., 100.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tku_py", "tku_py us cut", [-100., 100.], False,"data_plots/"),
                (ref, systematic, sys_abs, "chi2_tku", "chi2 us cut", [0., 5.],      True, "data_plots/"),
                (ref, systematic, sys_abs, "p_res", "p_res us cut",   [-25., 25.],   True, "data_plots/"),
                (ref, systematic, sys_abs, "tkd_max_r*", "tkd_max_r", [0., 200.], True,"cut_plots/"),
            ]
        for systematic in ["tku_base", "tkd_scale_SSDE2_plus", "tkd_scale_SSDC_plus",
                           ["tkd_scale_SSDC_plus", "tkd_scale_SSDC_neg"],
                           "tkd_scale_SSDE1_plus", "tkd_pos_plus", "tkd_rot_plus",
                           "tkd_density_plus"]:
            sys_list += [
                (ref, systematic, sys_abs, "chi2_tkd", "chi2 ds cut", [0., 10.],     True, "data_plots/"),
                (ref, systematic, sys_abs, "p_res", "p_res ds cut",   [-50., 50.],   True, "data_plots/"),
                (ref, systematic, sys_abs, "tkd_p", "tkd_p ds cut",   [120., 160.],  False,"data_plots/"),
                (ref, systematic, sys_abs, "tkd_x", "tkd_x ds cut",   [-150., 150.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tkd_y", "tkd_y ds cut",   [-100., 100.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tkd_px", "tkd_px ds cut", [-100., 100.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tkd_py", "tkd_py ds cut", [-100., 100.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tkd_max_r*", "tkd_max_r", [0., 200.], True,"cut_plots/"),
            ]
        self.run_systematics(dir_lists, sys_list)


    def lh2_full(self, dir_lists):
        ref = "mc_base"
        sys_abs = "lH2_full"
        sys_list = []
        for systematic in [["mc_beam_offset_plus", "mc_beam_offset_minus"],
                          "mc_base", 
                          "mc_fc_plus", "mc_ssu_match_plus", "mc_ssd_match_plus",
                          #"mc_lh2_plus"
                          ]:
            sys_list += [
                (ref, systematic, sys_abs, "p_res", "p_res us cut",   [-10., 40.],   True, "data_plots/"),
                (ref, systematic, sys_abs, "tku_p", "tku_p",   [130., 150.],  False,"data_plots/"),
                (ref, systematic, sys_abs, "tku_x", "tku_x",   [-150., 150.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tku_y", "tku_y",   [-150., 150.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tku_px", "tku_px", [-100., 100.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tku_py", "tku_py", [-100., 100.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tkd_p", "tkd_p",   [105., 145.],  False,"data_plots/"),
                (ref, systematic, sys_abs, "tkd_x", "tkd_x",   [-150., 150.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tkd_y", "tkd_y",   [-100., 100.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tkd_px", "tkd_px", [-100., 100.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tkd_py", "tkd_py", [-100., 100.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tkd_max_r*", "tkd_max_r", [0., 200.], True,"cut_plots/"),
            ]
        self.run_systematics(dir_lists, sys_list)

    def lih(self, dir_lists):
        ref = "mc_base"
        sys_abs = "LiH"
        sys_list = []
        for systematic in ["mc_base", "mc_lih_plus"]:
            sys_list += [
                (ref, systematic, "LiH",   "p_res",  "p_res",  [-50., 50.],   True, "data_plots/"),
                (ref, systematic, sys_abs, "tkd_p",  "tkd_p",  [100., 160.],  False,"data_plots/"),
                (ref, systematic, sys_abs, "tkd_x",  "tkd_x",  [-100., 100.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tkd_y",  "tkd_y",  [-100., 100.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tkd_px", "tkd_px", [-100., 100.], False,"data_plots/"),
                (ref, systematic, sys_abs, "tkd_py", "tkd_py", [-100., 100.], False,"data_plots/"),
            ]
        self.run_systematics(dir_lists, sys_list)

    def run_systematics(self, dir_lists, sys_list):
        rows = len(dir_lists)
        cols = min([len(sub_list) for sub_list in dir_lists])
        dir_list = []
        for sub_list in dir_lists:
            dir_list += sub_list

        first = True
        for ref, systematic, sys_abs, fname, hist, x_range, normalise, anal_dir in sys_list:
            conglomerate_list = []
            for beam in dir_list:
                [beam, absorber] = beam.split("140_")
                beam += "140"
                try:
                    config = CompareCutsSystematicConfig(beam, absorber, sys_abs, systematic, ref)
                    config.set_dirs(self.target_dir, self.target_dir, self.sys_run, 
                                    self.target_dir, self.dir_name, anal_dir)
                    if first:
                        self.mkdirs(config.output_dir)
                    config.set_hist_data(self.top_labels, self.right_labels,
                                         fname, hist, x_range, normalise)
                    config.finalise_setup()
                    cong = ConglomerateContainer(config)
                    cong.conglomerate()
                    conglomerate_list.append(cong)
                except Exception:
                    sys.excepthook(*sys.exc_info())
            first = False
            try:
                if len(dir_list) > 1:
                    print "Merging"
                    merge = ConglomerateMerge(conglomerate_list)
                    merge.merge_all(rows, cols)
            except Exception:
                sys.excepthook(*sys.exc_info())

    def mkdirs(self, output_dir):
        print output_dir
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir)



def main():
    target_dir = "output/2017-02-7-v11/"
    systematics_source_dir = "output/2017-02-7-Systematics-v5/"
    systematics_cut_summary(systematics_source_dir)
    root_style.setup_gstyle()
    ROOT.gROOT.SetBatch(True)
    top_labels = ["4-140", "6-140", "10-140"]

    dir_name = "compare_recon_systematics"
    lh2_empty_dir_list = [["2017-2.7_4-140_lH2_empty", "2017-2.7_6-140_lH2_empty", "2017-2.7_10-140_lH2_empty",],]
    right_labels = ["Empty\nLH2",]
    sys_conglomerate = SystematicsConglomerate(top_labels, right_labels, target_dir, systematics_source_dir, dir_name)
    sys_conglomerate.lh2_empty(lh2_empty_dir_list)

    dir_name = "compare_performance_systematics"
    lh2_full_dir_list = [["2017-2.7_4-140_lH2_full", "2017-2.7_6-140_lH2_full", "2017-2.7_10-140_lH2_full",],]
    right_labels = ["Full\nLH2",]
    sys_conglomerate = SystematicsConglomerate(top_labels, right_labels, target_dir, systematics_source_dir, dir_name)
    sys_conglomerate.lh2_full(lh2_full_dir_list)

    lih_full_dir_list = [["2017-2.7_4-140_LiH", "2017-2.7_6-140_LiH", "2017-2.7_10-140_LiH",],]
    right_labels = ["LiH",] 
    sys_conglomerate = SystematicsConglomerate(top_labels, right_labels, target_dir, systematics_source_dir, dir_name)
    #sys_conglomerate.lih(lih_full_dir_list)



if __name__ == "__main__":
    main()
    if not ROOT.gROOT.IsBatch():
        raw_input("Finished - press <CR> to end")
