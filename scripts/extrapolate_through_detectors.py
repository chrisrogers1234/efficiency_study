import copy
import sys
import os
import site
import json
import datetime

import numpy
import ROOT
analysis_path = os.environ["MAUS_ROOT_DIR"]+"/bin/user/fit_fields/scripts/"
site.addsitedir(analysis_path)
import xboa.common
from xboa.hit import Hit
from xboa.bunch import Bunch
import Configuration
import maus_cpp.global_error_tracking
import maus_cpp.globals
import maus_cpp.field
import libxml2
import xboa.common

class ExtrapolateTrackPoints(object):
    def __init__(self, config, config_anal, data_loader):
        self.config = config
        self.config_anal = config_anal
        self.current_space_points = [] # list of space points in current event
        self.fitted_track_points = [] # list of fitted space points in current event
        self.best_guess = None
        self.iteration = 0
        self.mu_mass = xboa.common.pdg_pid_to_mass[13]
        self.setup_maus_config()
        self.residual_dicts = {}
        self.miss_lists = {}
        self.weighted_residual_dicts = {}
        self.data_loader = data_loader
        self.vb_out = None
       
    def setup_maus_config(self):
        self.tracking = maus_cpp.global_error_tracking.GlobalErrorTracking()
        self.tracking.set_min_step_size(self.config.global_min_step_size)
        self.tracking.set_max_step_size(self.config.global_max_step_size)
        self.tracking.set_energy_loss_model("bethe_bloch_forwards")
        self.tracking.set_scattering_model("moliere_backwards")
        self.tracking.set_geometry_model("axial_lookup")
        for material in ["AIR", "POLYSTYRENE", "POLYCARBONATE", "POLYVINYL_TOLUENE",
                        "POLYVINYL_ACETATE", "AEROGEL_112a", "AEROGEL_107a",
                        "Al", "CELLULOSE_CELLOPHANE", "POLYVINYL_CHLORIDE_LOWD",
                        "POLYVINYL_CHLORIDE", "He"]:
            maus_cpp.global_error_tracking.enable_material(material)

    def get_point(self, event, detector):
        if len(event) == 0:
            raise ValueError("Event was empty")
        for rec_point in event:
            if rec_point["detector"] == detector:
                break
        if rec_point["detector"] != detector:
            detectors = str([rec["detector"] for rec in event])
            raise ValueError("Failed to find rec point for "+str(detector)+" from "+ detectors)
        return rec_point


    def make_rec(self, value, ellipse, detector):
        hit = Hit.new_from_dict( {
            "t":value[0],
            "x":value[1],
            "y":value[2],
            "z":value[3],
            "px":value[5],
            "py":value[6],
            "pz":value[7],
            "pid":-13,
            "mass":xboa.common.pdg_pid_to_mass[13],
          }, "energy"
        )
        rec_point = {
            "detector":self.global_key(detector),
            "covariance":ellipse,
            "hit":hit
        }
        return rec_point

    def global_key(self, detector):
        return "global_"+detector

    def extrapolate_event(self, event):
        event = sorted(event, key = lambda tp: tp["hit"]["z"]) # check that it is sorted
        tku_rec = self.get_point(event, "tku_tp")
        energy = (tku_rec["hit"]["px"]**2 + tku_rec["hit"]["py"]**2 + tku_rec["hit"]["pz"]**2 + self.mu_mass**2)**0.5
        seed = [0.,     tku_rec["hit"]["x"],  tku_rec["hit"]["y"],  tku_rec["hit"]["z"],
                tku_rec["hit"]["energy"], tku_rec["hit"]["px"], tku_rec["hit"]["py"], tku_rec["hit"]["pz"],]
        ellipse = tku_rec["covariance"]
        # walking upstream
        try:
            tof1_point, tof1_error = self.tracking.propagate_errors(seed, ellipse, self.config.z_tof1)
            tof1_point[0] = 0.
            event.append(self.make_rec(tof1_point, tof1_error, "tof1"))
            tof0_point, tof0_error = self.tracking.propagate_errors(tof1_point, tof1_error, self.config.z_tof0)
            event.append(self.make_rec(tof0_point, tof0_error, "tof0"))
        except RuntimeError:
            pass #sys.excepthook(*sys.exc_info())

        # walking downstream
        try:
            tkd_point, tkd_error = self.tracking.propagate_errors(seed, ellipse, self.config.z_tkd_1)
            event.append(self.make_rec(tkd_point, tkd_error, "tkd_tp"))
            tof2_point, tof2_error = self.tracking.propagate_errors(tkd_point, tkd_error, self.config.z_tof2)
            event.append(self.make_rec(tof2_point, tof2_error, "tof2"))
        except RuntimeError:
            pass #sys.excepthook(*sys.exc_info())
            
        return event

    def alt_min_max(self, float_list, n_sigma):
        length = len(float_list)+1
        while len(float_list) < length:
            mean = numpy.mean(float_list)
            sigma = numpy.std(float_list)
            min_max = [
                mean-n_sigma*sigma,
                mean+n_sigma*sigma,
            ]
            length = len(float_list)
            float_list = [x for x in float_list if x < min_max[1] and x > min_max[0]]
        return min_max

    def append_misses(self, event, detector):
        if detector not in self.miss_lists:
            self.miss_lists[detector] = []
        event = event["data"]
        try:
            det_rec = self.get_point(event, detector)
            return
        except ValueError:
            # we successfully extrapolated the track, but it is not observed in
            # the detector... I wonder why?
            pass

        try:
            glob_rec = self.get_point(event, self.global_key(detector))["hit"]
        except ValueError:
            # we failed to extrapolate a track, nothing more to do
            return
        self.miss_lists[detector].append(glob_rec)


    def append_residual(self, event, detector, axis):
        event = event["data"]
        if detector not in self.residual_dicts:
            self.residual_dicts[detector] = {}
            self.weighted_residual_dicts[detector] = {}
        if axis not in self.residual_dicts[detector]:
            self.residual_dicts[detector][axis] = []
            self.weighted_residual_dicts[detector][axis] = []
        
        residual_list = self.residual_dicts[detector][axis]
        weighted_residual_list = self.weighted_residual_dicts[detector][axis]

        glob_rec, det_rec = None, None
        try:
            glob_rec = self.get_point(event, self.global_key(detector))
        except ValueError:
            # we failed to extrapolate a track, nothing more to do
            return

        try:
            det_rec = self.get_point(event, detector)
        except ValueError:
            # we successfully extrapolated the track, but it is not observed in
            # the detector... I wonder why?
            #sys.excepthook(*sys.exc_info())
            return

        # we have an extrapolated track and a detector hit. How close is the
        # extrapolated track to the measured particle?
        residual = det_rec["hit"][axis]-glob_rec["hit"][axis]
        residual_list.append(residual)
        if axis in self.cov_key_list:
            key = self.cov_key_list.index(axis)
            sigma = (glob_rec["covariance"][key][key]+det_rec["covariance"][key][key])**0.5 # -ve sqrt exception?
            if sigma < 1e-12:
                raise ValueError("Zero RMS")
            weighted_residual = residual/sigma
            weighted_residual_list.append(weighted_residual)

    def print_canvas(self, canvas, name):
        canvas.Update()
        plot_name = os.path.join(self.config_anal["plot_dir"], name.replace(' ', '_'))
        for format in "png", "pdf", "root":
            canvas.Print(plot_name+"."+format)

    def get_geometry(self):
        path_to_info_file = self.config.info_file
        info = libxml2.parseFile(path_to_info_file)
        path = "gdml/MICE_Information/Configuration_Information/GeometryID"
        geo_id = info.xpathEval(path)[0].prop("value")
        info.freeDoc()
        return geo_id

    text_boxes = []
    def get_text_box(self, residual_list):
        text_box = ROOT.TPaveText(0.6, 0.4, 0.9, 0.9, "NDC")
        text_box.SetFillColor(0)
        text_box.SetBorderSize(0)
        text_box.SetTextSize(0.04)
        text_box.SetTextAlign(12)
        text_box.SetTextSize(0.03)
        text_box.AddText(self.config_anal["name"])
        text_box.AddText("Recon: "+self.data_loader.maus_version)
        text_box.AddText("Geometry: "+str(self.get_geometry()))
        text_box.AddText("All events (black)")
        text_box.AddText("Number: "+str(len(residual_list)))
        text_box.AddText("Mean:   "+str(numpy.mean(residual_list)))
        text_box.AddText("Std:    "+str(numpy.std(residual_list)))
        text_box.SetBorderSize(1)
        text_box.Draw()
        self.text_boxes.append(text_box)
        return text_box

    def plot_residuals(self, detector, axis):
        if detector not in self.residual_dicts:
            print "No residuals for", detector
            return
        residual_list = self.residual_dicts[detector][axis]
        weighted_residual_list = self.weighted_residual_dicts[detector][axis]
        if False: #axis == "pz":
            print "Residuals", residual_list
            print "Weighted Residuals", weighted_residual_list
            print "Det", detector, "axis", axis
        if len(residual_list) == 0:
            raise RuntimeError("No residuals were found? Something screwy?")

        nbins = self.config.residuals_plots_nbins
        name = "residuals - "+detector+" "+axis
        label = "Res("+axis+") ["+self.units[axis]+"]"
        min_max = self.alt_min_max(residual_list, 10)
        canvas = xboa.common.make_root_canvas(name)
        hist = xboa.common.make_root_histogram(name, residual_list, label, nbins, xmin = min_max[0], xmax = min_max[1])
        hist.SetTitle(detector+": "+axis)
        hist.Draw()
        self.get_text_box(residual_list)
        canvas.Update()
        self.print_canvas(canvas, name)

        name = "normalised "+name
        canvas = xboa.common.make_root_canvas(name)
        min_max = self.alt_min_max(weighted_residual_list, 10)
        if len(weighted_residual_list) > 0:
            hist = xboa.common.make_root_histogram(name, weighted_residual_list, "Res("+axis+")/#sigma("+axis+")", nbins, xmin = min_max[0], xmax = min_max[1])
            hist.SetTitle(detector+": "+axis)
            hist.Draw()
        self.get_text_box(weighted_residual_list)
        self.print_canvas(canvas, name)

    def get_miss_text_box(self, bunch, detector):
        text_box = ROOT.TPaveText(0.6, 0.4, 0.9, 0.9, "NDC")
        text_box.SetFillColor(0)
        text_box.SetBorderSize(0)
        text_box.SetTextSize(0.04)
        text_box.SetTextAlign(12)
        text_box.SetTextSize(0.03)
        text_box.AddText("Preliminary")
        text_box.AddText("Recon: "+self.data_loader.maus_version)
        text_box.AddText("Geometry: "+str(self.get_geometry()))
        text_box.AddText("All events (black)")
        text_box.AddText("Number missed: "+str(len(bunch)))
        n_rec = len(self.residual_dicts[detector].values()[0])
        text_box.AddText("Number recorded: "+str(n_rec))
        text_box.Draw()
        self.text_boxes.append(text_box)
        return text_box

    def draw_box(self, x0, y0, x1, y1):
        x = [x0, x0, x1, x1, x0]
        y = [y0, y1, y1, y0, y0]
        hist, graph = xboa.common.make_root_graph("", x, y, "", "")
        graph.Draw()
        return graph

    def plot_misses(self):
        for detector in "tof0", "tof1", "tof2":
            if detector not in self.miss_lists or len(self.miss_lists[detector]) == 0:
                print "No misses for", detector
                continue
            print detector
            bunch = Bunch.new_from_hits(self.miss_lists[detector])
            canvas, hist = bunch.root_histogram("x", "mm", nx_bins=20, xmin=-500., xmax=+500.)
            hist.SetTitle("Misses - "+detector)
            self.get_miss_text_box(bunch, detector)
            self.print_canvas(canvas, "misses_"+detector+"_x")
            canvas, hist = bunch.root_histogram("y", "mm", nx_bins=20, xmin=-500., xmax=+500.)
            hist.SetTitle("Misses - "+detector)
            self.get_miss_text_box(bunch, detector)
            self.print_canvas(canvas, "misses_"+detector+"_y")
            canvas, hist = bunch.root_histogram("x", "mm", "y", "mm", nx_bins=20, ny_bins=20, xmin=-500., xmax=+500., ymin=-500., ymax=+500.)
            hist.SetTitle("Misses - "+detector)
            self.get_miss_text_box(bunch, detector)
            self.print_canvas(canvas, "misses_"+detector+"_x-y")
        for detector in "tkd_tp",:
            if detector not in self.miss_lists or len(self.miss_lists[detector]) == 0:
                print "No misses for", detector
                continue
            print detector
            bunch = Bunch.new_from_hits(self.miss_lists[detector])
            canvas, hist = bunch.root_histogram("r", "mm", nx_bins=10, xmin=0., xmax=+300.)
            hist.SetTitle("Misses - "+detector)
            self.get_miss_text_box(bunch, detector)
            self.print_canvas(canvas, "misses_"+detector+"_r")
            canvas, hist = bunch.root_histogram("pt", "MeV/c", nx_bins=10, xmin=0., xmax=+100.)
            hist.SetTitle("Misses - "+detector)
            self.get_miss_text_box(bunch, detector)
            self.print_canvas(canvas, "misses_"+detector+"_pt")
            canvas, hist = bunch.root_histogram("r", "mm", "pt", "MeV/c", nx_bins=10, ny_bins=10, xmin=0., xmax=+300., ymin=0., ymax=+100.)
            hist.SetTitle("Misses - "+detector)
            self.get_miss_text_box(bunch, detector)
            self.print_canvas(canvas, "misses_"+detector+"_r-pt")
            canvas, hist = bunch.root_histogram("x", "mm", "y", "mm", nx_bins=10, ny_bins=10, xmin=-300., xmax=+300., ymin=-300., ymax=+300.)
            hist.SetTitle("Misses - "+detector)
            self.get_miss_text_box(bunch, detector)
            self.print_canvas(canvas, "misses_"+detector+"_x-y")

    units = {"t":"ns", "x":"mm", "y":"mm", "z":"mm", "energy":"MeV", "px":"MeV/c", "py":"MeV/c", "pz":"MeV/c"}
    cov_key_list = ["t", "x", "y", "energy", "px", "py"]
    mu_mass = xboa.common.pdg_pid_to_mass[13]

    def print_hits_csv(self, detector_list, headers):
        if self.vb_out == None:
            self.vb_out = open("data_for_vb.csv", "w")
            for head in headers:
                if head in self.units:
                    units = " ["+self.units[head]+"]"
                else:
                    units = ""
                print >> self.vb_out, str(head+units).ljust(15),
            print >> self.vb_out
        print "    writing data...",
        sys.stdout.flush()
        for event in self.data_loader.events:
            for detector in detector_list:
                for hit in event["data"]:
                    if hit["detector"] != detector:
                        continue
                    for head in headers:
                        if head in hit.keys():
                            print >> self.vb_out, str(hit[head]).ljust(15),
                        elif head in event.keys():
                            print >> self.vb_out, str(event[head]).ljust(15),
                        else:
                            raise KeyError("Couldn't find key "+head+" from event keys: "+str(event.keys())+" hit keys: "+str(hit.keys()))
                    print >> self.vb_out
                    break
        self.vb_out.flush()
        print "written"

    def maus_globals(config):
        str_conf = Configuration.Configuration().\
                                          getConfigJSON(command_line_args=False)
        json_conf = json.loads(str_conf)
        json_conf["simulation_geometry_filename"] = config.geometry
        json_conf["max_step_length"] = 1.
        json_conf["physics_processes"] = "mean_energy_loss"
        json_conf["verbose_level"] = config.verbose_level
        maus_conf = json.dumps(json_conf)
        maus_cpp.globals.birth(maus_conf)

def is_cut(event, config):
    for cut in config.extrapolation_cuts:
        if config.extrapolation_cuts[cut] and event["will_cut"][cut]:
            return True
    return False

def do_plots(extrapolate):
        xboa.common.clear_root()
        extrapolate.plot_misses()
        for detector in "tof0", "tof1", "tof2":
            for axis in "x", "y", "t":
                extrapolate.plot_residuals(detector, axis)      
        for detector in "tkd_tp",:
            for axis in "x", "y", "px", "py", "pz":
                extrapolate.plot_residuals(detector, axis)      

def do_extrapolation(config, config_anal, data_loader):
    print "Doing extrapolation"
    index = 0
    pass_cuts = 0
    pass_extrapolate = 0
    step = 500
    will_do_plots = False
    extrapolate = ExtrapolateTrackPoints(config, config_anal, data_loader)
    while index < len(data_loader.events):
        try:
            while will_do_plots == False:
                index += 1
                event = data_loader.events[index]
                if is_cut(event, config):
                    continue
                pass_cuts += 1
                try:
                    event["data"] = extrapolate.extrapolate_event(event["data"])  
                    pass_extrapolate += 1
                except ValueError:
                    pass #sys.excepthook(*sys.exc_info())
                for detector in "tof1", "tof0", "tof2":
                    extrapolate.append_misses(event, detector)
                    for axis in "x", "y", "t":
                        extrapolate.append_residual(event, detector, axis)
                for detector in "tkd_tp",:
                    extrapolate.append_misses(event, detector)
                    for axis in "x", "y", "px", "py", "pz":
                        extrapolate.append_residual(event, detector, axis)
                will_do_plots = pass_cuts % step == 0
        except IndexError:
            index = len(data_loader.events) # force to finish
        print "event", index, "of", len(data_loader.events), "with", pass_cuts, "passing cuts and", pass_extrapolate, "extrapolating okay"
        do_plots(extrapolate)
        will_do_plots = False
        

