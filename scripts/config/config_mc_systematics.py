import copy
import xboa.common
from xboa.bunch import Bunch

def rogers_mc_file_names(datasets):
    file_list = ["/work/ast/cr67/"+datasets+"/maus_reconstruction.root"]
    return file_list

def reco_file_names(run_number_list, maus, do_globals):
    file_list = []
    for run in run_number_list:
        run = str(run).rjust(5, '0')
        a_file = "/home/cr67/work/reco/"+maus+"/"+run+"/"+run+"_recon"
        if not do_globals:
            a_file += ".root"
        else:
            a_file += "_global.root"
        file_list.append(a_file)
    print file_list
    return file_list

def get_analysis(datasets, name, tof01_min_max, data_dir, p_bins, tkd_cut, do_globals):
    plot_dir = data_dir+"/plots_"+name+"/"
    plot_dir = plot_dir.replace(" ", "_")
    min_p = min([min(a_bin) for a_bin in p_bins])
    max_p = max([max(a_bin) for a_bin in p_bins])
    cov_us = Bunch.build_penn_ellipse(
                        1., xboa.common.pdg_pid_to_mass[13],
                        250., 0., (min_p+max_p)/2., 0., 3.e-3, 1. 
            ).tolist()
    cov_ds = Bunch.build_penn_ellipse(
                        1., xboa.common.pdg_pid_to_mass[13],
                        250., 0., (min_p+max_p)/2., 0., 2.e-3, 1. 
            ).tolist()
    return {
            "plot_dir":plot_dir, # makedirs and then put plots in this directory. Removes any old plots there!!!
            "tof0_n_sp":1, # number of space points required in TOF0
            "tof1_n_sp":1, # number of space points required in TOF1
            "tof12_cut_low":32., # TOF12 cut lower bound
            "tof12_cut_high":39., # TOF12 cut upper bound
            "delta_tof01_lower":-1., # Delta TOF01 cut lower bound 
            "delta_tof01_upper":1.5, # Delta TOF01 cut upper bound 
            "delta_tof12_lower":-5., # Delta TOF01 cut lower bound 
            "delta_tof12_upper":5., # Delta TOF01 cut upper bound 
            "tof01_cut_low":tof01_min_max[0], # TOF01 cut lower bound
            "tof01_cut_high":tof01_min_max[1], # TOF01 cut upper bound
            "tof01_tramline_lower":-15.+70, # p_tof01 - p_tku
            "tof01_tramline_upper":+15.+70, # p_tof01 - p_tku 
            "p_bins":p_bins, # set of momentum bins; for now really it is just a lower and upper bound
            "p_bins_alt":[[100, 180]], # alternative momentum cut
            "p_tot_ds_low":tkd_cut[0], # downstream momentum cut lower bound
            "p_tot_ds_high":tkd_cut[1], # downstream momentum cut upper bound
            "reco_files":rogers_mc_file_names(datasets), # list of strings to be handed to glob
            "name":name, # appears on plots
            "color":4, # not used
            "pid":-13, # assume pid of tracks following TOF cut
            "pvalue_threshold":0.02, # minimum allowed pvalue for pvalue cut
            "tku_chi2_threshold":8.0, # maximum allowed chi2/dof for chi2 cut
            "tkd_chi2_threshold":4.0, # BUGGGGGGGGG # maximum allowed chi2/dof for chi2 cut
            "amplitude_corrections":None,
            "amplitude_systematics":{},
            "field_uncertainty":0.02,
            "amplitude_chi2":False,
            "amplitude_mc":True,
            "weight_tof01_source":"output/2017-02_mc_weighting/plots_3-140-empty/tof01_weights",
            "weight_tof01_target":"output/2017-02_reco_weighting/plots_3-140-10069/tof01_weights",
            "weight_tof01_mode":"sample_using_distribution",
            "amplitude_systematic_reference":None,
            "amplitude_inefficiency_cutoff":50.,
            "amplitude_recon_systematic_sources":{},
            "amplitude_performance_systematic_sources":{},
            "tku_fiducial_radius":150.,
            "tkd_fiducial_radius":150.,
            "cov_fixed_us":None, #cov_us,
            "cov_fixed_ds":None, #cov_ds,
            "amplitude_algorithm":"binned",

            "fractional_emittance_mc":True,
            "fractional_emittance_corrections":False,
            "fractional_emittance_systematics":{},
            "fractional_emittance_corrections_draw":True,
            "fractional_emittance_systematics_draw":True,

            "density_mc":True,                  # True if pure Monte Carlo
            "density_corrections_cutoff":.5,    # Cutoff above which correction is averaged
            "density_corrections":None,         # Path to the correction file
            "density_systematics":{},           # Density systematic uncertainty dictionary
            "density_corrections_draw":True,    # True if density correctoins are to be drawn
            "density_systematics_draw":True,    # True if density systematics are to be drawn
            "density_sections":False,           # True if density sections are to be printed

            "do_magnet_alignment":False,
            "do_efficiency":False, #True,
            "do_fractional_emittance":False, #True,
            "do_amplitude":False, #True,
            "do_density":True,
            "do_extrapolation":False,
            "do_globals":False, #do_globals,
            "do_mc":False, #True,
            "do_plots":False, #True,
            "do_cuts_plots":False, #True,
            "do_tof01_weighting":False,
            "do_optics":False,
            "do_data_recorder":False,
        }

class Config(object):
    # location to which data and plots will be dumped following analysis
    info_file = "geometry_08681/Maus_Information.gdml"
    will_require_tof1 = False # require at least one TOF1 Space point to even load the data
    will_require_tof2 = False # require at least one TOF2 Space point to even load the data
    tk_station = 1 # load data from a particular tracker station
    tk_plane = 0
    # prerequisite for space point cut
    will_require_triplets = False #True # require triplet space points
    upstream_cuts = { # Set to true to make data_plotter and amplitude_analysis use these cuts; False to ignore the cut
          "any_cut":None,
          "scifi_space_clusters":False,
          "scifi_space_points":False,
          "scifi_tracks_us":True,
          "scifi_nan_us":False,
          "scifi_track_points_us":False,
          "scifi_fiducial_us":True,
          "aperture_us":False,
          "pvalue_us":False,
          "chi2_us":True,
          "aperture_ds":False,
          "scifi_tracks_ds":False,
          "scifi_nan_ds":False,
          "scifi_track_points_ds":False,
          "scifi_fiducial_ds":False,
          "pvalue_ds":False,
          "chi2_ds":False,
          "tof01":False,
          "tof01_tramlines":False,
          "tof12":False,
          "p_tot_us":True,
          "p_tot_us_alt":False,
          "p_tot_ds":False,
          "tof_0_sp":False,
          "tof_1_sp":False,
          "tof_2_sp":False,
          "upstream_aperture_cut":False,
          "downstream_aperture_cut":False,
          "delta_tof01":False, #extrapolatedtof01 compared to recon tof01
          "delta_tof12":False, #extrapolatedtof12 compared to recon tof12
          "global_through_tof0":False,
          "global_through_tof1":False,
          "global_through_us_apertures":False,
          "global_through_tku_tp":False,
          "global_through_tkd_tp":False,
          "global_through_tof2":False,
          "tof01_selection":False,
          "mc_p_us":False,
          "mc_stations_us":False,
          "mc_scifi_fiducial_us":False,
          "mc_p_ds":False,
          "mc_stations_ds":False,
          "mc_scifi_fiducial_ds":False,
    }
    data_recorder_cuts = copy.deepcopy(upstream_cuts)
    downstream_cuts = copy.deepcopy(upstream_cuts)
    downstream_cuts["p_tot_ds"] = True
    downstream_cuts["tof2_sp"] = False
    downstream_cuts["pvalue_ds"] = False
    downstream_cuts["chi2_ds"] = True
    downstream_cuts["scifi_fiducial_ds"] = True
    downstream_cuts["tof12"] = False # if TOF12 is out of range chuck it (but ignore "no TOF2" events)
    downstream_cuts["scifi_tracks_ds"] = True
    extrapolation_cuts = copy.deepcopy(downstream_cuts)
    extrapolation_cuts["downstream_aperture_cut"] = True
    extrapolation_cuts["tof2_sp"] = True
    extrapolation_cuts["global_through_tkd_tp"] = True
    extrapolation_cuts["global_through_tof2"] = True
    mc_true_us_cuts = copy.deepcopy(upstream_cuts)
    mc_true_us_cuts["mc_stations_us"] = True
    mc_true_us_cuts["mc_scifi_fiducial_us"] = True
    mc_true_ds_cuts = copy.deepcopy(mc_true_us_cuts)
    mc_true_ds_cuts["mc_stations_ds"] = True
    mc_true_ds_cuts["mc_scifi_fiducial_ds"] = True
    mc_true_ds_cuts["mc_p_ds"] = True
    cut_report = [[], [], [], [], [], []]

    cut_report[0]  = ["hline", "all events", "hline",]
    cut_report[0] += ["scifi_tracks_us", "chi2_us", "scifi_fiducial_us", "hline",]
    cut_report[0] += ["p_tot_us", "hline",]
    cut_report[0] += ["upstream_cut", "hline",]

    cut_report[1] += ["hline", "upstream_cut", "hline",]
    cut_report[1] += ["scifi_tracks_ds", "chi2_ds", "scifi_fiducial_ds", "hline",]
    cut_report[1] += ["downstream_cut", "hline",]

    cut_report[2] = ["hline", "upstream_cut", "hline",]
    cut_report[2] += ["mc_stations_us", "mc_scifi_fiducial_us",]
    cut_report[2] += ["hline", "mc_true_us_cut", "hline"]

    cut_report[3] = ["hline", "mc_true_us_cut", "hline",]
    cut_report[3] += ["mc_stations_ds", "mc_scifi_fiducial_ds"]
    cut_report[3] += ["hline", "mc_true_ds_cut", "hline"]

    cut_report[4] = ["hline", "mc_true_ds_cut", "hline",]
    cut_report[4] += ["scifi_tracks_ds", "chi2_ds", "scifi_fiducial_ds", "p_tot_ds",]

    cut_report[5] = ["hline", "downstream_cut", "hline",]
    cut_report[5] += ["mc_stations_ds", "mc_scifi_fiducial_ds", "mc_p_ds",]

    data_dir = "output/2017-02-7-Systematics-test"
    analyses = []


    files = "000?"
    lih_systematics_list = [
      "mc_base", "mc_lih_plus"
    ]
    lh2_systematics_list = [
      "mc_base", "mc_beam_offset_plus", "mc_lh2_plus",
      "mc_ssu_match_plus", "mc_beam_offset_minus", "mc_fc_plus",
      "mc_ssd_match_plus",
    ]
    empty_systematics_list = [
      "tku_base", "tku_full-p",
      "tku_pos_plus", "tku_rot_plus", "tku_density_plus",
      "tku_scale_SSUC_neg", "tku_scale_SSUC_plus", "tku_scale_SSUE1_plus", "tku_scale_SSUE2_plus",
      "tkd_pos_plus", "tkd_rot_plus", "tkd_density_plus",
      "tkd_scale_SSDC_neg", "tkd_scale_SSDC_plus", "tkd_scale_SSDE1_plus", "tkd_scale_SSDE2_plus",
    ]
    cuts = {}
    suffix = None
    vers = "v8"
    run_list = [
        ["10", "10052", [1.5, 4.5], vers, "lH2 empty", suffix, cuts, files, empty_systematics_list],
        ["4",  "10064", [1.5, 6.0], vers, "lH2 empty", suffix, cuts, files, empty_systematics_list],
        ["6",  "10051", [1.5, 5.5], vers, "lH2 empty", suffix, cuts, files, empty_systematics_list],
        ["4",  "9962", [1.5, 6.0], vers, "lH2 full", suffix, cuts, files, lh2_systematics_list],
        ["6",  "9966", [1.5, 5.5], vers, "lH2 full", suffix, cuts, files, lh2_systematics_list],
        ["10", "9970", [1.5, 4.5], vers, "lH2 full", suffix, cuts, files, lh2_systematics_list],
        ["4",  "10484", [1.5, 6.0], vers, "LiH", suffix, cuts, files, lih_systematics_list],
        ["6",  "10485", [1.5, 5.5], vers, "LiH", suffix, cuts, files, lih_systematics_list],
        ["10", "10486", [1.5, 4.5], vers, "LiH", suffix, cuts, files, lih_systematics_list],
        #["3",  "10483", [1.5, 6.5], vers, "LiH", suffix, cuts, files, lih_systematics_list],
        #["3",  "9971", [1.5, 6.5], vers, "lH2 full", suffix, cuts, files, lh2_systematics_list],
        #["3",  "10069", [1.5, 6.5], vers, "lH2 empty", suffix, cuts, files, empty_systematics_list],
    ]
    tmp_run_list = []
    for a_run in run_list:
        for i in range(10):
            files = "*"+str(i)
            if "tku_base" not in a_run[-1]:
                continue
            sys_list = ["tku_base"]
            new_run = copy.deepcopy(a_run)
            new_run[-4] = str(i)
            new_run[-2] = files
            new_run[-1] = sys_list
            tmp_run_list.append(new_run)
    selection_bias = {
        "tkd_chi2_threshold":4.3,
        "tkd_fiducial_radius":148., # r2 ~ pt/bz/c_light ~ pt [mm]
    }
    for a_run in run_list:
        if "lH2 empty" not in a_run:
            continue
        sys_list = ["tku_base"]
        for key, value in selection_bias.iteritems():
            cuts = {key:value}
            new_run = copy.deepcopy(a_run)
            new_run[-1] = sys_list
            new_run[-3] = cuts
            new_run[-4] = str(key)
            tmp_run_list.append(new_run)
    run_list += tmp_run_list

    for run_items in run_list:
        [emit, run, tof, vers, absorber, suffix, cuts, files, systematics_list] = run_items
        for sys in systematics_list:
            name = str(sys)
            if suffix != None:
                name += "_"+suffix
            my_analysis = get_analysis(run+"_systematics_"+vers+"/"+sys+"/"+files,
                            "Simulated 2017-2.7 "+emit+"-140 "+absorber+" Systematics "+name, 
                            tof, data_dir, [[135, 145]], [90, 170], False)
            for cut, value in cuts.iteritems():
                my_analysis[cut] = value
            analyses.append(my_analysis)

    print "Planned", len(analyses), "analyses"
    amplitude_bin_width = 5
    amplitude_max = 25

    required_trackers = [0, 1] # for space points
    required_number_of_track_points = 12 # doesnt do anything
    global_min_step_size = 1. # for extrapolation, set the extrapolation step size
    global_max_step_size = 100. # for extrapolation, set the extrapolation step size
    will_load_tk_space_points = True # determines whether data loader will attempt to load tracker space points
    will_load_tk_track_points = True # determines whether data loader will attempt to load tracker track points
    number_of_spills = None #100 # if set to an integer, limits the number of spills loaded for each sub-analysis
    # NOTE on memory usage: looks like about 10% + 5-10% per 200 spills with 16 GB RAM
    preanalysis_number_of_spills = 50 # number of spills to analyse during "pre-analysis"
    analysis_number_of_spills = 50 # number of spills to analyse during each "analysis" step
    momentum_from_tracker = True # i.e. not from TOFs
    time_from = "tof1"

    residuals_plots_nbins = 100 # used for track extrapolation plots
    extrapolation_does_apertures = True # set to True in order to include apertures in track extrapolation
    maus_verbose_level = 5

    amplitude_bin_width = 5
    amplitude_min_events = 100
    amplitude_min_bin = 0

    fractional_emittance_bins = [0., 5., 10., 15., 20., 30., 50.]
    fractional_emittance_fraction = 0.09
    fractional_emittance_uncertainty = 0

    density_nthreads = 1
    density_knn_rotate = True # rotate to eigenvector system
    density_uncertainty = False # assume Gaussian for errors; True - use subsampling for errors
    density_npoints = 100
    density_graph_scaling = 1e9

    magnet_alignment = {
        "n_events":10,
        "max_iterations":100,
        "resolution":1.,
    }

    tof0_offset = 0.18+25.4
    tof1_offset = 0.
    tof2_offset = -27.7

    # z position of central absorber (used for offsetting
    z_apertures = 0.
    # z position of detectors (used for track extrapolation) (z, name)
    # 13692: diffuser us edge
    # 13782: diffuser ds edge
    # 16639: lH2 window mount us (us edge)
    # 16710: lH2 window mount us (ds edge)
    # 16953.5: LiH centre
    # 17166: lH2 window mount ds (us edge)
    # 17242: lH2 window mount ds (ds edge)
    # 17585: SSU aperture
    # 18733: SSU He window
    tkd_pos = 18836.8+8.
    detectors = [
        (5287.2, None, "tof0"),
        (12929.6-25., None, "tof1_us"),
        (12929.6, None, "tof1"),
        (12929.6+25., None, "tof1_ds"),
        (13968.0, None, "tku_5"),
        (14318.0, None, "tku_4"),
        (14618.0, None, "tku_3"),
        (14867.0, None, "tku_2"),
        (15068.0, None, "tku_tp"),
        (tkd_pos+0., None, "tkd_tp"),
        (tkd_pos+200., None, "tkd_2"),
        (tkd_pos+450., None, "tkd_3"),
        (tkd_pos+750., None, "tkd_4"),
        (tkd_pos+1100., None, "tkd_5"),
        (21114.4, None, "tof2"),
        (21139.4, None, "tof2"),
        (21159.4, None, "tof2"),
        (21208., None, "tof2"),
        (21214.4, None, "tof2"),
        (21220.4, None, "tof2"),

    ]
    # Note this is disabled for hybrid MC as we don't go upstream of tracker
    upstream_aperture_cut = {
        "global_through_virtual_diffuser_us":90.,
        "global_through_virtual_diffuser_ds":90.,
    }
    downstream_aperture_cut = {
        "global_through_virtual_lh2_us_window_flange_1":100.,
        "global_through_virtual_lh2_us_window_flange_2":100.,
        "global_through_virtual_absorber_centre":100.,
        "global_through_virtual_lh2_ds_window_flange_1":100.,
        "global_through_virtual_lh2_ds_window_flange_2":100.,
        "global_through_virtual_ssd_aperture":100.,
        "global_through_virtual_ssd_he":100.,
        "global_ds_virtual_ds_pry":100.,
    }

    virtual_detectors = [
        (12904., None, "virtual_tof1_us"),
        (12929., None, "virtual_tof1"),
        (12954., None, "virtual_tof1_ds"),
        (13205., None, "virtual_us_pry_us"),
        (13305., None, "virtual_us_pry_ds"),
        (13403., None, "virtual_tku_butterfly_us"),
        (13486., None, "virtual_tku_butterfly_ds"),
        (13578., None, "virtual_diffuser_shield"),
        (13698., None, "virtual_diffuser_us"),
        (13782., None, "virtual_diffuser_ds"),
        (16639., None, "virtual_lh2_us_window_flange_1"),
        (16710., None, "virtual_lh2_us_window_flange_2"),
        (16953.5, None, "virtual_absorber_centre"),
        (17166., None, "virtual_lh2_ds_window_flange_1"),
        (17242., None, "virtual_lh2_ds_window_flange_2"),
        (17585., None, "virtual_ssd_aperture"),
        (18733., None, "virtual_ssd_he"),
        (20447., None, "virtual_ds_pry"),
    ]
    plot_virtual_stations = [station for dummy, dummy, station in virtual_detectors]
    virtual_detectors += [(500.*i, None, "virtual_"+str(i)) for i in range(51)]
    virtual_detectors += [(z, None, "virtual_"+det) for z, dummy, det in detectors]
    virtual_detectors = sorted(virtual_detectors)
    #for z, dummy, plane in virtual_detectors:
    #    print z, plane

    mc_plots = { # Used for virtual_cuts as well as plots
        "mc_stations" : { # one virtual plane for each tracker view; first must be tracker reference plane
            "tku_tp":["mc_virtual_tku_tp", "mc_virtual_tku_2", "mc_virtual_tku_3", "mc_virtual_tku_4", "mc_virtual_tku_5",],
            "tkd_tp":["mc_virtual_tkd_tp", "mc_virtual_tkd_2", "mc_virtual_tkd_3", "mc_virtual_tkd_4", "mc_virtual_tkd_5",],
            #"tof1":["mc_virtual_tof1"],
            #"tof01":["mc_virtual_tof0", "mc_virtual_tof1"],
            #"tof12":["mc_virtual_tof1", "mc_virtual_tof2"],
            #"global_through_virtual_diffuser_us":["virtual_diffuser_us"],
            #"global_through_virtual_diffuser_ds":["virtual_diffuser_ds"],
        }
    }

    bz_tku = 3e-3
    bz_tkd = -2e-3
    z_afc = 16955.74
    # z position of apertures (z, maximum radius, name)
    # Notes from Jason: 209.6 to fixed flange
    # 282.5 to movable flange
    # what about tracker butterfly us and ds?
    apertures = sorted(
      [(-3325.74+z_afc, 100., "diffuser_us"),
       (-3276.74+z_afc, 100., "diffuser_mid"),
       (-3227.74+z_afc, 100., "diffuser_ds"),]+\
      [(float(z)+z_afc, 150., "tku_"+str(z)) for z in range(-2918, -1817, 100)]+\
      [(float(z)+z_afc, 150., "tkd_"+str(z)) for z in range(1818, 2919, 100)]+\
      [(float(z)+z_afc, 200., "pipe_"+str(z)) for z in range(-1800, 1801, 100)]+\
      [(+209.6+z_afc, 160., "afc_209.6")],
    )


