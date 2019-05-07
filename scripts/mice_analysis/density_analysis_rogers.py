import tempfile
import json
import ROOT
import copy
import os
import numpy as np

from mice_analysis.density.density_data import DensityData
from mice_analysis.density.density_plotter import DensityPlotter
from mice_analysis.density.knn_density_estimator import kNNDensityEstimator
from mice_analysis.analysis_base import AnalysisBase

class DensityAnalysis(AnalysisBase):

    def __init__(self, config, config_anal, data_loader):
        """
        Initialise the DensityAnalysis class for the nonparametric density estimation analysis
        * config and config_anal are the configuration files
        * data_loader extracts the data from the MAUS ROOT files
        """
        super(DensityAnalysis, self).__init__(config, config_anal, data_loader)

        # Initialize the configuration and the data loader
        self.config = config
        self.config_anal = config_anal
        self.data_loader = data_loader

        # Initialize the density estimator parameters
        self.nthreads = self.config.density_nthreads
        self.rotate = self.config.density_knn_rotate
        self.uncertainty = self.config.density_uncertainty
        self.npoints = self.config.density_npoints
        self.density_max = 50*1e-9

        # Initialize the data containers
        self.a_dir = tempfile.mkdtemp()
        file_name = self.a_dir+"/density_data_"
        self.data_types = ("all_mc", "reco_mc", "reco")
        self.locations = ("us", "ds")
        self.data = {}
        self.density_data = {}
        for typ in self.data_types:
            self.data[typ] = {}
            self.density_data[typ] = {}
            for loc in self.locations:
                self.data[typ][loc] = DensityData(file_name+"_"+typ+"_"+loc)
                self.density_data[typ][loc] = {}

        # Calculate corrections if required
        self.calculate_corrections = self.config_anal["density_corrections"] == None

        # Initialize the systematics graphs if necessary
        if self.config_anal["density_systematics_draw"]:
            self.syst_graphs = {}
            for typ in self.data_types:
                self.syst_graphs[typ] = {}
                for loc in self.locations:
                    self.syst_graphs[typ][loc] = {}

    def birth(self):
        """
        Sets the output directory for the density plots
        Resets all the data containers and loads the uncertainties
        """
        # Create directory, clear the data containers
        self.set_plot_dir("density")
        for typ in self.data_types:
            for loc in self.locations:
                self.data[typ][loc].clear()

        self.clear_density_data()

        # Load the systematic corrections/uncertainties from the file if they are provided
        self.load_errors()

        # Load the data
        self.append_data()

        # Create the subdirectories
        try:
            os.mkdir(self.plot_dir+"/phase_space")
        except OSError:
            pass
        try:
            os.mkdir(self.plot_dir+"/corrections")
        except OSError:
            pass
        try:
            os.mkdir(self.plot_dir+"/systematics")
        except OSError:
            pass

    def process(self):
        """
        Imports the data to the density profile format
        """
        self.append_data()

    def death(self):
        """
        Called when all the data has been loaded
        """
        # First, build the profiles, store the corresponding arrays
        # If requested, draw Poincare sections of phase space in density/phase_space
        self.make_pdfs("reco")
        if self.config_anal["density_mc"]:
            self.make_pdfs("reco_mc")
            self.make_pdfs("all_mc")

        # Evaluate corrections if necessary
        # If requested, draw the corrections in density/corrections
        if self.calculate_corrections:
            self.set_corrections()
        if self.config_anal["density_corrections_draw"]:
            self.draw_corrections()

        # Apply the corrections, evaluate the systematic uncertainties
        # If requested, draw the systematics in density/systematics
        self.corrections_and_uncertainties("reco")
        if self.config_anal["density_mc"]:
            self.corrections_and_uncertainties("reco_mc")
            self.corrections_and_uncertainties("all_mc")

        for loc in ["us", "ds"]:
            print loc
            print "  RECO:", self.density_data["reco"][loc]["pdf"][0:10]
            print "  MC:  ", self.density_data["all_mc"][loc]["pdf"][0:10]
            print "  CORR:", self.density_data["reco"][loc]["corrected_pdf"][0:10]

        self.draw(self.config_anal["density_mc"])

        if self.config_anal["density_systematics_draw"]:
            self.draw_systematics()

        # Save everything to a json file
        self.save()

    def append_data(self):
        """
        Add data to the fractional amplitude calculation (done at death time)
        
        If amplitude_mc is false, we take 'tku' data from upstream_cut sample 
        and 'tkd' data from downstream_cut sample
        
        If density_mc is true, then we build a couple of additional samples
        1. all_mc is for MC truth of all events that should have been included 
           in the sample (some of which might have been missed by recon; some in
           recon sample maybe should have been excluded)
        We use "mc_true_us_cut"/"mc_true_ds_cut" for the all_mc sample
        """
        hits = {}
        for typ in self.data_types:
            hits[typ] = {}
            for loc in self.locations:
                hits[typ][loc] = []

        if self.config_anal["density_mc"]:
            station_us = self.config.mc_plots["mc_stations"]["tku_tp"][0]
            station_ds = self.config.mc_plots["mc_stations"]["tkd_tp"][0]

        for event in self.data_loader.events:
            if event['upstream_cut']:
                continue
            hits["reco"]["us"].append(event['tku'])
            if not event['downstream_cut']:
                hits["reco"]["ds"].append(event['tkd'])

            if self.config_anal["density_mc"]:
                hit_mc_us, hit_mc_ds = None, None
                for detector_hit in event["data"]:
                    if detector_hit["detector"] == station_us:
                        hit_mc_us = detector_hit["hit"]
                    if detector_hit["detector"] == station_ds:
                        hit_mc_ds = detector_hit["hit"]
                if not event['mc_true_us_cut']:
                    hits["all_mc"]["us"].append(hit_mc_us)# no inefficiency upstream
                    hits["reco_mc"]["us"].append(hit_mc_us)
                if not event['mc_true_ds_cut']:
                    hits["all_mc"]["ds"].append(hit_mc_ds)
                    if not event['downstream_cut']:
                        hits["reco_mc"]["ds"].append(hit_mc_ds)

        for typ in self.data_types:
            for loc in self.locations:
                self.data[typ][loc].append_hits(hits[typ][loc])

        for loc in self.locations:
            print "Loaded %s (density):" % loc
            for typ in self.data_types:
                print "    %s:   " % typ, len(hits[typ][loc])

        return

    def make_pdfs(self, typ):
        """
        Make pdfs, i.e. number of events vs density
        
        Pdf is always self.npoints long; the last bin is an overflow bin.
        """
        if not self.data[typ]["us"].get_n_events():
            return

        transmission = float(self.data[typ]["ds"].get_n_events())/self.data[typ]["us"].get_n_events()
        for loc in self.locations:
            # Extract the data as a numpy array
            n_events = self.data[typ][loc].get_n_events()
            id_data = [None]*n_events
            ps_data = np.ndarray((n_events, 4))
            for i, event in enumerate(self.data[typ][loc].retrieve()):
                ps_data[i] = event[-1].tolist()
                id_data[i] = (event[0], event[1], event[2])

            # Use the kNN module to get a density profile
            norm = 1.
            if loc == "ds":
                norm = transmission
            density_estimator = kNNDensityEstimator(ps_data, self.rotate, self.nthreads, norm)
            density_estimator.set_levels()
            levels = density_estimator.levels
            key_value_pairs = [(id_data[i], levels[i]) for i in range(n_events)]
            density_dict = dict(key_value_pairs)
            self.density_data[typ][loc]["densities"] = density_dict
            pdf, bx = np.histogram(density_estimator.levels, self.npoints-1, (0., self.density_max))
            pdf = pdf.tolist()
            pdf.append(len(levels)-sum(pdf))
            self.density_data[typ][loc]["pdf"] = pdf
            print "Found density pdf with max entry", max(density_dict), \
                  "max bin", self.density_max, "n(overflow)", pdf[-1], \
                  "and", len(density_dict), "=", sum(pdf), "entries"
            self.density_data[typ][loc]["cdf"] = self.make_cdf(pdf)

    @classmethod
    def make_cdf(cls, pdf):
        cdf = [0. for prob in pdf]
        count = 0.
        for i, n_events in enumerate(reversed(pdf)):
            count += n_events
            cdf[i] = count
        return cdf

    def corrections_and_uncertainties(self, typ):
        """
        Calculate corrected profiles and uncertainties
        * typ specifies the type of data (all_mc, reco)
        """
        # Set the upstream statistical uncertainty to 0 as it is the given
        # profile to which to compare the downstream profile. The downstream
        # statistical uncertainties are set in self.make_profiles()
        data = self.density_data[typ]
        # Do the correction for each of the tracker locations
        for loc in self.locations:
            print "Doing density level correction for", typ, loc
            corr_pdf = self.do_corrections(data, typ, loc, data)
            corr_cdf = self.make_cdf(corr_pdf)
            data[loc]["corrected_pdf"] = corr_pdf
            data[loc]["corrected_cdf"] = corr_cdf
        self.density_data[typ] = data
        return

        data["us"]["levels_stat_errors"] = [0. for bin in range(self.npoints)]
        # Evaluate the systematic uncertainties for each of the tracker locations
        for loc in self.locations:
            print "Finding density systematic errors for", typ, loc
            reco_syst_list = self.calculate_detector_systematics(typ, loc)
            perf_syst_list = self.calculate_performance_systematics(typ, loc)
            syst_error_list = [(reco_syst_list[i]**2+perf_syst_list[i]**2)**0.5 \
                                                  for i in range(self.npoints)]
            data[loc]["levels_syst_errors"] = syst_error_list

        self.density_data[typ] = data

    def do_corrections(self, ref_source, typ, loc, source):
        """
        Applies the corrections to the requested density profile
        Only applies response correction to the reconstructed sample
        * typ specifies the type of data (all_mc, reco)
        * loc specifies the location of the tracker (us, ds)
        * source specifies the source of the corrections to be used
        * Use capped corrections if use_capped is True
        """
        print "Doing corrections", typ, loc
        pdf = source[loc]["pdf"]
        mig = ref_source[loc]["migration_matrix"]
        eff = ref_source[loc]["pdf_ratio"]
        corr_pdf = [0. for i in range(self.npoints)]
        for i in range(self.npoints):
            for j in range(self.npoints):
                corr_pdf[i] += mig[i][j]*pdf[j]*eff[i]
                if i < 5 and j < 5:
                    print "  correction", i, j, ":", mig[i][j], "*", pdf[j], "=", corr_pdf[i]
            if i < 5:
                print
        return corr_pdf

    def calculate_detector_systematics(self, typ, loc):
        """
        Calculate the systematic errors in the reconstruction of tracks. The uncertainty
        on each level corresponds to the sum in quadrature of the residuals between the reference
        reconstruction set and the data sets that are shifted from the reference set
        * typ specifies the type of data (all_mc, reco)
        * loc specifies the location of the tracker (us, ds)
        """
        # If there is no reference specified, skip
        data = self.density_data[typ]
        syst_error_list = [0. for i in range(self.npoints)]
        if data["detector_reference"] == None:
            return syst_error_list

        print "\nEvaluating density reconstruction systematic errors", loc

        # Correct the density profile with the reference corrections
        ref_source = data["detector_reference"]
        ref_levels = self.do_corrections(ref_source, typ, loc, ref_source)
        print "Reference density sytematic error"
        print "   ", ref_levels[:5], "...", ref_levels[-5:]
        response = ref_source["response"][loc]["level_ratio_capped"]
        print "   ", response[:5], "...", response[-5:]
        # Loop over the detector systematics list
        systematics_list = data[loc]["detector_systematics"]
        for source in systematics_list:
            # Evaluate the levels with the corresponding systematic shift
            syst_levels = self.do_corrections(ref_source, typ, loc, source)
            print "Systematic density sytematic error for", source['source']
            print "   sys levels:", syst_levels[:5], "...", syst_levels[-5:]
            response = source["response"][loc]["level_ratio_capped"]
            print "   response:  ", response[:5], "...", response[-5:]
            inefficiency = source["inefficiency"][loc]["level_ratio_capped"]
            print "   inefficncy:", inefficiency[:5], "...", inefficiency[-5:]
            err_levels = source[typ][loc]["levels"]
            print "   err levels:", err_levels[:5], "...", err_levels[-5:]
            # Initialize a graph that contains the deviation from the reference
            name = self.get_syst_name(source["source"])
            if self.config_anal["density_systematics_draw"]:
                self.syst_graphs[typ][loc][name] = ROOT.TGraph(self.npoints)

            # Add in quadrature an uncertainty that corresponds to the level shift due
            # to the use of a different set of corrections
            scale = source["scale"]
            for j in range(self.npoints):
                err = (syst_levels[j] - ref_levels[j])*scale
                syst_error_list[j] = (syst_error_list[j]**2+err**2)**0.5

                if self.config_anal["density_systematics_draw"]:
                    alpha = (float(j+1.)/(self.npoints+1.))
                    val = 0.
                    if ref_levels[j] > 0:
                        val = err/ref_levels[j]
                    self.syst_graphs[typ][loc][name].SetPoint(j, alpha, val)
            print

        return syst_error_list

    def calculate_performance_systematics(self, typ, loc):
        """
        Calculate the systematic errors in the channel performance. The experiment measures
        ratios of downstream density levels over upstream density levels.
        The uncertainty is evaluated as the shifts in the downstream density profile for
        variable deviations from the expected cooling channel performance.
        * typ specifies the type of data (all_mc, reco)
        * loc specifies the location of the tracker (us, ds)
        """
        # If there is no reference specified, skip
        data = self.density_data[typ]
        syst_error_list = [0. for i in range(self.npoints)]
        if data["performance_reference"] == None:
            return syst_error_list

        print "\nEvaluating density performance systematic errors", loc

        # Get the reference ratio array
        source = data["performance_reference"]
        ref_ratio = np.array(source[typ]["ds"]["levels"])/np.array(source[typ]["us"]["levels"])
        ref_ratio = ref_ratio.tolist()

        # Loop over the performance systematics list
        systematics_list = data[loc]["performance_systematics"]
        ratio_error_list = [0. for i in range(self.npoints)]
        for source in systematics_list:
            # Evaluate the ratio with the corresponding systematic shift
            syst_ratio = np.array(source[typ]["ds"]["levels"])/np.array(source[typ]["us"]["levels"])
            syst_ratio = syst_ratio.tolist()

            # Initialize a graph that contains the deviation from the reference
            name = self.get_syst_name(source["source"])
            if self.config_anal["density_systematics_draw"]:
                self.syst_graphs[typ][loc][name] = ROOT.TGraph(self.npoints)

            # Add in quadrature an uncertainty that corresponds to the ratio shift due
            # to the use of a different cooling channel
            scale = source["scale"]
            for j in range(self.npoints):
                err = (syst_ratio[j] - ref_ratio[j])*scale
                ratio_error_list[j] = (ratio_error_list[j]**2+err**2)**0.5

                if self.config_anal["density_systematics_draw"]:
                    alpha = (float(j+1.)/(self.npoints+1.))
                    self.syst_graphs[typ][loc][name].SetPoint(j, alpha, err)

        # Convert the uncertainties in terms of density
        ref_levels_us = data["us"]["levels"]
        for i in range(self.npoints):
            syst_error_list[i] = ratio_error_list[i] * ref_levels_us[i]

        return syst_error_list

    def get_syst_name(self, path):
        """
        Convert systematic path to a systematic name
        """
        suffix = path.split("Systematics_",1)[1]
        name = suffix.split("/")[0]
        return name

    def draw_systematics(self):
        """
        Draws the systematic errors. The uncertainty on each level corresponds to the 
        residuals between the reference reconstruction set and the data sets that 
        are shifted from the reference set
        """
        # Feed the systematics graphs to the drawer
        for typ in self.data_types:
            for loc in self.locations:
                if len(self.syst_graphs[typ][loc]):
                    plotter = DensityPlotter(self.plot_dir, typ+"_"+loc)
                    plotter.plot_systematics(self.syst_graphs[typ][loc])

    def draw(self, is_mc):
        """
        Produce plots that compare the density profiles upstream and downstream
        of the absorber. Produce one for each category of data
        """
        # Initialize the graphs
        graphs = {}
        y_norm = 1./self.density_data["reco"]["us"]["corrected_cdf"][-1]
        for loc in self.locations:
            if is_mc:
                graphs[loc+"_all_mc_pdf"] = self.make_graph("all_mc", loc, "pdf", None, 1e9, y_norm)
            graphs[loc+"_reco_pdf"] = self.make_graph("reco", loc, "pdf", None, 1e9, y_norm)
            graphs[loc+"_reco_corrected_pdf"] = self.make_graph("reco", loc, "corrected_pdf", None, 1e9, y_norm)

        # Print up/down comparison
        canvas_name = 'density_pdf_correction'
        canvas = self.get_plot(canvas_name)["pad"]
        mg = self.make_multigraph(graphs, "density_pdf_correction")
        leg = self.make_multigraph_legend(graphs, [.6, .65, .8, .85])
        for fmt in ["pdf", "png", "root"]:
            canvas.Print(self.plot_dir+"/"+canvas_name+"."+fmt)

        # Initialize the graphs
        graphs = {}
        for loc in self.locations:
            if is_mc:
                graphs[loc+"_all_mc_cdf"] = self.make_graph("all_mc", loc, "cdf", None, 1e9, y_norm)
            graphs[loc+"_reco_cdf"] = self.make_graph("reco", loc, "cdf", None, 1e9, y_norm)
            graphs[loc+"_reco_corrected_cdf"] = self.make_graph("reco", loc, "corrected_cdf", None, 1e9, y_norm)

        # Print up/down comparison
        canvas_name = 'density_cdf_correction'
        canvas = self.get_plot(canvas_name)["pad"]
        mg = self.make_multigraph(graphs, "density_cdf_correction")
        leg = self.make_multigraph_legend(graphs, [.2, .65, .4, .85])
        for fmt in ["pdf", "png", "root"]:
            canvas.Print(self.plot_dir+"/"+canvas_name+"."+fmt)


    def make_graph(self, typ, loc, level_type, error_type, x_norm, y_norm):
        """
        Builds a TGraphErrors for the requested data type and location
        * typ specifies the type of data (all_mc, reco)
        * loc specifies the location of the tracker (us, ds)
        * include_corr is True if the corrected levels are to be represented
        * include_syst is True if the systematic uncertainty is to be includes
        """
        scale_factor = self.config.density_graph_scaling
        print "Plotting density", loc, typ
        graph = ROOT.TGraphErrors(self.npoints-1)
        x_points = [self.density_max*i/(self.npoints-1) for i in range(self.npoints-1)]
        y_points = self.density_data[typ][loc][level_type]
        x_points = [x*x_norm for x in x_points]
        print y_points, y_norm
        y_points = [y*y_norm for y in y_points]
        for i in range(self.npoints-1): # not overflow bin
            graph.SetPoint(i, x_points[i], y_points[i])
            if error_type and value > 0.:
                err = self.density_data[typ][loc][error_type][i]
                graph.SetPointError(i, 0., err)

        color = {"us":1,"ds":4}[loc]
        style = 2 # dashed line if we don't have correction or plot mc truth
        if typ == "reco" and "corrected" in level_type:
            style = 1 # full line for the corrected reco i.e. the basic plot
        width = 1
        if typ == "all_mc": # thicker if we do all_mc
            width = 3
        graph.SetLineWidth(width)
        graph.SetLineStyle(style)
        graph.SetLineColor(color)
        graph.SetFillColorAlpha(color, .25)

        return graph

    def make_multigraph(self, graphs, name):
        """
        Initializes a multigraph, draws it
        * typ specifies the type of data (all_mc, reco)
        * graphs is a dictionary containing the up and downstream graphs with stat errors
        * graphs_full is a dictionary containing the up and downstream graphs with full errors
        * name of the type of graph being output
        """

        mg = ROOT.TMultiGraph(name, ";#rho_{#alpha} [mm^{-2}(MeV/c)^{-2}];Fraction of upstream sample")
        for key in graphs:
            mg.Add(graphs[key], "LE3")
            self.plots[name]["graphs"][key] = graphs[key]

        mg.Draw("A")
        return mg

    def make_multigraph_legend(self, graphs, pos):
        """
        Initializes a multigraph legend, draws it
        * graphs is a dictionary containing the up and downstream graphs
        """

        leg = ROOT.TLegend(pos[0], pos[1], pos[2], pos[3])
        for key in sorted(graphs.keys()):
            leg.AddEntry(graphs[key], key, "LF")
        leg.Draw("SAME")
        return leg

    def make_ratio(self, typ, graphs, graphs_full, name):
        """
        Initializes a graph ratio, draws it
        * typ specifies the type of data (all_mc, reco)
        * graphs is a dictionary containing the up and downstream graphs with stat errors
        * graphs_full is a dictionary containing the up and downstream graphs with full errors
        * name of the type of graph being output
        """
        gratio_multi = ROOT.TMultiGraph()
        gratio = ROOT.TGraphErrors(self.npoints)
        gratio_full = ROOT.TGraphErrors(self.npoints)
        gratio_multi.SetTitle(";Fraction #alpha;#rho_{#alpha}^{d} /#rho_{#alpha}^{u}")
        for i in range(self.npoints):
            us, ds = graphs["us"].GetY()[i], graphs["ds"].GetY()[i]
            use, dse = graphs["us"].GetEY()[i], graphs["ds"].GetEY()[i]
            ratio = ds/us
            gratio.GetX()[i] = graphs["us"].GetX()[i]
            gratio.GetEX()[i] = graphs["us"].GetEX()[i]
            gratio.GetY()[i] = ratio
            gratio.GetEY()[i] = dse/us

            gratio_full.GetX()[i] = gratio.GetX()[i]
            gratio_full.GetEX()[i] = gratio.GetEX()[i]
            gratio_full.GetY()[i] = ratio
            gratio_full.GetEY()[i] = 0.
            us_rel_err = dse/us
            ds_rel_err = ratio*use/us
            gratio_full.GetEY()[i] = (us_rel_err**2 + ds_rel_err**2)**0.5

        self.plots[name+"_"+typ]["graphs"]["ratio"] = gratio
        self.plots[name+"_"+typ]["graphs"]["ratio_full"] = gratio_full

        gratio.SetLineColor(1)
        gratio.SetFillColorAlpha(1, .25)
        gratio.SetName("stats error")
        gratio_full.SetLineColor(1)
        gratio_full.SetFillColorAlpha(1, .25)
        gratio_full.SetName("sys error")

        gratio_multi.Add(gratio)
        gratio_multi.Add(gratio_full)
        gratio_multi.SetName(name+"_"+typ)
        gratio_multi.Draw("A LE3")
        return gratio_multi

    def save(self):
        """
        Saves the data dictionary to a json file
        """
        fout = open(self.plot_dir+"/density.json", "w")
        density_data = copy.deepcopy(self.density_data)
        for typ in self.data_types:
            for loc in self.locations:
                density_data[typ][loc]["densities"] = {}
        print >> fout, json.dumps(density_data, sort_keys=True, indent=4)

    def set_corrections(self):
        """
        Calculate the profile corrections, i.e. the inefficiency and response function
        * uses the Monte Carlo to generate corrections
        """
        # I realised there is a problem; I can't calculate the migration for
        # events that do not have a reconstruction (inefficiency) because there
        # is no reconstruction ... by definition...
        # this needs more thought
        for loc in self.locations:
            print "setting correction for", loc
            # Initialize the data
            all_mc_dens = self.density_data["all_mc"][loc]["densities"]
            reco_mc_dens = self.density_data["reco_mc"][loc]["densities"]
            reco_dens = self.density_data["reco"][loc]["densities"]
            all_mc_hist = self.density_data["all_mc"][loc]["pdf"]
            reco_mc_hist = self.density_data["reco_mc"][loc]["pdf"]
            reco_hist = self.density_data["reco"][loc]["pdf"]
            # use npoints-1 because last bin is  overflow
            step = self.density_max/(self.npoints-1)
            hist_2d = [[0. for i in range(self.npoints)] for j in range(self.npoints)]
            efficiency_corr = [1. for i in range(self.npoints)]
            for key in sorted(reco_mc_dens.keys()):
                if key not in reco_mc_dens:
                    continue
                bin_1 = int(reco_mc_dens[key]/step)
                bin_2 = int(reco_dens[key]/step)
                # deal with overflow
                if bin_1 >= self.npoints:
                    bin_1 = self.npoints-1
                if bin_2 >= self.npoints:
                    bin_2 = self.npoints-1
                hist_2d[bin_1][bin_2] += 1.
            for i in range(self.npoints):
                for j in range(self.npoints):
                    if reco_hist[j] > 0.5:
                        hist_2d[i][j] = hist_2d[i][j]/reco_hist[j]
                    else:
                        hist_2d[i][j] = 0.
                    if i < 5 and j < 5:
                        print "  migration", i, j, ":", hist_2d[i][j], reco_hist[j]
                if reco_mc_hist[i] > 0.5:
                    efficiency_corr[i] = 1.0*all_mc_hist[i]/reco_mc_hist[i]

            self.density_data["reco"][loc]["migration_matrix"] = hist_2d
            self.density_data["reco"][loc]["pdf_ratio"] = efficiency_corr
            self.density_data["reco_mc"][loc]["migration_matrix"] = \
                  np.identity(self.npoints).tolist()
            self.density_data["reco_mc"][loc]["pdf_ratio"] = efficiency_corr
            self.density_data["all_mc"][loc]["migration_matrix"] = \
                  np.identity(self.npoints).tolist()
            self.density_data["all_mc"][loc]["pdf_ratio"] = \
                  [1. for i in range(self.npoints)]

    def draw_corrections(self):
        """
        Draw the correction factors used
        """
        for loc in self.locations:
            inefficiency = self.density_data["inefficiency"][loc]["level_ratio"]
            response = self.density_data["response"][loc]["level_ratio"]
            plotter = DensityPlotter(self.plot_dir, loc)
            plotter.plot_corrections(inefficiency, response)

            inefficiency = self.density_data["inefficiency"][loc]["level_ratio_capped"]
            response = self.density_data["response"][loc]["level_ratio_capped"]
            plotter = DensityPlotter(self.plot_dir, "capped_"+loc)
            plotter.plot_corrections(inefficiency, response)

    def clear_density_data(self):
        """
        Initializes the dictionary that is used to store
        data and make corrections down the line
        """
        self.density_data = {
          "inefficiency":{
              "us":{
                  "level_ratio":[1. for i in range(self.npoints)],
                      "level_ratio_capped":[1. for i in range(self.npoints)]
              },
              "ds":{
                  "level_ratio":[1. for i in range(self.npoints)],
                      "level_ratio_capped":[1. for i in range(self.npoints)]
              },
          },
          "response":{
              "us":{
                  "level_ratio":[1. for i in range(self.npoints)],
                      "level_ratio_capped":[1. for i in range(self.npoints)]
              },
              "ds":{
                  "level_ratio":[1. for i in range(self.npoints)],
                      "level_ratio_capped":[1. for i in range(self.npoints)]
              },
          },
          "source":"",
          "scale":1.,
          "npoints":0
        }

        level_data = {
            "performance_reference":None,
            "detector_reference":None,
            "us":{
                "pdf":[],
                "corrected_pdf":[],
                "levels_stat_errors":[],
                "levels_syst_errors":[],
                "detector_systematics":[],
                "performance_systematics":[],
            },
            "ds":{
                "pdf":[],
                "corrected_pdf":[],
                "levels_stat_errors":[],
                "levels_syst_errors":[],
                "detector_systematics":[],
                "performance_systematics":[],
            }
        }

        for typ in self.data_types:
            self.density_data[typ] = copy.deepcopy(level_data)

    def load_errors(self):
        """
        Two "classes" of systematic errors;
        * systematic errors on the reconstruction are contained in the
          correction factors. For these we store the correction factors and 
          compare to the reference correction factors
        * systematic errors on the performance are contained in the actual
          density profile. For these we store the point-by-point fractional
          difference between the density profile and reference.
        """
        # If the corrections are to calculated in this analysis, skip this step
        if self.calculate_corrections:
            return

        # Set base correction factors
        self.load_corrections(self.config_anal["density_corrections"])

        # Load systematic uncertainties
        systematics = self.config_anal["density_systematics"]
        for typ in systematics:
            print "Loading density systematic errors for", typ
            if typ not in self.density_data:
                self.density_data[typ] = {}
            for ref_key in ["detector_reference", "performance_reference"]:
                ref_src = systematics[typ][ref_key]
                if ref_src == None:
                    self.density_data[typ][ref_key] = None
                else:
                    self.density_data[typ][ref_key] = \
                                              self.load_one_error(ref_src, None)
                print "  Loaded reference", typ, ref_key, ref_src, \
                                          type(self.density_data[typ][ref_key])
            for loc in ["us", "ds"]:
                if loc not in self.density_data[typ]:
                    self.density_data[typ][loc] = {}
                for key in ["detector_systematics", "performance_systematics"]:
                    err_src_dict = systematics[typ][loc][key]
                    self.density_data[typ][loc][key] = [
                        self.load_one_error(err_src, scale) \
                              for err_src, scale in err_src_dict.iteritems()
                    ]
                    print "  Loaded", len(self.density_data[typ][loc][key]), loc, key

    def load_corrections(self, file_name):
        """
        Load the density corrections to be applied during this density
        analysis. Loads the correction factors
        """
        fin = open(file_name)
        density_str = fin.read()
        src_density = json.loads(density_str)
        src_density["source"] = file_name
        self.density_data["inefficiency"] = src_density["inefficiency"]
        self.density_data["response"] = src_density["response"]

    def load_one_error(self, file_name, scale):
        """
        Load the density analysis output for a given uncertainty source
        """
        fin = open(file_name)
        density_str = fin.read()
        density = json.loads(density_str)
        density["source"] = file_name
        density["scale"] = scale
        return density
