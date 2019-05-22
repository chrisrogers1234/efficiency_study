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

        # Initialize the data containers
        self.a_dir = tempfile.mkdtemp()
        file_name = self.a_dir+"/density_data_"
        self.data_types = ("all_mc", "reco_mc", "reco")
        self.locations = ("us", "ds")
        self.data = {}
        for typ in self.data_types:
            self.data[typ] = {}
            for loc in self.locations:
                self.data[typ][loc] = DensityData(file_name+"_"+typ+"_"+loc)

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
        self.make_profiles("reco")
        if self.config_anal["density_mc"]:
            self.make_profiles("all_mc")
            self.make_profiles("reco_mc")

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
            self.corrections_and_uncertainties("all_mc")
            self.corrections_and_uncertainties("reco_mc")
        if self.config_anal["density_systematics_draw"]:
            self.draw_systematics()

        # Draw the density profiles
        self.draw_profiles()

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
        2. reco_mc is for MC truth of all events that should have been included
           in the sample and were reconstructed
        We use "mc_true_us_cut"/"mc_true_ds_cut" for the all_mc sample and we
        use ("mc_true_us_cut" AND "upstream_cut") for the reco_mc sample
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

    def make_profiles(self, typ):
        """
        Produces density profiles. Extract a numpy ndarray that contains all of
        the phase space vectors, initialize the kNN density estimator and extract
        the density profile with its uncertainty
        * typ specifies the type of data (all_mc, reco_mc, reco)
        """
        # Skip if no data
        if not self.data[typ]["us"].get_n_events():
            return

        transmission = float(self.data[typ]["ds"].get_n_events())/self.data[typ]["us"].get_n_events()
        for loc in self.locations:
            # Extract the data as a numpy array
            ps_data = np.ndarray((self.data[typ][loc].get_n_events(), 4))
            for i, event in enumerate(self.data[typ][loc].retrieve()):
                ps_data[i] = event[-1].tolist()

            # Use the kNN module to get a density profile
            norm = 1.
            if loc == "ds":
                norm = transmission
            density_estimator = kNNDensityEstimator(ps_data, self.rotate, self.nthreads, norm)
            levels, errors = density_estimator.profile(self.npoints, self.uncertainty)

            # Store the profiles and their statistical uncertainties
            self.density_data[typ][loc]["levels"] = levels
            self.density_data[typ][loc]["levels_stat_errors"] = errors
            print "DensityAnalysis", loc
            print "   levels", levels
            print "   errors", errors

            # Plot the Poincare sections if requested
            if self.config_anal["density_sections"]:
                plotter = DensityPlotter(self.plot_dir, typ+"_"+loc)
                plotter.plot_phase_space(density_estimator)

    def corrections_and_uncertainties(self, typ):
        """
        Calculate corrected profiles and uncertainties
        * typ specifies the type of data (all_mc, reco_mc, reco)
        * corrected profiles is given by rho(alpha) = c(alpha)*i(alpha)*rho(alpha) where c(alpha)
            is the response function and i(alpha) is the inefficiency function
        * statistical errors are provided by the density estimator
        * total errors are given by the sum in quadrature of statistical errors and
          systematic errors
        """
        # Set the upstream statistical uncertainty to 0 as it is the given
        # profile to which to compare the downstream profile. The downstream
        # statistical uncertainties are set in self.make_profiles()
        data = self.density_data[typ]
        data["us"]["levels_stat_errors"] = [0. for bin in range(self.npoints)]

        # Do the correction for each of the tracker locations
        for loc in self.locations:
            print "Doing density level correction for", typ, loc
            source = self.density_data
            levels = self.do_corrections(source, typ, loc, source)
            data[loc]["corrected_levels"] = levels

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
        * typ specifies the type of data (all_mc, reco_mc, reco)
        * loc specifies the location of the tracker (us, ds)
        * source specifies the source of the corrections to be used
        """
        print "density_analysis.do_corrections", ref_source[typ][loc].keys()
        levels = np.array(ref_source[typ][loc]["levels"])
        corr_key = "level_ratio"
        if self.config_anal["density_use_capped"]:
            corr_key = "level_ratio_capped"
        if typ == "reco":
            response = np.array(source["response"][loc][corr_key])
            levels = (1.+levels)*response-1.
        if typ == "reco" or typ == "reco_mc":
            inefficiency = np.array(source["inefficiency"][loc][corr_key])
            levels = (1.+levels)*inefficiency-1.

        return levels.tolist()

    def calculate_detector_systematics(self, typ, loc):
        """
        Calculate the systematic errors in the reconstruction of tracks. The uncertainty
        on each level corresponds to the sum in quadrature of the residuals between the reference
        reconstruction set and the data sets that are shifted from the reference set
        * typ specifies the type of data (all_mc, reco_mc, reco)
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
        * typ specifies the type of data (all_mc, reco_mc, reco)
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

    def draw_profiles(self):
        """
        Produce plots that compare the density profiles upstream and downstream
        of the absorber. Produce one for each category of data
        """
        for typ in self.data_types:
            # Skip if no data
            if not len(self.density_data[typ]["us"]["levels"]):
                continue

            # Initialize the graphs
            graphs = {}
            graphs_full = {}
            for loc in self.locations:
                graphs[loc] = self.make_graph(typ, loc, True, False, False)
                graphs_full[loc] = self.make_graph(typ, loc, True, True, False)

            # Print up/down comparison
            canvas_name = 'density_profile_%s' % typ
            canvas = self.get_plot(canvas_name)["pad"]

            mg = self.make_multigraph(typ, graphs, graphs_full, "density_profile")
            leg = self.make_multigraph_legend(graphs)
            for fmt in ["pdf", "png", "root"]:
                canvas.Print(self.plot_dir+"/"+canvas_name+"."+fmt)

            # Initialize the graphs
            graphs = {}
            graphs_full = {}
            for loc in self.locations:
                graphs[loc] = self.make_graph(typ, loc, True, False, True)
                graphs_full[loc] = self.make_graph(typ, loc, True, True, True)

            # Print up/down comparison
            canvas_name = 'normalised_density_profile_%s' % typ
            canvas = self.get_plot(canvas_name)["pad"]

            mg = self.make_multigraph(typ, graphs, graphs_full, "normalised_density_profile")
            leg = self.make_multigraph_legend(graphs)
            for fmt in ["pdf", "png", "root"]:
                canvas.Print(self.plot_dir+"/"+canvas_name+"."+fmt)


            # Print ratios
            canvas_name = 'density_ratio_%s' % typ
            canvas = self.get_plot(canvas_name)["pad"]

            gratio = self.make_ratio(typ, graphs, graphs_full, "density_ratio")
            for fmt in ["pdf", "png", "root"]:
                canvas.Print(self.plot_dir+"/"+canvas_name+"."+fmt)

    def make_graph(self, typ, loc, include_corr, include_syst, normalised):
        """
        Builds a TGraphErrors for the requested data type and location
        * typ specifies the type of data (all_mc, reco_mc, reco)
        * loc specifies the location of the tracker (us, ds)
        * include_corr is True if the corrected levels are to be represented
        * include_syst is True if the systematic uncertainty is to be includes
        """
        scale_factor = self.config.density_graph_scaling
        level_type = "levels"
        if include_corr:
            level_type = "corrected_levels"
        if normalised:
            scale_factor = 1./max(self.density_data[typ]["us"][level_type])

        print "Plotting density", loc, typ
        graph = ROOT.TGraphErrors(self.npoints)
        for i in range(self.npoints):
            alpha = (float(i+1.)/(self.npoints+1.))
            value = self.density_data[typ][loc][level_type][i]*scale_factor
            graph.SetPoint(i, alpha, value)
            print "    ", i, format(value, "4.2g"),
            all_err = self.density_data[typ][loc]["levels_stat_errors"][i]*scale_factor
            print "stt:", format(all_err, "4.2g"), "sys:",
            if include_syst and value > 0.:
                syst_err = self.density_data[typ][loc]["levels_syst_errors"][i]
                if normalised:
                    syst_err *= scale_factor/self.config.density_graph_scaling
                if abs(syst_err) < 1e-6 and abs(syst_err) > 1e-12:
                    print """
                  Warning - I did something hacky; I ran the systematics calculation
                  with a 1e9 scale factor and then the main analysis without the scale factor; so
                  now there is a discrepancy. If I rerun the systematics they will come out 1e9
                  too low and I will need to scale them (e.g. *scale_factor) like everywhere else.
                  Systematic error was """+str(syst_err)
                all_err = (syst_err**2+all_err**2)**0.5
                print format(syst_err, "4.2g"),
            print "all:", format(all_err, "4.2g")
            graph.SetPointError(i, 0., all_err)

        color = {"us":1,"ds":4}[loc]
        graph.SetLineColor(color)
        graph.SetFillColorAlpha(color, .25)

        return graph

    def make_multigraph(self, typ, graphs, graphs_full, name):
        """
        Initializes a multigraph, draws it
        * typ specifies the type of data (all_mc, reco_mc, reco)
        * graphs is a dictionary containing the up and downstream graphs with stat errors
        * graphs_full is a dictionary containing the up and downstream graphs with full errors
        * name of the type of graph being output
        """

        mg = ROOT.TMultiGraph(name+"_"+typ, ";Fraction #alpha;#rho_{#alpha} [mm^{-2}(MeV/c)^{-2}]");
        for loc in self.locations:
            mg.Add(graphs[loc], "LE3")
            mg.Add(graphs_full[loc], "LE3")
            self.plots[name+"_"+typ]["graphs"][loc] = graphs[loc]
            self.plots[name+"_"+typ]["graphs"][loc+"_full"] = graphs_full[loc]

        mg.Draw("A")
        return mg

    def make_multigraph_legend(self, graphs):
        """
        Initializes a multigraph legend, draws it
        * graphs is a dictionary containing the up and downstream graphs
        """

        leg = ROOT.TLegend(.6, .65, .8, .85)
        leg.AddEntry(graphs["us"], "Upstream", "LF")
        leg.AddEntry(graphs["ds"], "Downstream", "LF")
        leg.Draw("SAME")
        return leg

    def make_ratio(self, typ, graphs, graphs_full, name):
        """
        Initializes a graph ratio, draws it
        * typ specifies the type of data (all_mc, reco_mc, reco)
        * graphs is a dictionary containing the up and downstream graphs with stat errors
        * graphs_full is a dictionary containing the up and downstream graphs with full errors
        * name of the type of graph being output
        """
        gratio_multi = ROOT.TMultiGraph()
        gratio = ROOT.TGraphErrors(self.npoints)
        gratio_full = ROOT.TGraphErrors(self.npoints)
        gratio_multi.SetTitle(";Fraction #alpha;#rho_{#alpha}^{d} /#rho_{#alpha}^{u}")
        print "Making the ratio plot for", typ
        for i in range(self.npoints):
            us, ds = graphs["us"].GetY()[i], graphs["ds"].GetY()[i]
            use, dse = graphs["us"].GetEY()[i], graphs["ds"].GetEY()[i]
            ratio = ds/us
            gratio.GetX()[i] = graphs["us"].GetX()[i]
            gratio.GetEX()[i] = graphs["us"].GetEX()[i]
            gratio.GetY()[i] = ratio
            gratio.GetEY()[i] = dse/us

            use, dse = graphs_full["us"].GetEY()[i], graphs_full["ds"].GetEY()[i]
            gratio_full.GetX()[i] = gratio.GetX()[i]
            gratio_full.GetEX()[i] = gratio.GetEX()[i]
            gratio_full.GetY()[i] = ratio
            gratio_full.GetEY()[i] = 0.
            us_rel_err = dse/us
            ds_rel_err = ratio*use/us
            gratio_full.GetEY()[i] = (us_rel_err**2 + ds_rel_err**2)**0.5
            print i, ratio, us_rel_err, ds_rel_err, gratio_full.GetEY()[i]

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
        print >> fout, json.dumps(self.density_data, sort_keys=True, indent=4)

    def set_corrections(self):
        """
        Calculate the profile corrections, i.e. the inefficiency and response function
        * uses the Monte Carlo to generate corrections
        """
        for loc in self.locations:
            # Initialize the data
            all_mc_levels = self.density_data["all_mc"][loc]["levels"]
            reco_mc_levels = self.density_data["reco_mc"][loc]["levels"]
            reco_levels = self.density_data["reco"][loc]["levels"]

            # Calculate the corrections
            inefficiency, response = [], []
            for i in range(len(all_mc_levels)):
                # Inherent detector inefficiency
                inefficiency.append(float(1.+all_mc_levels[i])/(1.+reco_mc_levels[i]))

                # Detector response function
                response.append(float(1.+reco_mc_levels[i])/(1.+reco_levels[i]))

            # Produce a capped version of the corrections
            cutoff = self.config_anal["density_corrections_cutoff"]
            cutoff_index = int(cutoff*(self.npoints+1.))
            ineff_cap = inefficiency[cutoff_index]
            resp_cap = response[cutoff_index]

            inefficiency_capped = copy.deepcopy(inefficiency)
            response_capped = copy.deepcopy(response)
            for i in range(cutoff_index, len(inefficiency_capped)):
                inefficiency_capped[i] = ineff_cap
                response_capped[i] = resp_cap

            # Store the correction factors
            self.density_data["inefficiency"][loc] = {
                "level_ratio":inefficiency,
                "level_ratio_capped":inefficiency_capped,
            }
            self.density_data["response"][loc] = {
                "level_ratio":response,
                "level_ratio_capped":response_capped,
            }

    def draw_corrections(self):
        """
        Draw the correction factors used
        """
        for loc in self.locations:
            levels = np.array(self.density_data["reco"][loc]["levels"])

            response = self.density_data["response"][loc]["level_ratio"]
            levels_corr = (1.+levels)*response-1.
            response = [levels_corr[i]/levels[i] if levels[i] > 0. else 1. for i in range(len(levels))]
            inefficiency = self.density_data["inefficiency"][loc]["level_ratio"]
            levels_corr_2 = (1.+levels_corr)*inefficiency-1.
            inefficiency = [levels_corr_2[i]/levels_corr[i] if levels[i] > 0. else 1. for i in range(len(levels))]
            plotter = DensityPlotter(self.plot_dir, loc)
            plotter.plot_corrections(inefficiency, response)

            response = self.density_data["response"][loc]["level_ratio_capped"]
            levels_corr = (1.+levels)*response-1.
            response = [levels_corr[i]/levels[i] if levels[i] > 0. else 1. for i in range(len(levels))]
            inefficiency = self.density_data["inefficiency"][loc]["level_ratio_capped"]
            levels_corr_2 = (1.+levels_corr)*inefficiency-1.
            inefficiency = [levels_corr_2[i]/levels_corr[i] if levels[i] > 0. else 1. for i in range(len(levels))]
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
                "levels":[],
                "corrected_levels":[],
                "levels_stat_errors":[],
                "levels_syst_errors":[],
                "detector_systematics":[],
                "performance_systematics":[],
            },
            "ds":{
                "levels":[],
                "corrected_levels":[],
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
        print "Loading", file_name
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
