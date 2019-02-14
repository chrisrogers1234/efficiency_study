import tempfile
import json
import ROOT
import xboa.common

from mice_analysis.amplitude.amplitude_data_binned import AmplitudeDataBinned
from mice_analysis.fractional_emittance.amplitude_quantile import AmplitudeQuantile
from mice_analysis.analysis_base import AnalysisBase

class FractionalEmittance(AnalysisBase):

    def __init__(self, config, config_anal, data_loader):
        """
        Initialise the FractionalEmittance class for the fractional emittance analysis
        * config and config_anal are the configuration files
	* data_loader extracts the data from the MAUS ROOT files
        """
        super(FractionalEmittance, self).__init__(config, config_anal, data_loader)

	# Initialize the configuration and the data loader
        self.config = config
        self.config_anal = config_anal
        self.data_loader = data_loader

	# Initialize the fractional emittance parameters
        self.bin_edges = self.config.fractional_emittance_bins
        self.fraction = self.config.fractional_emittance_fraction
        self.uncertainty = self.config.fractional_emittance_fraction

	# Initialize the data containers
        self.a_dir = tempfile.mkdtemp()
        self.data_types = ("all_mc", "reco")
	self.data = {}
	for typ in self.data_types:
	    self.data[typ] = {}

        self.dict_list = []

	# Calculate corrections if required
        self.calculate_corrections = self.config_anal["density_corrections"] == None

    def birth(self):
        """
        Sets the output directory for the fractional emittance plots
	Resets all the data containers and loads the uncertainties
        """
	# Create directory, clear the data containers
        self.set_plot_dir("fractional_emittance")
	for typ in self.data_types:
	    self.data[typ].clear()

        self.get_planes()
        self.clear_feps_data()

	# Load the data
        self.append_data()

    def process(self):
        """
        Imports the data to the fractional emittance format
        """
        self.append_data()

    def death(self):
        """
        Called when all the data has been loaded
        """
	# First, evaluate the amplitude quantiles and their statistical uncertainties
        self.make_amplitude_quantiles()

	# Draw the fractional emittance evolution graphs
        self.draw_plot()

	# Save everything to a json file
        self.save()

    def get_planes(self):
        """
        Initialize a data container for each of the data type.
	* all_mc has a container for each virtual plane between the two ref. planes included
	* reco and reco has a container for each ref. plane
        """
	# Initialize z position container, muon mass
        self.z_pos = {}
        mass = xboa.common.pdg_pid_to_mass[13]

	# Initialize data container for each reco plane
        for name in ['tku', 'tkd']:
            for z, dummy, z_name in self.config.detectors:
                if z_name == name+"_tp":
                    self.z_pos[name] = z
                    break

            file_name = self.a_dir+"/amp_data_"+name+"_"
            self.data["reco"][name] = AmplitudeDataBinned(file_name, self.bin_edges, mass, False)
            self.data["reco"][name].clear()

        print "Set up reco planes"

	# If MC is present, initialize their data containers as well
	# Skip if the virtual plane is before tku or after tkd
        if self.config_anal["amplitude_mc"]:
            found_tku = False
            found_tkd = False
            for z, dummy, name in self.config.virtual_detectors:
                if "virtual_tku_tp" in name:
                    found_tku = True
                if not found_tku or found_tkd:
                    continue
                if "virtual_tkd_tp" in name:
                    found_tkd = True

                name = "mc_"+name
                file_name = self.a_dir+"/amp_data_"+name+"_"
                self.data["all_mc"][name] = AmplitudeDataBinned(file_name, self.bin_edges, mass, False)
                self.data["all_mc"][name].clear()
                self.z_pos[name] = z

            print "Set up", len(self.data["all_mc"]), "mc planes."


    def append_data(self):
        """
        Add data to the density calculation (done at death time)
        
        If density_mc is false, we take 'tku' data from upstream_cut sample 
        and 'tkd' data from downstream_cut sample
        
        If density_mc is true, then we build an additional samples
        1. all_mc is for MC truth of all events that should have been included 
           in the sample (some of which might have been missed by recon; some in
           recon sample maybe should have been excluded)
        We use "mc_true_us_cut"/"mc_true_ds_cut" for the all_mc sample
        """
	hits = {}
        for typ in self.data_types:
	    hits[typ] = {}
	    for loc in self.data[typ].keys():
		hits[typ][loc] = []

        for event in self.data_loader.events:
            if event['upstream_cut']:
                continue
            hits["reco"]["tku"].append(event['tku'])
            if not event['downstream_cut']:
                hits["reco"]["tkd"].append(event['tkd'])

            if self.config_anal["amplitude_mc"]:
                for detector_hit in event["data"]:
                    det = detector_hit["detector"]
                    if det in self.data["all_mc"].keys():
                        hits["all_mc"][det].append(detector_hit["hit"])

        for typ in self.data_types:
	    for loc in self.data[typ].keys():
		self.data[typ][loc].append_hits(hits[typ][loc])

    def make_amplitude_quantiles(self):
        """
        Calculates an amplitude array for each of the datum point.
	Feeds the amplitude arrays to the AmplitudeQuantile object
	which computes the quantiles and their stat. uncertainties
        """
	# Build the amplitude dictionaries
	self.build_amplitude_dictionaries()

	# Get the number of particles in the upstream tracker
        if self.dict_list[0][0] != 'tku':
            raise RuntimeError("First dict must be tku to do fractional emittance")
        us_events = len(self.dict_list[0][1])

	# Loop over the planes, get the quantiles
        for name, amps in self.dict_list:
	    # Initialize the quantile object
            amp_values = amps.values()
            n_events = len(amp_values)		
	    norm = float(n_events)/us_events
	    quantile = AmplitudeQuantile(4, norm, self.uncertainty)

	    # Get the quantile
            feps = -1.
	    feps_stat = 0.
	    if self.fraction < norm:
		quantile.set_quantile(amp_values, self.fraction)
                feps = quantile.value()
                feps_stat = quantile.error()

	    # Add a datum point for the current plane
            datum = {
              'name':name,
              'n_events':n_events,
              'fraction':self.fraction,
              'amplitude_quantile':feps,
              'stat_error':feps_stat,
              'z_pos':self.z_pos[name],
            }
            json.dumps(datum)
            self.feps_data.append(datum)

    def build_amplitude_dictionaries(self):
        """
        Temporarily build amplitude dictionaries
        """
        self.dict_list = [
            ('tku', self.data["reco"]["tku"].fractional_amplitude()),
            ('tkd', self.data["reco"]["tkd"].fractional_amplitude())
        ]
	for typ in self.data_types:
	    if typ == "reco":
		continue
            for name, amp_data in self.data[typ].iteritems():
                try:
                    self.dict_list.append((name, amp_data.fractional_amplitude()))
                except Exception:
                    print "Failed to get amplitudes for", name

    def draw_plot(self):
        """
        Produce a plot that shows the evolution of the fractional emittance
	across the MICE cooling channel. Displays the different types of data on
	shared canvas using the ROOT TMultiGraph class
        """
	# Initialize the canvas and the TMultiGraph
        canvas_name = 'fractional_emittance'
        canvas = self.get_plot(canvas_name)["pad"]
	mg = ROOT.TMultiGraph("mg_feps", ";z [m];#epsilon_{%d}  [mm]" % int(1e2*self.fraction))

	# Initialize the reco graph
        reco_predicate = lambda datum: datum['name'] == "tku" or datum['name'] == "tkd"
        name = "reco"
        point_list = [(datum['z_pos'], datum['amplitude_quantile'], datum['stat_error']) \
                                for datum in self.feps_data if reco_predicate(datum)]
        reco_graph = self.make_graph(point_list, 1, 20, name)
	mg.Add(reco_graph, "p")

	# Initialize the MC truth graph
        name = "all_mc"
        point_list = [(datum['z_pos'], datum['amplitude_quantile'], datum['stat_error']) \
                                for datum in self.feps_data if not reco_predicate(datum)]
        mc_graph = self.make_graph(point_list, 4, 24, name)
	mg.Add(mc_graph, "ple3")

	# Initialize a legend
	leg = ROOT.TLegend(.6, .65, .8, .85)
	leg.AddEntry(reco_graph, reco_graph.GetName(), "p")
	leg.AddEntry(mc_graph, mc_graph.GetName(), "plf")

	# Draw
	mg.Draw("A")
	leg.Draw("SAME")
        for fmt in ["pdf", "png", "root"]:
            canvas.Print(self.plot_dir+"/fractional_emittance."+fmt)

    def make_graph(self, point_list, color, style, name):
        """
        Produce a single graph associated with a list of points
	* point_list is a list of points to be respresented, array of (z, val, err)
	* color specifies the color of the graph and markers
	* style specifies the marker style
	* name specifies the name of the graph
        """
        n_points = len(point_list)
        graph = ROOT.TGraphErrors(n_points)
	graph.SetTitle(";z [m];#epsilon_{%d}  [mm]" % int(1e2*self.fraction))
        graph.SetLineColor(color)
        graph.SetFillColorAlpha(color, .25)
        graph.SetMarkerColor(color)
        graph.SetMarkerStyle(style)
        graph.SetName(name)

        point_list = sorted(point_list)
        for i, point in enumerate(point_list):
            z = point[0]/1e3
            val = point[1]
            err = point[2]
            graph.SetPoint(i, z, val)
            graph.SetPointError(i, 0, err)

        self.plots["fractional_emittance"]["graphs"][name] = graph

        return graph

    def save(self):
        fout = open(self.plot_dir+"/fractional_emittance.json", "w")
        print >> fout, json.dumps(self.feps_data, sort_keys=True, indent=4)

    def clear_feps_data(self):
        """
        Initializes the dictionary that is used to store
	data and make corrections down the line
        """
	self.feps_data = []
