import sys
import time


import numpy
import ROOT

import Configuration
import maus_cpp.globals
import maus_cpp.field as field # MAUS field map calls
from maus_cpp.global_error_tracking import GlobalErrorTracking
import maus_cpp.polynomial_map
import xboa.common
from xboa.bunch import Bunch
from xboa.hit import Hit


def initialise_maus():
    configuration = Configuration.Configuration().\
                                          getConfigJSON(command_line_args=True)
    maus_cpp.globals.birth(configuration)


class OpticsTool(object):
    def __init__(self):
        self.particles_in = []
        self.particles_out = []
        self.is_scraped = []
        self.tracker_z = 15068.0
        self.diffuser_z = 13750.
        self.matrix = None
        self.emittance_list = [20]#, 40, 80]
        self.n_events = 1000

    def setup_matrix(self):
        delta = [1., 0.1, 1., 0.1]
        events_in = [[0., 0., 0., 0.]]
        for i in range(4):
            an_event = [0., 0., 0., 0.]
            an_event[i] = delta[i]
            events_in.append(an_event)
        for event in events_in:
            hit = self.event_from_four_vector(event)
            self.particles_in.append(hit)
        self.track()
        events_out = [[hit[var] for var in self.phase_space] for hit in self.particles_out]
        self.matrix = maus_cpp.polynomial_map.PolynomialMap.least_squares_fit(
                                              events_in,
                                              events_out,
                                              1)

    def print_matrix(self):
        coeff = self.matrix.get_coefficients_as_matrix()
        for row in coeff:
            for cell in row:
                print format(cell, "12.4g"),
            print
        print

    def event_from_four_vector(self, event):
        mass = xboa.common.pdg_pid_to_mass[13]
        event_dict = {"x":event[0], "px":event[1], "y":event[2], "py":event[3],
                      "mass":mass, "z":self.tracker_z, "charge":1.}
        pz = (140**2-event[1]**2-event[3]**2)
        if pz > 0.:
            event_dict["pz"] = pz**0.5
        else:
            return None
        return Hit.new_from_dict(event_dict, "energy")

    def four_vector_from_event(self, hit):
        event = [hit[var] for var in self.phase_space]
        return event

    def setup_beam(self):
        self.particles_in = []
        for emittance in self.emittance_list:
            bz = field.get_field_value(0., 0., self.tracker_z, 0.)[2]
            if abs(bz) < 1e-3:
                bz = 3e-3
            p = 140.
            mass = xboa.common.pdg_pid_to_mass[13]
            beta = p/150./bz
            alpha = 0.
            q = 1.
            ellipse = Bunch.build_penn_ellipse(emittance, mass, beta,
                                                          alpha, p, 0., bz, q)
            #events = Bunch.new_hit_shell(7, ellipse, ["x", "px", "y", "py"], "", defaults)
            mean = numpy.array([0.]*4)
            events = numpy.random.multivariate_normal(mean, ellipse, self.n_events)
            for event in events:
                hit = self.event_from_four_vector(event)
                if hit != None:
                    self.particles_in.append(hit)
                else:
                    continue
        print "Generated", len(self.particles_in), "events"

    def track(self):
        print "Tracking..."
        self.particles_out = []
        var_list = ['t', 'x', 'y', 'z', 'energy', 'px', 'py', 'pz']
        tracking = GlobalErrorTracking()
        start_time = time.time()
        for hit_n, hit in enumerate(self.particles_in):
            print "\rProcessed ", hit_n, "hits in time", round(time.time() - start_time, 1), "s",
            sys.stdout.flush()
            centroid = [hit[var] for var in var_list]
            cov = numpy.identity(6).tolist()
            centroid, cov = tracking.propagate_errors(centroid, cov, self.diffuser_z)
            hit_dict = dict([(var, centroid[i]) for i, var in enumerate(var_list)])
            hit_dict['mass'] = xboa.common.pdg_pid_to_mass[13]
            hit = Hit.new_from_dict(hit_dict, "energy")
            self.particles_out.append(hit)
        print "\n...tracked"

    def matrix_transport(self):
        print "Doing matrix transport"
        self.particles_out = []
        for hit in self.particles_in:
            point_in = self.four_vector_from_event(hit)
            point_out = self.matrix.evaluate(point_in)
            self.particles_out.append(self.event_from_four_vector(point_out))

    def plot(self, prefix):
        self.set_scraping()
        if prefix != "":
            prefix = prefix+"_"
        #self.plot_one(self.particles_in, prefix+"tku_ref")
        #self.plot_one(self.particles_out, prefix+"diffuser")
        self.amplitude_weighted_scatter(True)
        self.amplitude_weighted_scatter(False)

    def set_scraping(self):
        self.is_scraped = []
        for hit in self.particles_out:
            self.is_scraped.append(hit['r'] > 90.)

    def calc_amplitude(self):
        events = [self.four_vector_from_event(hit) for hit in self.particles_in]
        events = numpy.array(events)
        var = numpy.cov(events.T)
        var_inv = numpy.linalg.inv(var)
        amplitude_list = [numpy.dot(event.T, numpy.dot(var_inv, event)) for event in events]
        return amplitude_list

    def plot_one(self, particles, plot_name):
        scraped_list = []
        not_scraped_list = []
        for i, hit in enumerate(particles):
            if self.is_scraped[i]:
                scraped_list.append(hit)
            else:
                not_scraped_list.append(hit)
        all_bunch = Bunch.new_from_hits(particles)
        scraped_bunch = Bunch.new_from_hits(scraped_list)
        not_scraped_bunch = Bunch.new_from_hits(not_scraped_list)
        for axes in [('y', 'px'), ('x', 'y')]:
            canvas, hist, graph = all_bunch.root_scatter_graph(axes[0], axes[1], xmin=-200., xmax=200., ymin=-200, ymax=200)
            if len(not_scraped_bunch) > 0:
                canvas, hist, graph = not_scraped_bunch.root_scatter_graph\
                                              (axes[0], axes[1], canvas=canvas)
                graph.SetMarkerStyle(7)
            if len(scraped_bunch) > 0:
                canvas, hist, graph = scraped_bunch.root_scatter_graph\
                                              (axes[0], axes[1], canvas=canvas)
                graph.SetMarkerColor(ROOT.kRed)
                graph.SetMarkerStyle(7)
            canvas.Update()
            canvas.Print(plot_name+"_"+axes[0]+"_vs_"+axes[1]+".png")
            canvas, hist = not_scraped_bunch.root_histogram(axes[0], "", axes[1], "", xmin=-200., xmax=200., ymin=-200, ymax=200)
            canvas.Print(plot_name+"_hist_"+axes[0]+"_vs_"+axes[1]+".png")


    def amplitude_weighted_scatter(self, will_cut):
        amplitude_list = self.calc_amplitude()
        amplitude_list = [(amp, i) for i, amp in enumerate(amplitude_list)]
        amplitude_list = reversed(sorted(amplitude_list))
        canvas = xboa.common.make_root_canvas("coloured scatter")
        graph2d = ROOT.TGraph2D(len(self.particles_in))
        index = 0
        for amp, i in amplitude_list:
            hit = self.particles_in[i]
            if will_cut and self.is_scraped[i]:
                continue
            graph2d.SetPoint(index, hit['y'], hit['px'], amp)
            index += 1
        ROOT.gStyle.SetPalette(1)
        graph2d.SetMarkerStyle(20)
        graph2d.Draw('pcol')
        self.root_objects.append(graph2d)
        canvas.SetTheta(-90)
        canvas.SetPhi(90)
        if will_cut:
            will_cut_str = "_scraped"
        else:
            will_cut_str = "_not_scraped"
        canvas.Print("tku_ref_coloured_y_vs_px"+will_cut_str+".png")

    root_objects = []
    phase_space = ['x', 'px', 'y', 'py']

def set_defaults():
    optics = OpticsTool()
    optics.n_events = 10000
    return optics


def main():
    initialise_maus()
    optics = set_defaults()
    optics.setup_matrix()
    print "Setup transfer matrix"
    optics.print_matrix()
    optics.setup_beam()
    optics.matrix_transport()
    optics.plot("matrix")
    #optics.track()
    #optics.plot("tracking")

if __name__ == "__main__":
    main()