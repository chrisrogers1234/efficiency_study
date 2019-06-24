import glob
import json
import copy

import numpy
import scipy.stats
from scipy.stats import chi2

class ExtractCoolingStats(object):
    def __init__(self):
        self.table = []
        self.rows = []
        self.headings = []
        self.job_name = ""
        self.all_data = {}

    def extract(self, source_glob):
        self.source_glob = source_glob
        self.load_json()

    @classmethod
    def get_job_name(cls, glob_name):
        job_name = glob_name.split("plots_")[1]
        job_name = job_name.split("/ampl")[0]
        job_name = job_name.replace("2017-2.7_", "")
        job_name = job_name.replace("_", " ")
        #job_name = job_name.ljust(26)
        return job_name

    def load_json(self):
        #print "Loading from", self.source_glob
        glob_list = sorted(glob.glob(self.source_glob))
        if len(glob_list) > 1:
            raise RuntimeError("Glob list too long")
        job_name = self.get_job_name(glob_list[0])
        print job_name
        self.data = json.loads(open(glob_list[0]).read())
        self.all_data[job_name] = copy.deepcopy(self.data)
        return job_name

    def print_one_dataset(self):
        for sample in "all_upstream", "all_downstream":
            print "   ", sample
            for key in "pdf",: #, "corrected_pdf", "pdf_sys_errors", "pdf_stats_errors":
                print "       ", key,
                for value in self.data['reco'][sample][key]:
                    print format(value),
                print
        print

    def stats_mangle(self, pdf_1, pdf_2):
        norm = sum(pdf_2)*1./sum(pdf_1)
        pdf_1_alt = [p for i, p in enumerate(pdf_1) if pdf_1[i] > 5 and pdf_2[i] > 5]
        pdf_2_alt = [p/norm for i, p in enumerate(pdf_2) if pdf_1[i] > 5 and pdf_2[i] > 5]
        return pdf_1_alt, pdf_2_alt

    def pdf_ds_chi2(self, us_ds, n_bins):
        chi2_table = []
        bl_list = ["4-140", "6-140", "10-140"]
        abs_list = ["LiH", "None", "lH2 full", "lH2 empty"]
        for i in range(4):
            for j in range(i+1, 4):
                abs_pair = (abs_list[i], abs_list[j])
                chi2_table.append([abs_pair[0]+" vs "+abs_pair[1]])
                for bl in bl_list:
                    conf_1 = bl+" "+abs_pair[0]
                    conf_2 = bl+" "+abs_pair[1]
                    pdf_1 = self.all_data[conf_1]["reco"][us_ds]["pdf"]
                    pdf_2 = self.all_data[conf_2]["reco"][us_ds]["pdf"]
                    pdf_1, pdf_2 = self.stats_mangle(pdf_1, pdf_2)
                    pdf_1 = numpy.array(pdf_1)
                    pdf_2 = numpy.array(pdf_2)
                    chisq, p = scipy.stats.chisquare(pdf_1[:n_bins], pdf_2[:n_bins])
                    chi2_table[-1].append([chisq, p])
                    #print "    chi square test "+conf_1+ " vs "+conf_2+" "+us_ds
                    #print "        Chi^2", chisq, "p", p
        us_ds = us_ds.replace("all_", "")
        print "Doing chi2 test for", us_ds, "with", n_bins, "bins"
        print "".ljust(25),
        for bl in bl_list:
            print "&", "\multicolumn{2}{c}{"+bl+"}",
        print "\\\\"
        print "".ljust(25),
        for i in range(3):
            print "& $\chi^2$ &  p-value",
        print "\\\\"
        for row in chi2_table:
            for element in row:
                if type(element) == type([]):
                    print "&", format(element[0], "8.2g"), "&", format(element[1], "8.2g"),
                else:
                    print element.ljust(25),
            print "\\\\"

    def pdf_us_ds_chi2(self, n_bins):
        bin_list = range(n_bins)
        chi2_table = []
        bl_bin_list = [("4-140", 1), ("6-140", 4), ("10-140", 6)]
        for an_abs in ["LiH", "None", "lH2 full", "lH2 empty"]:
            chi2_table.append([an_abs])
            for bl, cdf_bin in bl_bin_list:
                conf = bl+" "+an_abs
                pdf_1 = self.all_data[conf]["reco"]["all_upstream"]["corrected_pdf"]
                pdf_2 = self.all_data[conf]["reco"]["all_downstream"]["corrected_pdf"]
                err_1 = self.all_data[conf]["reco"]["all_upstream"]["pdf_sum_errors"]
                err_2 = self.all_data[conf]["reco"]["all_downstream"]["pdf_sum_errors"]
                self.chi2(pdf_1, err_1, pdf_2, err_2, bin_list)
                cdf_1 = self.all_data[conf]["reco"]["all_upstream"]["corrected_cdf"]
                cdf_2 = self.all_data[conf]["reco"]["all_downstream"]["corrected_cdf"]
                err_1 = self.all_data[conf]["reco"]["all_upstream"]["cdf_sum_errors"]
                err_2 = self.all_data[conf]["reco"]["all_downstream"]["cdf_sum_errors"]
                chi2_val, chi2_p = self.chi2(cdf_1, err_1, cdf_2, err_2, [cdf_bin])
                chi2_table[-1].append([chi2_val, chi2_p])
        print "Chi2 comparing upstream to downstream distributions"
        print "".ljust(10),
        for bl, cdf_bin in bl_bin_list:
            print "&", "\multicolumn{2}{c}{"+bl+"}",
        print "\\\\"
        print "".ljust(10),
        for i in range(3):
            print "& $\chi^2$ &  p-value",
        print "\\\\"
        for row in chi2_table:
            for element in row:
                if type(element) == type([]):
                    print "&", format(element[0], "8.2g"), "&", format(element[1], "8.2g"),
                else:
                    print element.ljust(10),
            print "\\\\"



    @classmethod
    def chi2(self, pdf_1, err_1, pdf_2, err_2, bin_list)  :
        # assume the pdf is already normalised as d/s is subset of u/s
        delta_list = []
        for i in bin_list:
            delta_list.append((pdf_1[i]-pdf_2[i])**2/((err_1[i])**2+err_2[i]**2))
        #print delta_list
        chi2_val = sum(delta_list)
        dof = len(bin_list)
        chi2_p = scipy.stats.chi2.cdf(chi2_val, dof)
        #print "    ", chi2_val, 1-chi2_p
        return chi2_val, 1-chi2_p

def main():
    prefix = "output/2017-02-7-v12/"
    stats = ExtractCoolingStats()
    stats.rows = ["None", "lH2_empty", "lH2_full", "LiH",]
    bin_headings = [("4-140"), ("6-140"), ("10-140")]
    stats.headings = [head[1] for head in bin_headings]
    table = []
    for absorber in stats.rows:
        for emittance in bin_headings:
            for real in ["_"]:
                suffix = "/amplitude/amplitude.json"
                a_glob = prefix+"plots"+real+"2017-2.7_"+emittance+"*"+absorber+"*"+suffix
                job_name = stats.extract(a_glob)
                stats.print_one_dataset()
    print "\n\n\n"
    stats.pdf_ds_chi2("all_upstream", 15)
    print
    stats.pdf_ds_chi2("all_downstream", 15)
    print "\n\n\n"
    stats.pdf_us_ds_chi2(6)

if __name__ == "__main__":
    main()