import sys
import ROOT
import utilities

class ConglomerateMerge(object):
    """
    Merge individual plots of the same data into one uber-canvas

    Each plot is intended to have the same axes, so we can suppress axes to
    save white space
    """
    def __init__(self, conglomerate_list):
        """
        conglomerate_list: list of ConglomerateBase (each of which is a merge
        of MC and data that pertains to several plots for a particular beam 
        setting)
        """
        self.conglomerate_list = conglomerate_list
        self.merge_list = []
        self.rows = 0
        self.cols = 0
        self.cong_options = conglomerate_list[0].conglomerations[0]
        
    def merge_all(self, rows, cols):
        self.rows = rows
        self.cols = cols
        unique_id = set()
        for cong_base in self.conglomerate_list:
            for cong in cong_base.conglomerations:
                unique_id.add(cong.options["unique_id"])
        for name in unique_id:
            self.merge_one(name)

    def do_outer_box(self):
        if "outer_box" not in self.cong_options.options["merge_options"]:
            return
        options = self.cong_options.options["merge_options"]["outer_box"]
        hist_x = [0.1422, 0.3813, 0.6223, 0.8637]
        hist_y = [0.1345, 0.319, 0.5064, 0.6938, 0.8811]
        canvas_xmax = hist_x[-1]+0.115
        canvas_ymax = hist_y[-1]+0.05
        fill_color = 10003 #ROOT.kBlue-10 #ROOT.kTeal-9 #ROOT.kGray #
        transparency = 1.0 #0.2 # 0.3
        box = ROOT.TBox(hist_x[0], hist_y[-1], canvas_xmax, canvas_ymax)
        box.SetFillColor(fill_color) #'kBlue-10
        box.Draw()
        self.root_objects.append(box)
        box = ROOT.TBox(hist_x[-1], hist_y[0], canvas_xmax, hist_y[-1])
        box.SetFillColor(fill_color) #'kBlue-10
        box.Draw()
        self.root_objects.append(box)
        for i in range(self.cols+1):
            line = ROOT.TLine(hist_x[i], hist_y[-1], hist_x[i], canvas_ymax)
            line.Draw("l")
            self.root_objects.append(line)
        for i in range(self.rows+1):
            line = ROOT.TLine(hist_x[-1], hist_y[i], canvas_xmax, hist_y[i])
            line.Draw("l")
            self.root_objects.append(line)
        line = ROOT.TLine(canvas_xmax, hist_y[0], canvas_xmax, canvas_ymax)
        line.Draw("l")
        self.root_objects.append(line)
        line = ROOT.TLine(hist_x[0], canvas_ymax, canvas_xmax, canvas_ymax)
        line.Draw("l")
        self.root_objects.append(line)
        line = ROOT.TLine(hist_x[-1], hist_y[0], hist_x[-1], canvas_ymax)
        line.Draw("l")
        self.root_objects.append(line)
        line = ROOT.TLine(hist_x[0], hist_y[-1], canvas_xmax, hist_y[-1])
        line.Draw("l")
        self.root_objects.append(line)

        text_box = ROOT.TPaveText(hist_x[-1], hist_y[-1], canvas_xmax, canvas_ymax-0.01, "NDC")
        text_box.SetTextSize(0.012)
        text_box.SetFillStyle(0)
        text_box.SetBorderSize(0)
        text_box.SetTextAlign(11)
        text_box.AddText("ISIS User Runs")
        text_box.AddText("2017/02 and 2017/03")
        text_box.AddText("MAUS-v3.3.2")
        text_box.Draw()
        self.root_objects.append(text_box)



    def do_legend(self, i, j, legend, legend_size, pad):
        if legend == None:
            return
        # weird axis bounds!
        if j == 0:
            # x goes from 0.0 to 1.0/0.9?
            legend.SetX1(legend_size[0]*1.0/0.9)
            legend.SetX2(legend_size[2]*1.0/0.9)
        else:
            # x goes from -0.1 to 1.0/0.9?
            legend.SetX1((legend_size[0]-0.1)*1.0/0.9)
            legend.SetX2((legend_size[2])*1.0/0.9)
        if i == self.rows-1:
            # y goes from 0.0 to 1.0/0.9
            legend.SetY1(legend_size[1]*1.0/0.9)
            legend.SetY2(legend_size[3]*1.0/0.9)
        else:
            # y goes from -0.1 to 1.0/0.9
            legend.SetY1(legend_size[1]*1.0/0.9)
            legend.SetY2(legend_size[3]*1.0/0.9)
        legend.Draw()

    def do_legend_2(self):
        if "legend" not in self.cong_options.options["merge_options"]:
            return
        options = self.cong_options.options["merge_options"]["legend"]
        xmin = 0.75
        ymax = 0.13
        print "LEGEND", options
        for i, item in enumerate(options):
            yc = ymax-(i+1)*0.05
            box = ROOT.TBox(xmin, yc-0.003, xmin+0.05, yc+0.003)
            box.SetLineWidth(0)
            box.SetFillColorAlpha(item["fill_color"], item["transparency"])
            box.Draw()
            self.root_objects.append(box)
            if item["marker_style"]:
                marker = ROOT.TMarker(xmin+0.025, yc, item["marker_style"])
                marker.SetMarkerColor(item["marker_color"])
                marker.Draw()
                self.root_objects.append(marker)
            text_box = ROOT.TPaveText(xmin+0.05, yc-0.025, xmin+0.2, yc+0.025, "NDC")
            text_box.SetTextSize(0.04)
            text_box.SetFillStyle(0)
            text_box.SetBorderSize(0)
            text_box.SetTextAlign(12)
            text_box.AddText(item["text"])
            text_box.Draw()
            self.root_objects.append(text_box)

    def do_legend_3(self):
        if "legend" not in self.cong_options.options["merge_options"]:
            return
        options = self.cong_options.options["merge_options"]["legend"]
        xmin = 0.75
        ymax = 0.89
        print "LEGEND", options
        for i, item in enumerate(options):
            yc = ymax-(i+1)*0.026
            box = ROOT.TBox(xmin, yc-0.003, xmin+0.025, yc+0.003)
            box.SetLineWidth(0)
            box.SetFillColorAlpha(item["fill_color"], item["transparency"])
            box.Draw()
            self.root_objects.append(box)
            if item["marker_style"]:
                marker = ROOT.TMarker(xmin+0.0125, yc, item["marker_style"])
                marker.SetMarkerColor(item["marker_color"])
                marker.Draw()
                self.root_objects.append(marker)
            text_box = ROOT.TPaveText(xmin+0.03, yc-0.02, xmin+0.1, yc+0.02, "NDC")
            text_box.SetTextSize(0.02)
            text_box.SetFillStyle(0)
            text_box.SetBorderSize(0)
            text_box.SetTextAlign(12)
            text_box.AddText(item["text"])
            text_box.Draw()
            self.root_objects.append(text_box)


    def do_fit_tables(self):
        pass

    def do_mice_logo(self, do_isis, row, col):
        y_min=0.87
        x_min=0.02 #0.7
        if col == 0:
            x_min = 0.16
        if do_isis:
            y_min = 0.65
        text_box = ROOT.TPaveText(x_min, y_min, x_min+0.17, 0.97, "NDC")
        text_box.SetTextSize(0.08)
        text_box.SetFillColor(0)
        text_box.SetFillStyle(0)
        text_box.SetBorderSize(0)
        text_box.SetTextAlign(11)
        text_box.AddText("MICE")
        if do_isis:
            text_box.AddText("ISIS User Runs")
            text_box.AddText("2017/02 and 2017/03")
            text_box.AddText("MAUS-v3.3.2")
        text_box.Draw()
        self.root_objects.append(text_box)
        
    def draw_cong(self, pad, cong, i, j):
        pad.cd()
        legend_size = None
        if cong.legend != None:
            legend_size = cong.options["legend"]["pos"]
        frame_color = cong.pad.GetFrameFillColor()
        if "row_fill" in cong.options["merge_options"]:
            try:
                frame_color = cong.options["merge_options"]["row_fill"][i]
            except Exception:
                print "Failed to fill column", i, j, cong.options["merge_options"]["row_fill"]
        elif "col_fill" in cong.options["merge_options"]:
            try:
                frame_color = cong.options["merge_options"]["col_fill"][j]
            except Exception:
                print "Failed to fill column", i, j, cong.options["merge_options"]["col_fill"]
        print "Setting frame color", frame_color
        pad.SetFrameFillColor(frame_color)
        if j != 0:
            for hist in cong.hist_list:
                hist.GetYaxis().SetLabelOffset(100)
                #hist.GetYaxis().SetRangeUser(0, 100000)

        cong.label_size = 0.1
        cong.murgle_many(pad, cong.hist_list, cong.graph_list)
        if False: #i == 1 and j == 3:
            self.do_legend(i, j, cong.legend, legend_size, pad)

    def extra_labels(self, cong_list):
        extra_labels = self.cong_options.options["merge_options"]
        if not extra_labels:
            return
        left_m = 0.13
        right_m = 0.125
        top_m = 0.11
        bottom_m = 0.12
        text_height = 0.04
        x_step = (1.0-left_m-right_m)/self.cols
        for i, item in enumerate(extra_labels["top_labels"]):
            lines = item.split("\n")
            row_btm = 1.0-top_m
            row_top = row_btm+text_height*len(lines)
            text_box = ROOT.TPaveText(left_m+x_step*i, row_btm,
                                      left_m+x_step*(i+1), row_top, "NDC")
            text_box.SetTextSize(0.03)
            text_box.SetFillStyle(0)
            text_box.SetBorderSize(0)
            print "Setting text", item
            for line in lines:
                text_box.AddText(line)
            text_box.Draw()
            self.root_objects.append(text_box)
        y_step = (1.0-top_m-bottom_m)/self.rows
        for i, item in enumerate(extra_labels["right_labels"]):
            lines = item.split("\n")
            row_btm = 1.0-top_m-y_step*(i+0.5)-len(lines)/2.*text_height
            row_top = row_btm + len(lines)*text_height
            text_box = ROOT.TPaveText(1.0-right_m, row_top,
                                      0.99, row_btm, "NDC")
            text_box.SetTextSize(0.03)
            text_box.SetTextAlign(12)
            text_box.SetFillStyle(0)
            text_box.SetBorderSize(0)
            print "Setting text", item
            for line in lines:
                text_box.AddText(line)
            text_box.Draw()
            self.root_objects.append(text_box)

    def write_canvas(self, canvas, options):
        canvas.Update()
        target_dir = self.cong_options.config.target_plot_dir
        if options["write_plots"]["file_name"]:
            canvas_name = options["write_plots"]["file_name"]
            print "OPTIONS"
        else:
            canvas_name = options["canvas_name"]
        print "CANVAS NAME", canvas_name
        format_list = options["write_plots"]["formats"]
        canvas_name = canvas_name.replace("*", "")
        canvas_name = canvas_name.replace("?", "")
        for fmt in format_list:
            canvas.Print(target_dir+"/"+canvas_name+"."+fmt)
 
    def merge_one(self, canvas_name):
        """
        Merge one plot
        
        Draw data from all elements in self.conglomerate_list; find the "ith"
        plot and merge into one uber-canvas
        """
        print "Merging", canvas_name
        merge_dict = {}
        source_hist_list = []
        source_graph_list = []
        source_legend_list = []

        merge_canvas = ROOT.TCanvas(canvas_name+"_c1_", canvas_name+"_c1_", 1400, 1000)
        merge_dict["canvas"] = merge_canvas
        merge_canvas.Draw()
        hist_index = 0
        merge_canvas.Divide(1, 1, 0.1, 0.1)
        pad = merge_canvas.cd(1)
        pad.Divide(self.cols, self.rows, 0., 0.)
        for i in range(self.rows):
            for j in range(self.cols):
                try:
                    cong_base = self.conglomerate_list[hist_index]
                except Exception:
                    print "Exception while looping over pads in ", i, j
                    print "Failed to find hist", hist_index, "in", canvas_name
                for cong in cong_base.conglomerations:
                    if cong.options["unique_id"] == canvas_name:
                        break
                source_x_axis = cong.x_axis
                source_y_axis = cong.y_axis
                self.draw_cong(pad.GetPad(hist_index+1), cong, i, j)
                hist_index += 1
                if cong.options["mice_logo"]:
                    self.do_mice_logo(False, i, j) #i == 0 and j == self.cols-1)
        merge_canvas.cd()
        if source_x_axis:
            source_x_axis.Draw()
        if source_y_axis:
            source_y_axis.Draw()
        self.do_legend_3()
        self.do_outer_box()
        self.extra_labels(self.conglomerate_list)
        self.merge_list.append(merge_dict)
        self.write_canvas(merge_canvas, cong.options)
        
    root_objects = []
    light_fill_color = ROOT.TColor(10003, 0.94, 0.94, 1.0)

