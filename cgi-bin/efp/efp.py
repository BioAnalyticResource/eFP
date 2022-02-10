#!/usr/bin/python3
import glob
import os.path
import tempfile
import MySQLdb
import PIL.Image
import PIL.ImageDraw
import PIL.ImageFont
import lxml.sax
import math
import re
from lxml import etree
from xml.sax.handler import ContentHandler
from . import efpBase
from . import efpDb

# set HOME environment variable to a directory the httpd server can write to (for matplotlib)
# Todo: this line could be deleted
os.environ['HOME'] = '/tmp/'
import matplotlib

# Todo: this line could be deleted
matplotlib.use('Agg')  # set image engine, needed for png creation
import pylab
import numpy as np


class Sample:
    def __init__(self, name, view):
        self.name = name
        self.view = view
        self.signals = {}

    def get_signal(self, gene_id, probeset_id):
        """
        Check to see if a gene_id has a signal.
        If not, use it's probeset it to get a signal
        @param gene_id:
        @param probeset_id
        @return:
        """
        if gene_id not in self.signals:
            self.signals[gene_id] = self.view.get_tissue_signal(probeset_id, self.name)

        return self.signals[gene_id]


class Tissue:
    def __init__(self, name, colorKey):
        self.name = name
        if colorKey == '':
            colorKey = '#EEEEEE'
        self.color_string = colorKey
        self.colorKey = efpBase.html_color_to_rgb(colorKey)
        self.samples = []  # List of project set_id strings
        self.url = ''
        self.coords = []
        self.control = None

    def add_url(self, url):
        self.url = url

    def add_coords(self, coords):
        self.coords.append(coords)

    def add_sample(self, sample):
        self.samples.append(sample)

    def set_control(self, ctrl):
        self.control = ctrl

    def get_mean_signal(self, gene):
        """
        For the given gene, returns the average signal among all samples
        stored in the tissue as well as the standard deviation
        @param gene:
        @return: The tissue's average signal strength and stddev
        """
        if len(self.samples) == 0:
            return 0.0

        values = []
        mean = 0.0
        i = 0
        num_samples = len(self.samples)

        while i < num_samples:
            signal = self.samples[i].get_signal(gene.gene_id, gene.probeset_id)
            if signal is None:
                return None, None
            else:
                values.append(signal)
                mean += signal
                i += 1

        mean /= num_samples
        stddev = math.sqrt(sum([(x - mean) ** 2 for x in values]) / num_samples)
        stddev = round(stddev, 2)
        return mean, stddev


class Group:
    def __init__(self, name):
        self.tissues = []
        self.name = name
        self.ctrl_samples = []

    def add_ctrl_sample(self, control_sample):
        self.ctrl_samples.append(control_sample)

    def add_tissue(self, tissue):
        self.tissues.append(tissue)

    def get_control_signal(self, gene):
        if len(self.ctrl_samples) == 0:
            return 0.0

        mean = 0.0
        i = 0
        num_samples = len(self.ctrl_samples)

        while i < num_samples:
            if self.ctrl_samples[i].get_signal(gene.gene_id, gene.probeset_id) is None:
                return None
            mean += self.ctrl_samples[i].get_signal(gene.gene_id, gene.probeset_id)
            i += 1

        mean /= num_samples
        mean = round(mean, 3)

        return mean


class Extra:
    def __init__(self, name, link, parameters, coords, check, check_column):
        self.name = name
        self.link = link
        self.button = False

        if check:
            self.check = check
        else:
            self.check = ""

        if check_column:
            self.checkColumn = int(check_column)

        if parameters == "Yes":
            self.parameters = True
        else:
            self.parameters = False

        self.coords = coords


class View:
    chart_space_per_char = 0.06

    def __init__(self, name, db, dbGroup, image, conf):
        self.groups = []
        self.name = name
        self.database = db
        self.colorMap = None
        if os.path.exists(image):
            self.colorMap = PIL.Image.open(image)  # Original color map
        self.image = image
        self.extras = []
        self.graph = (0, 0, 0, 0)
        self.legend = (15, 30)
        self.legendSize = 12
        self.table = ''
        self.signals = []
        if dbGroup is None:
            dbGroup = conf['dbGroupDefault']
        self.dbGroup = dbGroup
        self.conn = None
        self.conf = conf

    @staticmethod
    def _clamp(my_input, my_min, my_max):
        """
        Returns my_input clamped to [my_min, my_max]

        @param my_input Input value
        @param my_min Minimum value
        @param my_max Maximum value
        @return:
        """
        if my_input > my_max:
            return my_max
        if my_input < my_min:
            return my_min
        return my_input

    @staticmethod
    def _get_gene_list(file_name, n):
        """
        get_gene_list returns a list of genes contained in a tab delimited file, commented lines start with #
        begin at the second line of the file

        @param file_name the file name
        @param n the column that contains the gene agis
        """
        gene_list = []

        # Check if file name is empty
        if file_name is None:
            return gene_list

        try:
            gene_file = open(file_name, 'r')
            column = int(n)
            lines = gene_file.readlines()

            # appends each gene agi to the gene list if the line contains a gene name
            for i in range(1, len(lines)):
                lines[i] = lines[i].strip()
                s = lines[i].split('\t')

                if s[column] and (not s[column][0] == '#'):
                    gene_list.append(s[column])

            gene_file.close()
        except FileNotFoundError:
            print('Content-Type: text/html\n')
            print("A data (txt) file is not found. Please contact us.")
            exit()

        return gene_list

    def add_extra(self, extra):
        self.extras.append(extra)

    def add_group(self, group):
        self.groups.append(group)

    def add_graph_coords(self, graph):
        self.graph = graph

    def add_legend_coords(self, legend):
        self.legend = legend

    def get_database(self):
        return self.database

    def get_image_path(self):
        return self.image

    def get_view_max_signal(self, mode, gene, gene2=None):
        """
        This is used to calculate max/min in legend and color scales.
        modified Nov 4 2009 (HN): - Some sample signals in the Root DB are equal to 0 (no signal), and cause errors
        during calculations. Added condition to check if signal is zero. If signal/control != 0,
        signal = abs(math.log(signal/control)/math.log(2)), else signal = 0
        @param mode:
        @param gene:
        @param gene2:
        @return:
        """
        view_max_signal = 0.0   # This is the view max signal
        view_max_gene = 0.0     # This is for the little red line in graph on some eFPs
        view_max_gene2 = 0.0    # This is for little blue line in graph on some eFPs
        tissue_signal = 0.0
        log2 = math.log(2)

        for group in self.groups:
            control = group.get_control_signal(gene)

            if mode == "Compare":
                control2 = group.get_control_signal(gene2)
            else:
                control2 = None

            for tissue in group.tissues:
                # Common to all
                (signal, stddev) = tissue.get_mean_signal(gene)

                if mode == "Absolute":
                    # No change in Absolute mode
                    if signal is None:
                        tissue_signal = 0
                    else:
                        # Set gene signal
                        if signal > view_max_gene:
                            view_max_gene = signal

                        tissue_signal = signal

                elif mode == "Relative":
                    if (signal is None) or (control is None) or (signal == 0) or (control == 0):
                        tissue_signal = 0
                    else:
                        # Set gene signal
                        if signal > view_max_gene:
                            view_max_gene = signal

                        # Get log ratio
                        signal = math.log(signal / control) / log2

                        # Now it's +ve number (we need max abs)
                        tissue_signal = abs(signal)

                elif mode == "Compare":

                    (signal2, stddev2) = tissue.get_mean_signal(gene2)

                    if (signal is None) or (signal2 is None) or (control is None) or (control2 is None) \
                            or (signal == 0) or (signal2 == 0) or (control == 0) or (control2 == 0):
                        tissue_signal = 0
                    else:
                        # Set gene signal
                        if signal > view_max_gene:
                            view_max_gene = tissue_signal

                        # Set gene2 signal
                        if signal2 > view_max_gene2:
                            view_max_gene2 = signal2

                        # Get log ratio
                        signal = math.log((signal / control) / (signal2 / control2)) / log2
                        tissue_signal = abs(signal)
                else:
                    # Not sure how we would get here
                    pass

                # Here tissue signal must be positive
                if tissue_signal is not None:
                    if tissue_signal > view_max_signal:
                        view_max_signal = tissue_signal

        view_max_signal = round(view_max_signal, 2)
        return view_max_signal, view_max_gene, view_max_gene2

    def start_table(self, mode="Absolute"):
        """
        Set up for the Table of Expression Values
        @param mode:
        @return:
        """
        self.table += '<style type=text/css>\n'
        # Background color of the Rows Alternates
        self.table += 'tr.r0 {background-color:#FFFFDD}\n'
        self.table += 'tr.r1 {background-color:#FFFFFF}\n'
        self.table += 'tr.rt {background-color:#FFFF99}\n'
        self.table += 'tr.rg {background-color:#DDDDDD}\n'
        self.table += 'td {font-family:arial;font-size:8pt;}\n'
        self.table += '</style>\n'
        self.table += '<table cellspacing=0 border=1 cellpadding=2 align=center>\n'

        # Column Headings

        self.table += '<tr class=rt><td><B>Group #</B></td><td><B>Tissue</B></td>'

        if mode == "Relative":
            self.table += '<td><B>Sample signal</B></td><td><B>Control signal</B></td>'
            self.table += '<td><B>Log2 Ratio</B></td><td><B>Fold-Change</B></td>'
        elif mode == "Compare":
            self.table += '<td><B>Log2 Ratio</B></td><td><B>Fold-Change</B></td>'
        else:
            self.table += '<td><B>Expression Level</B></td><td><B>Standard Deviation</B></td>'

        self.table += '<td><B>Samples</B></td><td><B>Links</B></td></tr>\n'

    def append_table(self, mode, tissue, value, n, ratio, stddev, sample_sig, control_sig, color):
        """
        Produces a row in the Table for each tissue
        @param mode:
        @param tissue:
        @param value:
        @param n:
        @param ratio: boolean
        @param stddev:
        @param sample_sig:
        @param control_sig:
        @param color:
        @return:
        """
        signal_dict = {'group': n}

        # Add Group number and tissue name to the table
        self.table += '<tr class=r%s><td>%s</td><td>%s</td>' % ((n % 2), n, tissue.name)

        if value is None:
            val_floor = np.nan
        else:
            val_floor = round(value, 2)

        if control_sig is None:
            control_sig_floor = np.nan
        else:
            control_sig_floor = round(control_sig, 2)

        if sample_sig is None:
            signal_dict['sample_sig'] = val_floor
            sample_sig_floor = np.nan
        else:
            sample_sig_floor = round(sample_sig, 2)
            signal_dict['sample_sig'] = sample_sig_floor - control_sig_floor

        # Fold Change for Relative and Compare Modes
        if ratio:
            if value is None:
                fold = np.nan
            else:
                fold = round(math.pow(2, value), 2)

            signal_dict['ratio'] = fold

        else:
            if stddev is None:
                stddev = np.nan
            signal_dict['stddev'] = stddev

        # Append to table
        if mode == "Absolute":
            self.table += '<td align=right>%s</td>' % val_floor
            self.table += '<td align=right>%s</td>' % stddev
        elif mode == "Relative":
            self.table += '<td align=right>%s</td><td align=right>%s</td><td align=right>%s</td>' % (
                sample_sig_floor, control_sig_floor, val_floor)
            self.table += '<td align=right>%s</td>' % fold
        else:
            # Compare mode
            # The val_floor here is actually Log2Ratio
            self.table += '<td align=right>%s</td>' % val_floor
            self.table += '<td align=right>%s</td>' % fold

        self.table += '<td>'

        # Add sample name to the table
        first = True
        for sample in tissue.samples:
            if first:
                self.table += '%s' % sample.name
                first = False
            else:
                self.table += ', %s' % sample.name

        # Add signal and link to table
        self.table += '</td><td><A target="_blank" href=%s>To the Experiment</A></td></tr>\n' % tissue.url

        signal_dict['signal_sample_name'] = tissue.name
        signal_dict['sample_signal_color'] = color

        self.signals.append(signal_dict)

    def end_table(self):
        self.table += '</table>\n'

    def save_chart(self, filename, mode):
        """
        This functions prints the charts in eFP Browser
        @param filename:
        @param mode:
        @return:
        """
        data_count = len(self.signals)
        x = np.arange(data_count)
        y = []
        colors = []
        ratio = []
        color_arr = ('#FFFFFF', '#FFFFDD')
        bar_width = 0.9
        box_color = '#777777'
        grid_color = '#999999'
        samples = []
        stddev = []
        groups = {}
        bg_colors = []
        max_sample = 0

        for signal in self.signals:
            y.append(signal['sample_sig'])
            # collect sample names
            # insert ' ' upfront to get background color down to bottom
            # insert 'Iy' upfront to ensure equal height of text line for all samples
            # so with background color there are now white lines in between
            samples.append('Iy' + ' ' * 500 + signal['signal_sample_name'])

            if len(signal['signal_sample_name']) > max_sample:
                max_sample = len(signal['signal_sample_name'])

            colors.append(color_arr[signal['group'] % 2])

            if 'stddev' in signal:
                stddev.append(signal['stddev'])
            if 'ratio' in signal:
                ratio.append(signal['ratio'])
            if signal['group'] in groups:
                groups[signal['group']] = groups[signal['group']] + 1
            else:
                groups[signal['group']] = 1
            bg_colors.append(signal['sample_signal_color'])

        # initialize image with size depending on amount of values and length of sample names
        plot = pylab.figure(frameon=False, dpi=180)  # initialize image
        img_width = 1.6 + data_count * 0.17
        left_border = 0.8 / img_width
        img_height = 3 + max_sample * self.chart_space_per_char
        bottom_border = max_sample * self.chart_space_per_char / img_height
        plot.set_size_inches(img_width, img_height)
        pylab.subplots_adjust(bottom=bottom_border, left=left_border, right=(1 - left_border), top=0.95, wspace=0,
                              hspace=0)  # make room for long x-axis labels (tissue names)
        ax1 = pylab.subplot(111)

        # plot colored background for individual groups
        n = 0
        c = 1
        for group in list(groups.values()):
            ax1.axvspan(n, n + group, facecolor=color_arr[c % 2], linewidth=0.0, label=('Group %d' % c))
            c = c + 1
            n = n + group

        # plot chart depending on mode
        if mode == "Absolute":
            ax1.bar(x + bar_width / 2., y, bar_width, color=box_color, linewidth=0, yerr=stddev)
            ax1.yaxis.grid(True, linewidth=0.5, color=grid_color)
            ax1.set_ylabel("%s" % (self.conf['Y_AXIS'][self.database]), size='x-small')

        elif mode == "Relative":
            ax1.bar(x + bar_width / 2., y, bar_width, color=box_color, linewidth=0)
            ax1.yaxis.grid(True, linewidth=0.5, color=grid_color)
            ax1.set_ylabel("difference in %s" % (self.conf['Y_AXIS'][self.database]), size='x-small')

            ax2 = pylab.twinx()
            ax2.plot(x + bar_width / 2., ratio, color='b', linestyle='None', marker='.', label=None)  # , visible=False
            ax2.set_ylabel('fold change', size='x-small')

            # workaround for printing labels on x-axis twice (fixed in more current version of matplotlib
            for tl in ax2.get_xticklabels():
                tl.set_visible(False)
            for t_line in ax2.get_xticklines():
                t_line.set_visible(False)
        else:
            ax1.plot(x + bar_width / 2., ratio, color='b', linestyle='None', marker='.', label=None)
            ax1.set_ylabel('fold change', size='x-small')

        # label x axis
        ax1.set_xticks(x + bar_width / 1.2)
        for t_line in ax1.get_xticklines():
            t_line.set_visible(False)
        labels = ax1.set_xticklabels(samples, rotation=90, ha='right', size='xx-small')

        # set background color according to signal as in efp original image
        for i in range(len(labels)):
            labels[i].set_backgroundcolor(efpBase.rgb_to_html_color(bg_colors[i]))
            if efpBase.rgb_to_gray(bg_colors[i]) < 186:
                labels[i].set_color('#FFFFFF')  # set foreground color to white if background is dark

        ax1.set_xlim(left=0, right=data_count)
        pylab.savefig(filename, dpi=180)

    def get_image_map(self, mode, gene1, gene2, useT, threshold, datasource, grey_low, grey_stddev):
        """
        Forms a Map Covering the Image of Hyperlinks and Drop-Down Expression Values modified Nov 3 2009 (HN): -
        Some sample signals in the Root DB are equal to 0 (no signal), and cause errors during calculations.
        Added condition to check if signal is zero. If Absolute mode, if signal = 0, Level = 0; else
        level = math.floor(signal * 100). Else if mode is Relative/Compute: if signal != 0:
        signal = math.log(signal / control) / math.log(2)
        """
        out = '<map name="imgmap_%s">' % self.name
        log2 = math.log(2)
        # This is for eFP Human and eFP Tomato
        zero_gene = False

        for extra in self.extras:
            if extra.parameters:
                if useT is None:
                    threshold_switch = ""
                    out += '<area shape="polygon" coords="%s" title="%s" href="%s&mode=%s&primaryGene=%s&secondaryGene=%s&useThreshold=%s&threshold=%s&grey_low=%s&grey_stddev=%s">\n' % (
                        extra.coords, extra.name, extra.link, mode, gene1, gene2, threshold_switch, threshold, grey_low,
                        grey_stddev)
                else:
                    out += '<area shape="polygon" coords="%s" title="%s" href="%s&mode=%s&primaryGene=%s&secondaryGene=%s&useThreshold=%s&threshold=%s&grey_low=%s&grey_stddev=%s">\n' % (
                        extra.coords, extra.name, extra.link, mode, gene1, gene2, useT, threshold, grey_low,
                        grey_stddev)
            else:
                # not a heat map button
                if extra.check == "":
                    out += '<area shape="polygon" coords="%s" title="%s" href="%s">\n' % (
                        extra.coords, extra.name, extra.link)

                    # is a heat map button
                else:
                    gene_list = self._get_gene_list(extra.check, extra.checkColumn)
                    # if the searched gene is contained in the heatmap list, activate the link
                    # if extra.button = 1, then gene1 is contained in the list
                    # if extra.button = 2, then gene2 is contained in the list
                    # if extra.button = 3, then both genes are contained in the list
                    if gene1.get_gene_id() in gene_list:
                        extra.button = 1
                        if self.conf['species'] in ["TOMATO", "HUMAN"]:
                            zero_gene = True
                    if gene2.get_gene_id() in gene_list:
                        extra.button = 2
                    if (gene1.get_gene_id() in gene_list) and (gene2.get_gene_id() in gene_list):
                        extra.button = 3

                    # draw button img map
                if extra.button == 1:
                    out += '<area shape="polygon" coords="%s" title="%s" href="%s%s">\n' % (
                        extra.coords, extra.name, extra.link, gene1.get_gene_id())
                elif extra.button == 2:
                    out += '<area shape="polygon" coords="%s" title="%s" href="%s%s">\n' % (
                        extra.coords, extra.name, extra.link, gene2.get_gene_id())
                elif extra.button == 3:
                    out += '<area shape="polygon" coords="%s" title="%s" href="%s%s%%0A%s">\n' % (
                        extra.coords, extra.name, extra.link, gene1.get_gene_id(), gene2.get_gene_id())
                else:
                    out += '<area shape="polygon" coords="%s" title="%s" href="%s">\n' % (
                        extra.coords, extra.name, extra.link)

        for group in self.groups:
            # Get control for the group here:
            if mode == "Relative":
                control = group.get_control_signal(gene1)
            elif mode == "Compare":
                control = group.get_control_signal(gene1)
                control2 = group.get_control_signal(gene2)
            else:
                # For Absolute and bugged eFP
                control = None
                control2 = None

            for tissue in group.tissues:
                (signal, stddev) = tissue.get_mean_signal(gene1)

                if signal is None:
                    # Note: fold is not needed here
                    sig_string = "No data available"

                else:
                    # Signal is not None for here onwards
                    if mode == "Absolute":
                        sig_floor = round(signal, 2)
                        sig_string = "Level: %s, SD: %s" % (sig_floor, stddev)

                    elif mode == "Relative":
                        # We already have controls
                        if (control is None) or (signal == 0) or (control == 0):
                            sig_string = "No data available"
                        else:
                            # The actual math and science of Relative mode
                            signal = math.log(signal / control) / log2
                            sig_floor = round(signal, 2)
                            fold = round(math.pow(2, signal), 2)
                            sig_string = "Log2 Value: %s, Fold-Change: %s" % (sig_floor, fold)

                    elif mode == "Compare":
                        # Get Signals, We already have controls
                        (value2, stddev2) = tissue.get_mean_signal(gene2)

                        if (value2 is None) or (control is None) or (control2 is None) or \
                                (signal == 0) or (value2 == 0) or (control == 0) or (control2 == 0):
                            sig_string = "No data available"
                        else:
                            # No nones or zeros
                            signal = math.log(signal / control) / log2
                            signal2 = math.log(value2 / control2) / log2
                            signal = signal - signal2
                            sig_floor = round(signal, 2)
                            fold = round(math.pow(2, signal), 2)
                            sig_string = "Log2 Value: %s, Fold-Change: %s" % (sig_floor, fold)

                    else:
                        # So Mode is none of the 3. How did we get here?
                        sig_string = "An error has occurred. Please report to the BAR."

                for coords in tissue.coords:
                    out += '<area shape="polygon" coords="%s" title="%s \n%s" href="%s">\n' % (
                        coords, tissue.name, sig_string, tissue.url)

        # Determining coordinates of triangle for link from graph
        x = int(self.graph[0])
        y = int(self.graph[1])
        w = int(self.graph[2])
        h = int(self.graph[3])
        coords = '%i,%i,%i,%i,%i,%i' % (x, y - h, x + w, y, x, y)
        out += '<area shape="polygon" coords="%s" title="The red line indicates the maximum expression of your primary gene, while the blue line, if present, indicates the maximum expression of your secondary gene" href="http://bbc.botany.utoronto.ca/affydb/BAR_instructions.html#efp_distro_%s">\n' % (
            coords, self.database)
        out += '</map>'
        return out, zero_gene

    def render_legend(self, img, title, efp_max, efp_min, stages=11, less_than=True, greater_than=True,
                      is_relative=False):
        """
        Draws a color legend

        @param img: Target image
        @param title: Legend title, drawn at the top of the legend
        @param efp_max:
        @param efp_min: Minimum value
        @param stages:
        @param less_than:
        @param greater_than:
        @param is_relative:
        @return:
        """
        # if efp_max == efp_min draw only on gradient stage in legend
        if efp_max == efp_min:
            stages = 1
            signal_grad = 0
        else:
            # Change in signal per step, always descending from efp_max
            signal_delta = efp_min - efp_max
            signal_grad = abs(signal_delta) / (stages - 1)

        # Start out with maximum signal
        signal = efp_max

        # Height of each legend item, in pixels (hardcoded for now ...)
        height = self.legendSize

        draw = PIL.ImageDraw.Draw(img)

        # Load a PIL font - these fonts are available from http://effbot.org/downloads/
        font_sizes = [(12, '08'), (15, '10'), (18, '12'), (21, '14'), (27, '18'), (36, '24')]
        fontsize = '08'
        for h, s in font_sizes:
            if height >= h:
                fontsize = s
        font = PIL.ImageFont.load(self.conf['STATIC_FILES'] + ("/pilfonts/helvR%s.pil" % fontsize))

        # Get top left coordinates of legend
        left = self.legend[0]
        # bottom = img.size[1] - (stages + 1) * height - self.legend[1]
        bottom = self.legend[1]

        draw.text((left, bottom), title, font=font, fill=(0, 0, 0))
        bottom += height + 2

        for y in range(stages):
            signal_output = " "

            # If we're at either the top or the bottom row of the legend, edit the
            # output string to include less than or greater than signs
            if y == 0 and greater_than:
                signal_output += "> "
            elif y == stages - 1 and less_than:
                signal_output += "< "

            # Ensure True Minimum
            if y == stages - 1:
                signal = efp_min

            # Keep two decimal points accuracy
            signal_output += str(round(signal, 2))

            if efp_max == 0:
                intensity = 0
            else:
                intensity = int(signal * 255.0 / efp_max)

                # Clamp intensity to [-255, 255]
                intensity = self._clamp(intensity, -255, 255)

            # Draw the colored rectangle
            if signal > 0:
                # Yellow (neutral) to red
                color = (255, 255 - intensity, 0)
            else:
                # Yellow to blue
                color = (255 + intensity, 255 + intensity, - intensity)

            # Draw the colored box
            draw.rectangle((left, bottom + y * height, left + 12, bottom + (y + 1) * height), fill=color)

            # Explanation of Relative Scale
            if y == 0 and is_relative:
                fold = math.pow(2, signal)
                fold_dec = (fold % 1) * 10
                signal_output += "  (%i.%i-Fold)" % (fold, fold_dec)

            # Draw the signal value
            draw.text((left + 12, bottom + (y * height)), signal_output, font=font, fill=(0, 0, 0))

            # Go from max to min
            signal -= signal_grad

        draw.rectangle((left, bottom + ((y + 1) * height), left + 12, bottom + (y + 2) * height), fill=(204, 204, 204))
        draw.text((left + 12, bottom + ((y + 1) * height)), " Masked", font=font, fill=(0, 0, 0))

    def render_absolute(self, gene, threshold=0.0, grey_mask=False):
        """
        Renders tissue data on a scale of the maximum signal.
        :param gene:
        :param threshold:
        :param grey_mask:
        :return: A PIL Image object containing the final rendered data
        """
        out_image = self.colorMap.copy()
        max_signal, max_signal_gene, max_signal_gene2 = self.get_view_max_signal("Absolute", gene)
        max_greater = False

        if threshold >= self.conf['minThreshold_Absolute']:
            efp_max = threshold
            if max_signal > threshold:
                max_greater = True
        else:
            # If the user doesn't give us a reasonable value for threshold,
            # use the maximum signal from dbGroup for this gene
            efp_max = max_signal

        n = 1
        sd_alert = 0
        self.start_table("Absolute")

        for group in self.groups:
            for tissue in group.tissues:
                # If for some reason this tissue object doesn't have a color key
                # assigned (malformed XML data?), skip it
                if tissue.colorKey == (0, 0, 0):
                    continue

                (signal, stddev) = tissue.get_mean_signal(gene)
                if signal is None:
                    color = (221, 221, 221)
                    sd_alert = 1
                    stddev = None
                else:
                    if signal == 0:
                        color = (255, 255, 0)  # Yellow
                        stddev = 0  # If signal is zero, set stddev to zero as well

                    else:
                        # In case if Max value is zero
                        if efp_max == 0:
                            intensity = 0
                        else:
                            intensity = int(math.floor(signal * 255.0 / efp_max))

                        # Grey out expression levels with high standard deviations
                        if (stddev / signal) > 0.5:
                            if grey_mask == 'on':
                                color = (221, 221, 221)  # Masked: CCCCCC
                            else:
                                sd_alert = 1
                                color = (255, 255 - intensity, 0)
                        else:
                            # Otherwise, color appropriately
                            color = (255, 255 - intensity, 0)

                # Append to table
                self.append_table("Absolute", tissue, signal, n, False, stddev, None, None, tissue.colorKey)

                # Perform fast color replacement
                out_image.replaceFill(self.colorMap, tissue.colorKey, color)

            n += 1

        # Complete Table of Expression Values
        self.end_table()
        self.render_legend(out_image, "Absolute", efp_max, 0, less_than=False, greater_than=max_greater)
        return out_image, max_signal, max_signal_gene, max_signal_gene2, sd_alert

    def render_relative(self, gene, threshold=0.0, grey_mask=False):
        """
        Renders tissue data relative to the control signal on a scale of the
        maximum signal. Note that unlike the 'efp_max' signal, there is more than
        one control compared against. Groups of tissues share the same control
        signal against which we compare.
        :param gene:
        :param threshold:
        :param grey_mask:
        :return: A PIL Image object containing the final rendered data
        """

        out_image = self.colorMap.copy()
        max_signal, max_signal_gene, max_signal_gene2 = self.get_view_max_signal("Relative", gene)
        max_greater = False

        if threshold >= self.conf['minThreshold_Relative']:
            max_log2 = threshold
            if max_signal > threshold:
                max_greater = True
        else:
            # If the user doesn't give us a reasonable value for threshold,
            # use the maximum signal from dbGroup for this gene
            max_log2 = max_signal

        intensity = 0
        log2 = math.log(2)

        n = 1
        low_alert = 0
        self.start_table("Relative")
        for group in self.groups:
            # Get control for this group
            control = group.get_control_signal(gene)

            for tissue in group.tissues:
                # If for some reason this tissue object doesn't have a color key
                # assigned (malformed XML data?), skip it
                if tissue.colorKey == (0, 0, 0):
                    continue

                (signal, stddev) = tissue.get_mean_signal(gene)

                if (signal is None) or (control is None) or (signal == 0) or (control == 0):
                    ratio_log2 = None
                    color = (221, 221, 221)  # Masked
                    low_alert = 1
                    self.append_table("Relative", tissue, ratio_log2, n, True, None, signal, control, tissue.colorKey)
                else:
                    ratio_log2 = math.log(signal / control) / log2

                    if max_log2 == 0:
                        intensity = 0
                    else:
                        intensity = int(math.floor(255 * (ratio_log2 / max_log2)))

                    intensity = self._clamp(intensity, -255, 255)

                    if signal <= 20 and grey_mask == "on":
                        low_alert = 1
                        color = (221, 221, 221)  # Masked: #CCCCCC
                    elif intensity > 0:
                        # Otherwise, color appropriately
                        color = (255, 255 - intensity, 0)
                    else:
                        color = (255 + intensity, 255 + intensity, - intensity)

                    # Set up Grey Mask alert
                    if signal <= 20 and grey_mask != "on":
                        low_alert = 1

                    self.append_table("Relative", tissue, ratio_log2, n, True, None, signal, control, tissue.colorKey)

                # Finally run replaceFill
                out_image.replaceFill(self.colorMap, tissue.colorKey, color)

            n += 1

        # Complete Table of Expression Values
        self.end_table()
        self.render_legend(out_image, "Log2 Ratio", max_log2, -max_log2, less_than=max_greater,
                           greater_than=max_greater, is_relative=True)
        return out_image, max_signal, max_signal_gene, max_signal_gene2, low_alert

    def render_comparison(self, gene1, gene2, threshold=0.0):
        """
        Renders each tissue according to the ratio of the signal strength
        of the first gene to its control relative to ratio of the signal
        strength of the second gene to its control
        :param gene1: Gene we're evaluating
        :param gene2: Base Gene, used as control
        :param threshold:
        :return: A PIL Image object containing the final rendered data
        """
        out_image = self.colorMap.copy()
        max_signal, max_signal_gene, max_signal_gene2 = self.get_view_max_signal("Compare", gene1, gene2)
        max_greater = False

        if threshold >= self.conf['minThreshold_Compare']:
            efp_max = threshold
            if max_signal > threshold:
                max_greater = True
        elif max_signal == 0:
            efp_max = self.conf['minThreshold_Compare']
        else:
            # If the user doesn't give us a reasonable value for threshold,
            # use the maximum signal from dbGroup for this gene
            efp_max = max_signal

        intensity = 0  # Cast as int
        log2 = math.log(2)

        n = 1
        self.start_table("Compare")
        for group in self.groups:
            # Get control data
            control1 = group.get_control_signal(gene1)
            control2 = group.get_control_signal(gene2)

            for tissue in group.tissues:
                # If for some reason this tissue object doesn't have a color key
                # assigned (malformed XML data?), skip it
                if tissue.colorKey == (0, 0, 0):
                    continue

                # Get the signal
                (sig1, stddev1) = tissue.get_mean_signal(gene1)
                (sig2, stddev2) = tissue.get_mean_signal(gene2)

                if (sig1 is None) or (sig2 is None) or (control1 is None) or (control2 is None) or \
                        (sig1 == 0) or (sig2 == 0) or (control1 == 0) or (control2 == 0):
                    color = (221, 221, 221)
                    self.append_table("Compare", tissue, None, n, True, None, None, None, tissue.colorKey)

                else:
                    ratio1_log = math.log(sig1 / control1) / log2
                    ratio2_log = math.log(sig2 / control2) / log2

                    if efp_max == 0:
                        intensity = 0
                    else:
                        intensity = int((ratio1_log - ratio2_log) * 255.0 / efp_max)

                    intensity = self._clamp(intensity, -255, 255)

                    if intensity is None:
                        color = (221, 221, 221)
                    elif intensity > 0:
                        # Map values above equal point to [yellow, red]
                        color = (255, 255 - intensity, 0)
                    else:
                        # Map values below equal point to [blue, yellow]
                        color = (255 + intensity, 255 + intensity, - intensity)

                    self.append_table("Compare", tissue, ratio1_log - ratio2_log, n, True, None, None, None, tissue.colorKey)

                out_image.replaceFill(self.colorMap, tissue.colorKey, color)
            n += 1

        # Complete Table of Expression Values
        self.end_table()

        self.render_legend(out_image, "Log2 Ratio", efp_max, -efp_max, less_than=max_greater, greater_than=max_greater,
                           is_relative=True)
        return out_image, max_signal, max_signal_gene, max_signal_gene2

    def draw_line(self, img, signal, off_set_val, displace_val, top, bottom, color):
        draw = PIL.ImageDraw.Draw(img)
        offset_x = signal * off_set_val
        offset_x += displace_val
        offset_x = int(offset_x)
        draw.line((offset_x, top, offset_x, bottom), fill=color)

    def draw_image(self, mode, max_signal_in_datasource, max_signal_gene, max_signal_gene2, gene1, gene2, img):
        # save generated image in output file.
        # First clean up output folder if necessary
        file_counter = 50
        files = glob.glob(self.conf['OUTPUT_FILES'] + "/*")

        if len(files) > file_counter:
            os.system("rm -f " + self.conf['OUTPUT_FILES'] + "/*")

        # Create a named temporary file with global read permissions
        outfile = tempfile.mkstemp(suffix='.png', prefix='efp-', dir=self.conf['OUTPUT_FILES'])
        os.system("chmod 644 " + outfile[1])

        # colors
        red = (255, 0, 0)  # FF0000
        blue = (0, 0, 255)  # 0000FF
        font_sizes = [(12, '10'), (15, '12'), (18, '14'), (21, '18'), (27, '24')]
        font_sizes_sm = [(12, '08'), (15, '10'), (18, '12'), (21, '14'), (27, '18'), (36, '24')]
        fontsize = '10'
        font_sizes_m = '08'

        for h, s in font_sizes:
            if self.legendSize >= h:
                fontsize = s
        for h, s in font_sizes_sm:
            if self.legendSize >= h:
                font_sizes_m = s

        # Draw the AGI in the top left corner
        draw = PIL.ImageDraw.Draw(img)
        font = PIL.ImageFont.load(self.conf['STATIC_FILES'] + "/pilfonts/helvB%s.pil" % fontsize)
        font_small = PIL.ImageFont.load(self.conf['STATIC_FILES'] + "/pilfonts/helvB%s.pil" % font_sizes_m)
        font_oblique = PIL.ImageFont.load(self.conf['STATIC_FILES'] + "/pilfonts/helvBO%s.pil" % font_sizes_m)
        color = (153, 153, 153)  # grey, #999999

        # draw a red box around button links to heatmaps if the searched gene is contained in the heatmap
        for extra in self.extras:
            # extra.button is true when a red box should be drawn, ie. when the searched gene is in the heatmap list
            if extra.button:
                # split the coords into a list and cast to integers
                str_coords = extra.coords.split(',')
                coords = []

                for coord in str_coords:
                    coords.append(int(coord))

                # eFP Tomato uses a * instead of a box
                if self.conf['species'] == 'TOMATO':
                    draw.text((coords[0], coords[1]), "*", font=font, fill=(0, 0, 0))
                else:
                    # draw the box using the coords list
                    draw.polygon(coords, outline=(255, 0, 0))

        for group in self.groups:
            for tissue in group.tissues:
                for coords in tissue.coords:
                    str_coords = coords.split(',')
                    int_coord_list = []
                    for intCoord in str_coords:
                        int_coord_list.append(int(intCoord))

        draw.text(self.conf['GENE_ID1_POS'], gene1.get_gene_id(), font=font, fill=color)

        if gene1.get_alias():
            draw.text(self.conf['GENE_ALIAS1_POS'], gene1.get_alias(), font=font_oblique, fill=color)

        if gene1.get_probeset_id():
            draw.text(self.conf['GENE_PROBESET1_POS'], gene1.get_probeset_id(), font=font_small, fill=color)

        if gene1.get_lookup():
            draw.text(self.conf['GENE_ORTHO1_POS'], gene1.get_lookup(), font=font, fill=color)

        if mode == 'Compare':
            draw.text(self.conf['GENE_ID2_POS'], gene2.get_gene_id(), font=font, fill=color)
            if gene2.get_alias():
                draw.text(self.conf['GENE_ALIAS2_POS'], gene2.get_alias(), font=font_oblique, fill=color)

            # underline gene ids to distinguish bars in chart
            draw.line((self.conf['GENE_ID1_POS'][0], self.conf['GENE_ID1_POS'][1] + 14,
                       self.conf['GENE_ID1_POS'][0] + 8 * len(gene1.get_gene_id()), self.conf['GENE_ID1_POS'][1] + 14),
                      fill=red)
            draw.line((self.conf['GENE_ID2_POS'][0], self.conf['GENE_ID2_POS'][1] + 14,
                       self.conf['GENE_ID2_POS'][0] + 8 * len(gene2.get_gene_id()), self.conf['GENE_ID2_POS'][1] + 14),
                      fill=blue)

            if gene2.get_probeset_id():
                draw.text(self.conf['GENE_PROBESET2_POS'], gene2.get_probeset_id(), font=font_small, fill=color)

        displace_y = int(self.graph[1])
        height = int(self.graph[3])
        bottom = displace_y

        data_source = self.database
        if data_source not in self.conf['GRAPH_SCALE_UNIT']:
            data_source = 'default'

        # show where the maximum signal in any data source lies on the little graph
        for (signal, color) in zip((max_signal_in_datasource, max_signal_gene, max_signal_gene2),
                                   ('#cccccc', red, blue)):

            # in case view_max_signal2 is not defined
            if signal is None:
                break

            displace_x = int(self.graph[0])  # (re)set x base coordinate for bar
            for scale_param in self.conf['GRAPH_SCALE_UNIT'][data_source]:
                if signal <= scale_param[0]:  # if signal is within range draw bar
                    self.draw_line(img, signal, scale_param[1], displace_x,
                                   displace_y - height * (
                                           self.conf['GRAPH_SCALE_UNIT'][data_source].index(scale_param) + 1) / (
                                           1.0 * len(self.conf['GRAPH_SCALE_UNIT'][data_source])),
                                   bottom, color)
                    break  # ... and go to next signal
                else:  # else adjust x base coordinates to next segment in graph
                    displace_x += scale_param[0] * scale_param[1]

        img.save(outfile[1])
        return outfile[1]

    def get_tissue_signal(self, probeset_id, sample_id):
        """
        get_tissue_signal
        :param probeset_id:
        :param sample_id:
        :return: Returns the tissue signal for a given gene ID
        """
        signal = None

        if self.conn is None:
            self.connect()

        cursor = self.conn.cursor()
        cursor.execute(
            "SELECT data_signal FROM " + self.conf['DB_DATA_TABLE'] + " WHERE data_probeset_id=%s AND data_bot_id=%s",
            (probeset_id, sample_id))
        row = cursor.fetchone()

        if row:
            signal = row[0]

        cursor.close()

        return signal

    def get_max_in_data_source(self, gene):
        """
        Returns the efp_max signal and the datasource it occurs in across all datasources of assigned group for a
        particular gene
        :param gene:
        :return:
        """
        # overall efp_max signal across all data sources
        overall_max = 0
        max_signal = 0
        max_data_source = ''
        max_sample = ''

        for data_source in self.conf['groupDatasource'][self.dbGroup]:
            if data_source == self.name:
                max_signal, max_sample = self.get_max_signal(gene)
            else:
                spec = Specimen(self.conf)
                spec.load("%s/%s.xml" % (self.conf['dataDir'], data_source))
                view = list(spec.get_views().values())[0]  # take first view in account for efp_max signal
                if gene is not None:
                    max_signal, max_sample = view.get_max_signal(gene)

            # Now get the data source name for that sample.
            # This is needed if multiple XML files share the same data source! Example: eFP Human
            # Otherwise eFP shows wrong data source for  max signal
            if (max_sample is not None) and (max_signal >= overall_max):
                overall_max = max_signal
                with open("%s/%s.xml" % (self.conf['dataDir'], data_source)) as f:
                    if max_sample in f.read():
                        max_data_source = data_source

        return overall_max, max_data_source

    def get_max_in_datasource_dict(self, gene):
        """
        This if for the max table
        @param gene:
        @return:
        """

        output = {}
        for data_source in self.conf['groupDatasource'][self.dbGroup]:
            if data_source == self.name:
                max_signal = self.get_max_signal(gene)
            else:
                spec = Specimen(self.conf)
                spec.load("%s/%s.xml" % (self.conf['dataDir'], data_source))
                view = list(spec.get_views().values())[0]  # take first view in account for max signal
                gene = view.alter_gene(gene)
                max_signal = view.get_max_signal(gene)

            # If max signal is None, eFP Crashes!
            if max_signal[0] is None or max_signal[1] is None:
                output[data_source] = (0, "")
            else:
                output[data_source] = max_signal

        return output

    def get_max_signal(self, gene):
        if self.conn is None:
            self.connect()

        # select the efp_max data signal in this datasource
        cursor = self.conn.cursor()
        cursor.execute("SELECT data_signal, data_bot_id FROM sample_data WHERE (data_probeset_id=%s OR data_probeset_id=%s) \
         ORDER BY data_signal DESC LIMIT 1", (gene.get_probeset_id(), gene.get_gene_id()))
        row = cursor.fetchone()

        if row is not None:
            efp_max = row[0]
            sample = row[1]
        else:
            efp_max = 0
            sample = None

        cursor.close()
        return efp_max, sample

    def connect(self):
        try:
            self.conn = MySQLdb.connect(host=self.conf['DB_HOST'], user=self.conf['DB_USER'],
                                        passwd=self.conf['DB_PASSWD'],
                                        db=self.database)
        except MySQLdb.Error as e:
            print("Error %d: %s" % (e.args[0], e.args[1]))

    def create_gene(self, gene_id):
        return efpDb.Gene(gene_id, self.database, self.conf)

    def alter_gene(self, gene):
        """
        Only used in eFP Maize
        :param gene:
        :return:
        """
        if self.conf["species"] == "MAIZE":
            if gene.gene_id is not None:

                if self.database == "maize_RMA_linear":
                    maize_convert1 = re.match("^GRMZM(2|5)G[0-9]{6}$", gene.gene_id)
                    maize_convert2 = re.match("^AC[0-9]{6}\.[0-9]{1}_FG[0-9]{3}$", gene.gene_id)

                    if maize_convert1 is not None:
                        gene.retrieve_lookup_gene_data(gene.gene_id)

                    if maize_convert2 is not None:
                        gene.gene_id = gene.gene_id.replace("FG", "FGT")

                if self.database == "maize_iplant" or self.database == "maize_ears" or self.database == "maize_leaf_gradient" or self.database == "maize_rice_comparison":
                    maize_convert3 = re.match("GRMZM(2|5)G[0-9]{6}_T[0-9]{1,2}", gene.gene_id)
                    maize_convert4 = re.match("^AC[0-9]{6}\.[0-9]{1}_FGT[0-9]{3}$", gene.gene_id)

                    if maize_convert3 is not None:
                        gene.gene_id = re.sub("_T[0-9]{1,2}", "", gene.gene_id)

                    if maize_convert4 is not None:
                        gene.gene_id = gene.gene_id.replace("FGT", "FG")

                if self.database == "maize_gdowns":
                    return gene

        else:
            gene.webservice_gene = gene.gene_id

        return gene

    def alter_webservice_gene(self, gene):
        """
        Only used in eFP Maize
        :param gene:
        :return:
        """
        if self.conf["species"] == "MAIZE":
            if gene.gene_id is not None:
                maize_convert1 = re.match("^GRMZM(2|5)G[0-9]{6}$", gene.gene_id)
                maize_convert2 = re.match("^AC[0-9]{6}\.[0-9]{1}_FG[0-9]{3}$", gene.gene_id)
                maize_sp1 = "_T01"
                if maize_convert1 is not None:
                    gene.webservice_gene = gene.gene_id + maize_sp1
                elif maize_convert2 is not None:
                    gene.webservice_gene = gene.gene_id.replace("FG", "FGT")
                else:
                    gene.webservice_gene = gene.gene_id
        else:
            gene.webservice_gene = gene.gene_id

        return gene


class SpecimenHandler(ContentHandler):
    def __init__(self, specimen, conf):
        super().__init__()
        self.spec = specimen
        self.ctrlSample = ""
        self.sampleDict = {}
        self.conf = conf

    def startElementNS(self, efp_dict, qname, attrs):
        uri, name = efp_dict
        if name == 'view':
            view_class = attrs.get('class')
            if view_class is None:
                view_class = "View"
            self.sampleDict = {}
            exec(
                "self.currentView = %s(attrs.getValueByQName('name'), attrs.getValueByQName('db'), attrs.getValueByQName('dbGroup'), '%s/' + attrs.getValueByQName('img'), self.conf)" % (
                    view_class, self.conf['dataDir']))

        if name == 'coords':
            graph = (
                attrs.getValueByQName('graphX'), attrs.getValueByQName('graphY'), attrs.getValueByQName('graphWidth'),
                attrs.getValueByQName('graphHeight'))
            legend = (int(attrs.getValueByQName('legendX')), int(attrs.getValueByQName('legendY')))
            legend_size = 12
            try:
                legend_size = int(attrs.getValueByQName('legend_size'))
                if legend_size < 8:
                    legend_size = 8
            except KeyError:
                pass
            self.currentView.add_graph_coords(graph)
            self.currentView.add_legend_coords(legend)
            self.currentView.legendSize = legend_size

        if name == 'extra':
            # These two are optional and not in all XML files. So try first
            check_value = None
            check_column = None

            try:
                check_value = attrs.getValueByQName('check')
            except KeyError:
                check_value = None

            try:
                check_column = attrs.getValueByQName('checkColumn')
            except KeyError:
                check_column = None

            e = Extra(attrs.getValueByQName('name'), attrs.getValueByQName('link'), attrs.getValueByQName("parameters"),
                      attrs.getValueByQName('coords'), check_value, check_column)
            self.currentView.add_extra(e)

        if name == 'group':
            self.currentGroup = Group(attrs.getValueByQName('name'))

        if name == 'control':
            sample_name = attrs.getValueByQName('sample')
            ctrl_sample = self.sampleDict.get(sample_name)
            if ctrl_sample is None:
                ctrl_sample = Sample(sample_name, self.currentView)
                self.sampleDict[sample_name] = ctrl_sample
            self.currentGroup.add_ctrl_sample(ctrl_sample)

        if name == 'tissue':
            t = Tissue(attrs.getValueByQName('name'), attrs.getValueByQName('colorKey'))
            self.currentTissue = t

        if name == 'link':
            url = attrs.getValueByQName('url')
            self.currentTissue.add_url(url)

        if name == 'area':
            coords = attrs.getValueByQName('coords')
            self.currentTissue.add_coords(coords)

        if name == 'sample':
            sample_name = attrs.getValueByQName('name')
            sample = self.sampleDict.get(sample_name)
            if sample is None:
                sample = Sample(sample_name, self.currentView)
                self.sampleDict[sample_name] = sample
            self.currentTissue.add_sample(sample)

    def endElementNS(self, qname, name):
        if name == 'view':
            self.spec.add_view(self.currentView)

        if name == 'group':
            self.currentView.add_group(self.currentGroup)

        if name == 'tissue':
            self.currentGroup.add_tissue(self.currentTissue)


class Specimen:
    def __init__(self, conf):
        self.conf = conf
        self.views = {}  # Dictionary of views

    def add_view(self, view):
        self.views[view.name] = view

    def get_view(self, name):
        return self.views[name]

    def get_views(self):
        return self.views

    def load(self, file):
        # Create the handler
        handler = SpecimenHandler(self, self.conf)

        # Parse the file
        tree = etree.parse(file)
        lxml.sax.saxify(tree, handler)


if __name__ == '__main__':
    pass
