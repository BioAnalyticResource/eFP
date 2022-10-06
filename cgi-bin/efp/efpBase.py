#!/usr/bin/python3
"""
Created on Jan 19, 2017

Code in this file was extracted from efpWeb.cgi to remove business logic
from the presentation logic.

@author: Andrew Robinson
@author: Robert Breit
"""

import os
import re
import tempfile

from . import efp, efpService, efpDb


def process_request(dataSource, primaryGene, secondaryGene, threshold, ncbi_gi, mode, useThreshold, grey_low,
                    grey_stddev, nav_bar, conf):
    """
    Handles request for main web view and returns a dictionary with all values
    to be displayed
    @param conf: efp_dict object of configuration variables
    @return: efp_dict object with values to place on website
    """
    error = 0
    error_strings = []
    alert_strings = []
    test_alert_strings = []
    low_alert = 0
    sd_alert = 0

    img = None
    img_map = {}
    img_filename = {}
    table_file = {}
    views = []
    max_dict = {}

    webservice_gene1 = None
    webservice_gene2 = None
    gene1 = None
    gene2 = None

    # Default gene id
    if primaryGene is None and ncbi_gi:
        # TODO: This maybe a bug here.
        primaryGene = efpDb.ncbi_to_gene_id(ncbi_gi)
        if primaryGene is None:
            primaryGene = ncbi_gi
            error_str = 'The requested NCBI gi "%s" doesn\'t correspond to a given AGI.<br>' % ncbi_gi
            error_strings.append(error_str)
    elif primaryGene is None:
        primaryGene = conf['GENE_ID_DEFAULT1']
    if secondaryGene is None:
        secondaryGene = conf['GENE_ID_DEFAULT2']

    if useThreshold == "":
        useThreshold = None

    # Try Entered Threshold; if fails or threshold not checked use default threshold
    if useThreshold:
        try:
            threshold = float(threshold)  # Convert str to float
        except FloatingPointError:
            # Threshold string was malformed
            error = 1
            error_str = 'Invalid Threshold Value "%s"<br>' % threshold
            error_strings.append(error_str)
            useThreshold = None
    if useThreshold is None and threshold is None:
        # assign a default value for first calls
        if mode == "Relative" or mode == "Compare":
            threshold = 2.0
        else:  # Absolute or none
            threshold = 500
        first_call = 1
    else:
        threshold = float(threshold)
        first_call = 0

    if dataSource is None:
        dataSource = conf['defaultDataSource']

    # Serialize data from XML file into a Specimen object
    spec = efp.Specimen(conf)
    xml_name = "%s/%s.xml" % (conf['dataDir'], dataSource)
    spec.load(xml_name)

    # Right now the browser only has one view - "all"
    # In the future, there should be a drop down menu letting users
    # choose multiple views

    default_img_filename = "%s/%s.png" % (conf['dataDirWeb'], dataSource)
    if mode is None:
        # If no mode is selected (99% of the time this means the user just arrived
        # at the page), just show them a color map
        # Set Developmental_Map as default DataSource
        if dataSource is None:
            dataSource = conf['defaultDataSource']
    else:
        # all view
        for name, view in spec.get_views().items():
            # If either of these probe IDs are None (bad inputs), then we just
            # spit out the default image again
            gene1 = view.create_gene(primaryGene)
            gene2 = view.create_gene(secondaryGene)
            gene1 = view.alter_gene(gene1)
            gene2 = view.alter_gene(gene2)
            webservice_gene1 = view.alter_webservice_gene(gene1)
            webservice_gene2 = view.alter_webservice_gene(gene2)
            if gene1.get_gene_id() is None:
                # TODO: A bug maybe here
                # if gene1.get_probeset_id() == None:
                error_str = 'The requested Primary gene / probeset ID "%s" cannot be found in %s datasource ' % (
                    primaryGene, view.dbGroup)
                error = 1
                error_strings.append(error_str)
            elif mode == 'Compare' and gene2.get_gene_id() is None:
                # elif mode == 'Compare' and gene2.get_probeset_id() == None:
                error = 1
                error_str = 'The requested Secondary gene / probeset ID "%s" cannot be found in %s datasource <br>' % (
                    secondaryGene, view.dbGroup)
                error_strings.append(error_str)
            elif primaryGene == secondaryGene and mode == 'Compare':
                error = 1
                error_str = 'The requested Secondary gene / probeset ID "%s" must be different than the Primary ID<br>' \
                            % secondaryGene
                error_strings.append(error_str)
                view_max_signal = 2.0
            else:
                if mode == 'Absolute':
                    if useThreshold:
                        (img, view_max_signal, view_max_signal1, view_max_signal2, sd_alert) = view.render_absolute(
                            gene1, threshold, grey_mask=grey_stddev)
                    else:
                        (img, view_max_signal, view_max_signal1, view_max_signal2, sd_alert) = view.render_absolute(
                            gene1, grey_mask=grey_stddev)
                elif mode == 'Relative':
                    if useThreshold:
                        (img, view_max_signal, view_max_signal1, view_max_signal2, low_alert) = view.render_relative(
                            gene1, threshold, grey_mask=grey_low)
                    else:
                        (img, view_max_signal, view_max_signal1, view_max_signal2, low_alert) = view.render_relative(
                            gene1, grey_mask=grey_low)
                elif mode == 'Compare':
                    if useThreshold:
                        (img, view_max_signal, view_max_signal1, view_max_signal2) = view.render_comparison(
                            gene1, gene2, threshold)
                    else:
                        (img, view_max_signal, view_max_signal1, view_max_signal2) = view.render_comparison(gene1, gene2)

                # find the efp_max signal across all data sources and provide a link to that data source
                (max_signal_in_data_source, max_data_source) = view.get_max_in_data_source(gene1)

                max_signal_in_data_source = round(max_signal_in_data_source, 2)
                max_data_source = re.sub(r"_", " ", max_data_source)
                alert_str = "For %s data, this probe set reaches its maximum expression level (expression potential)" \
                            " of <b>%s</b> in the <b>%s</b> data source." % (
                                view.dbGroup, max_signal_in_data_source, max_data_source)
                alert_strings.append(alert_str)

                # find the max accross all sourse
                if nav_bar:
                    max_dict = view.get_max_in_datasource_dict(gene1)

                # alert the user that the scale has changed if no threshold is set
                if useThreshold is None and first_call != 1:
                    if view_max_signal > threshold:
                        use_threshold_flag = "on"
                        threshold_level_suggested = max_signal_in_data_source
                        if mode == 'Relative':
                            threshold_level_suggested = 4
                        if mode == 'Compare':
                            threshold_level_suggested = 4
                        alert_str = "For %s data, note the maximum signal value has increased to %s from %s. Use the " \
                                    "<a href='efpWeb.cgi?dataSource=%s&mode=%s&primaryGene=%s&secondaryGene=%s&useThreshold=%s&threshold=%s&grey_low=%s&grey_stddev=%s'>Signal Threshold option to keep it constant at %s</a>, or enter a value in the Signal Threshold box, such as <a href='efpWeb.cgi?dataSource=%s&mode=%s&primaryGene=%s&secondaryGene=%s&useThreshold=%s&threshold=%s&grey_low=%s&grey_stddev=%s'>%s</a>. The same color scheme will then be applied across all views.<br>" % (
                                        view.dbGroup, view_max_signal, threshold, dataSource, mode, primaryGene,
                                        secondaryGene,
                                        use_threshold_flag, threshold, grey_low, grey_stddev, threshold, dataSource, mode,
                                        primaryGene, secondaryGene, use_threshold_flag, threshold_level_suggested,
                                        grey_low, grey_stddev, threshold_level_suggested)
                        alert_strings.append(alert_str)
                    elif view_max_signal < threshold:
                        use_threshold_flag = "on"
                        threshold_level_suggested = max_signal_in_data_source

                        if mode == 'Relative':
                            threshold_level_suggested = 4

                        if mode == 'Compare':
                            threshold_level_suggested = 4

                        alert_str = "For %s data, note the maximum signal value has decreased to %s from %s. Use the " \
                                    "<a href='efpWeb.cgi?dataSource=%s&mode=%s&primaryGene=%s&secondaryGene=%s&useThreshold=%s&threshold=%s&grey_low=%s&grey_stddev=%s'>Signal Threshold option to keep it constant at %s</a>, or enter a value in the Signal Threshold box, such as <a href='efpWeb.cgi?dataSource=%s&mode=%s&primaryGene=%s&secondaryGene=%s&useThreshold=%s&threshold=%s&grey_low=%s&grey_stddev=%s'>%s</a>. The same color scheme will then be applied across all views.<br>" % (
                                        view.dbGroup, view_max_signal, threshold, dataSource, mode, primaryGene,
                                        secondaryGene, use_threshold_flag, threshold, grey_low, grey_stddev, threshold,
                                        dataSource, mode, primaryGene, secondaryGene, use_threshold_flag,
                                        threshold_level_suggested, grey_low, grey_stddev, threshold_level_suggested)
                        alert_strings.append(alert_str)
                    else:
                        alert_str = ""
                    threshold = view_max_signal
                elif useThreshold is None and first_call == 1:
                    threshold = view_max_signal

                # alert the user if SD filter or low filter should be activated
                if grey_stddev != "on" and sd_alert == 1 and mode == 'Absolute':
                    grey_stddev_flag = "on"
                    if useThreshold is None:
                        useThreshold = ""
                    alert_str = "Some samples exhibit high standard deviations for replicates. You can use " \
                                "<a href='efpWeb.cgi?dataSource=%s&mode=%s&primaryGene=%s&secondaryGene=%s&useThreshold=%s&threshold=%s&grey_low=%s&grey_stddev=%s'>standard deviation filtering</a> to mask those with a deviation greater than half their expression value.<br>" % (
                                    dataSource, mode, primaryGene, secondaryGene, useThreshold, threshold, grey_low,
                                    grey_stddev_flag)
                    alert_strings.append(alert_str)

                # alert the user if SD filter or low filter should be activated
                if grey_low != "on" and low_alert == 1 and mode == 'Relative':
                    grey_low_flag = "on"
                    if useThreshold is None:
                        useThreshold = ""
                    alert_str = "Some sample ratios were calculated with low values that exhibit higher variation, " \
                                "potentially leading to ratios that are not a good reflection of the biology. You can <a href='efpWeb.cgi?dataSource=%s&mode=%s&primaryGene=%s&secondaryGene=%s&useThreshold=%s&threshold=%s&grey_low=%s&grey_stddev=%s'>low filter below 20 units</a> to mask these.<br>" % (
                                    dataSource, mode, primaryGene, secondaryGene, useThreshold, threshold,
                                    grey_low_flag, grey_stddev)
                    alert_strings.append(alert_str)

                # Otherwise, we render and display the option
                img_map[view.name], zero_gene = view.get_image_map(mode, gene1, gene2, useThreshold, threshold, dataSource,
                                                        grey_low, grey_stddev)

                # The zero expression genes of eFP Human and eFP Tomato
                if zero_gene:
                    alert_strings.append('<font color="red">Some samples for this gene have expression values of 0.<br>')

            if img:
                img_filename[view.name] = view.draw_image(mode, max_signal_in_data_source, view_max_signal1,
                                                          view_max_signal2, gene1, gene2, img)
                # Create a table of Expression Values and save it in a temporary file
                exp_table = view.table
                table_file[view.name] = tempfile.mkstemp(suffix='.html', prefix='efp-', dir=conf["OUTPUT_FILES"])
                os.system("chmod 644 " + table_file[view.name][1])
                tf = open(table_file[view.name][1], 'w')
                tf.write(exp_table)
                tf.close()
                chart_file = table_file[view.name][1].replace(".html", ".png")
                view.save_chart(chart_file, mode)
                views.append(view.name)

    if mode and error == 0:
        # process links
        info = efpService.Info()
        if info.load("%s/efp_info.xml" % conf["dataDir"]) is None:
            services = []
            for name in (info.get_services()):
                s = {'name': name}
                service = info.get_service(name)
                s['service'] = service
                s['external'] = service.get_external()
                highlight1 = service.check_service(webservice_gene1.webservice_gene)
                highlight2 = None
                if mode == 'Compare':
                    highlight2 = service.check_service(webservice_gene2.webservice_gene)
                if highlight1 == 'error' or highlight2 == 'error':
                    continue
                elif highlight1:
                    link = service.get_link(webservice_gene1.webservice_gene)
                    gene = webservice_gene1.webservice_gene
                elif highlight2:
                    link = service.get_link(webservice_gene2.webservice_gene)
                    gene = webservice_gene1.webservice_gene
                else:
                    continue
                s['gene'] = gene
                s['link'] = link
                services.append(s)

        # process views
        view_no = 1
        views_proc = []
        for view_name in views:
            view = {'view_no': view_no, 'view_name': view_name, 'imgFile': img_filename[view_name]}
            view['imgFileWeb'] = "%s/%s" % (conf["OUTPUT_FILES_WEB"], os.path.basename(view['imgFile']))
            if view_name in img_map:
                view['img_map'] = img_map[view_name]
            view['last_element'] = os.path.basename(img_filename[view_name])
            match = re.search(r"^efp", view['last_element'])
            if match is not None:
                view['imgFile'] = conf["OUTPUT_FILES_WEB"] + '/%s' % (view['last_element'])
            view['tableFile_name'] = conf["OUTPUT_FILES_WEB"] + '/%s' % (os.path.basename(table_file[view_name][1]))
            view['tableChart_name'] = view['tableFile_name'].replace(".html", ".png")
            view['popup_width'] = '1000'
            view['bg_color'] = '#FFFFFF'

            views_proc.append(view)
            view_no += 1

    # make list of data sources
    data_sources = []
    xml = find_xml(conf["dataDir"])
    for x in sorted(xml):
        src = (x, x.replace('_', ' '))
        data_sources.append(src)

    return default_img_filename, error, error_strings, test_alert_strings, alert_strings, webservice_gene1, \
           webservice_gene2, views, img_filename, img_map, table_file, gene1, gene2, dataSource, primaryGene, \
           secondaryGene, threshold, ncbi_gi, mode, useThreshold, grey_low, grey_stddev, max_dict


def find_xml(data_dir):
    """
    Find xml files in data directory
    @param data_dir:
    @return:
    """
    xml = []
    files = os.listdir(data_dir)
    for f in files:
        if f.endswith(".xml") and not f.startswith("efp_"):
            xml.append(f[0:-4])
    return xml


def html_color_to_rgb(color_string):
    if color_string == "":
        return 0, 0, 0
    color_string = color_string.strip()
    if color_string[0] == '#':
        color_string = color_string[1:]
    if len(color_string) != 6:
        raise ValueError("input #%s is not in #RRGGBB format" % color_string)
    r, g, b = color_string[:2], color_string[2:4], color_string[4:]
    r, g, b = [int(n, 16) for n in (r, g, b)]
    return r, g, b


def rgb_to_html_color(rgb):
    if len(rgb) != 3:
        raise ValueError("input %s is not in rgb format" % rgb)
    color_string = '#%02X%02X%02X' % (rgb[0], rgb[1], rgb[2])
    return color_string


def rgb_to_gray(rgb):
    if len(rgb) != 3:
        raise ValueError("input %s is not in rgb format" % rgb)
    gray = rgb[0] * 0.299 + rgb[1] * 0.587 + rgb[2] * 0.114
    return gray


def clean_exit(message):
    print('Content-Type: text/html\n')
    print('<html><head></head><body>')
    print('An error has occurred: ' + message + '<br>')
    print('If your data is valid, please contact us.')
    print('</body></html')
    exit()
