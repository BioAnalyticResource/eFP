#!/usr/bin/python3
import cgi
import operator
import os
import re
import sys
import warnings
from efp import efpConfig, efpService, efpBase

warnings.filterwarnings("ignore", category=DeprecationWarning)

form = cgi.FieldStorage(keep_blank_values=1)

# Retrieve cgi inputs
dataSource = form.getvalue("dataSource")
primaryGene = form.getvalue("primaryGene")
secondaryGene = form.getvalue("secondaryGene")
ncbi_gi = form.getvalue("ncbi_gi")
threshold = form.getvalue("threshold")
mode = form.getvalue("mode")
mode_input = form.getvalue("modeInput")
useThreshold = form.getvalue("useThreshold")
grey_low = form.getvalue("grey_low")
grey_stddev = form.getvalue("grey_stddev")
nav_bar = form.getvalue("navbar")

if mode_input is not None:
    if mode_input != "":
        mode = mode_input

# Validate CGI inputs:
if dataSource and re.search(r"^[\d\D\s\-_]{0,48}$", dataSource) is None:
    efpBase.clean_exit("Data Source is invalid")

if primaryGene and re.search(efpConfig.inputRegEx, primaryGene, re.I) is None:
    efpBase.clean_exit("Primary Gene is invalid")

if secondaryGene and re.search(efpConfig.inputRegEx, secondaryGene, re.I) is None:
    efpBase.clean_exit("Secondary Gene is invalid")

if ncbi_gi and re.search(r"^\d{0,16}$", ncbi_gi) is None:
    efpBase.clean_exit("NCBI GI is invalid.")

if threshold and re.search(r"^\d{0,16}\.*\d*$", threshold) is None:
    efpBase.clean_exit("Threshold is invalid.")

if mode and re.search(r"^Absolute$|^Relative$|^Compare$", mode) is None:
    efpBase.clean_exit("Mode is invalid.")

if useThreshold and re.search(r"^on$", useThreshold) is None:
    efpBase.clean_exit("Use Threshold is invalid.")

if grey_low and re.search(r"^on$|^None$", grey_low) is None:
    efpBase.clean_exit("Grey low is invalid.")

if grey_stddev and re.search(r"^on$|^None$", grey_stddev) is None:
    efpBase.clean_exit("Grey low is invalid.")

if nav_bar == "0":
    nav_bar = False
else:
    nav_bar = True


# Fix soybean IDs
# Add a .1 at the end of the end of it is soybean_senescence, and remove if it doesn't
if efpConfig.species == "SOYBEAN":
    if dataSource == "soybean_senescence":
        if (primaryGene is not None) and (re.search(r"\.\d$", primaryGene) is None):
            primaryGene += ".1"
        if (secondaryGene is not None) and (re.search(r"\.\d$", secondaryGene) is None):
            secondaryGene = secondaryGene + ".1"
    else:
        if (primaryGene is not None) and re.search(r'\.\d$', primaryGene):
            primaryGene = re.sub(r'\.\d$', '', primaryGene)
        if (secondaryGene is not None) and re.search(r'\.\d$', secondaryGene):
            secondaryGene = re.sub(r'\.\d$', '', secondaryGene)

# Fix for Marchantia
if efpConfig.species == "MARCHANTIA":
    if (primaryGene is not None) and (re.search(r"\.\d$", primaryGene) is None):
        primaryGene += ".1"
    if (secondaryGene is not None) and (re.search(r"\.\d$", secondaryGene) is None):
        secondaryGene = secondaryGene + ".1"

# Fix eFP Wheat IDs
if efpConfig.species == "WHEAT":
    # This makes eFP Work for both v1.0 and v1.1 genes
    if primaryGene is not None:
        if dataSource in ["Wheat_Abiotic_Stress", "Wheat_Meiosis"]:
            # V1.1 only
            if re.search(r"(.+1G.+)", primaryGene):
                primaryGene = re.sub(r"1G", "2G", primaryGene)
            if re.search(r"(\d\d\d$)", primaryGene):
                primaryGene = primaryGene + ".1"
        else:
            # V1.0 only
            if re.search(r"(.+2G.+)", primaryGene):
                primaryGene = re.sub(r"2G", "1G", primaryGene)
            if re.search(r"(\.\d$)", primaryGene):
                primaryGene = re.sub(r"\.\d$", "", primaryGene)

    if secondaryGene is not None:
        if dataSource in ["Wheat_Abiotic_Stress", "Wheat_Meiosis"]:
            # V1.1 only
            if re.search(r"(.+1G.+)", secondaryGene):
                secondaryGene = re.sub(r"1G", "2G", secondaryGene)
            if re.search(r"(\d\d\d$)", secondaryGene):
                secondaryGene = secondaryGene + ".1"
        else:
            # V1.0 only
            if re.search(r"(.+2G.+)", secondaryGene):
                secondaryGene = re.sub(r"2G", "1G", secondaryGene)
            if re.search(r"(\.\d$)", secondaryGene):
                secondaryGene = re.sub(r"\.\d$", "", secondaryGene)

# Open test file here if required
# testing_efp = open(efpConfig.OUTPUT_FILES + "/testing_efpweb.txt", "w")

# create a dictionary from efpConfig module
CONF = {}
for k in dir(efpConfig):
    if not k.startswith('_'):
        CONF[k] = getattr(efpConfig, k)

CONF['webservice_gene1'] = None
CONF['webservice_gene2'] = None

# process request
(default_img_filename, error, error_strings, test_alert_strings, alert_strings, webservice_gene1, webservice_gene2,
 views, img_filename, img_map, table_file, gene1, gene2, dataSource, primaryGene, secondaryGene, threshold, ncbi_gi,
 mode, useThreshold, grey_low, grey_stddev, max_dict) = \
    efpBase.process_request(dataSource, primaryGene, secondaryGene, threshold, ncbi_gi, mode, useThreshold, grey_low,
                            grey_stddev, nav_bar,  CONF)

# HTML header
print('Content-Type: text/html\n')

# HTML code
print('<!DOCTYPE html>')
print('<html lang="en">')
print('<head>')
print('  <title>%s eFP Browser</title>' % efpConfig.spec_names[efpConfig.species])
print('  <meta http-equiv="Content-Type" content="text/html; charset=utf-8">')
print(
    '  <meta name="keywords" content="%s, genomics, expression profiling, mRNA-seq, Affymetrix, microarray, protein-protein interactions, protein structure, polymorphism, subcellular localization, proteomics, poplar, rice, Medicago, barley, transcriptomics, proteomics, bioinformatics, data analysis, data visualization, AtGenExpress, PopGenExpress, cis-element prediction, coexpression analysis, Venn selection, molecular biology">' % \
    efpConfig.spec_names[efpConfig.species])
print('  <link rel="stylesheet" type="text/css" href="%s/efp.css"/>' % efpConfig.STATIC_FILES_WEB)
print('  <link rel="stylesheet" type="text/css" href="%s/domcollapse.css"/>' % efpConfig.STATIC_FILES_WEB)
print('  <script src="%s/efp.js"></script>' % efpConfig.STATIC_FILES_WEB)
print('  <script>')
print('    regId = /%s/i;' % efpConfig.inputRegEx)
print('  </script>')
print('  <script src="%s/domcollapse.js"></script>' % efpConfig.STATIC_FILES_WEB)
print('</head>')
print('')
print('<body>')
print('  <form action="efpWeb.cgi" name="efpForm" method="POST" onSubmit="return checkAGIs()">')
print('  <table width="1100px" border="0" align="center" cellspacing="1" cellpadding="0">')
print('    <tr>')
print('      <td>')
print('        <span style="float:right; top:0px; left:538px; width:250px; height:75px;">')
print('          <script src="%s/popup.js"></script>' % efpConfig.STATIC_FILES_WEB)
print('        </span>')
print("        <h1 style='vertical-align:middle;'><a href='//bar.utoronto.ca'><img src='//bar.utoronto.ca/bbc_logo_small.gif' alt='To the Bio-Array Resource Homepage' border=0 align=absmiddle></a>&nbsp;<img src='//bar.utoronto.ca/bar_logo.gif' alt='The Bio-Array Resource' border=0 align=absmiddle>&nbsp;<img src='//bar.utoronto.ca/images/eFP_logo_large.png' align=absmiddle border=0>&nbsp;%s eFP Browser<br><img src='//bar.utoronto.ca/images/green_line.gif' width=98%% alt='' height='6px' border=0></h1>" % \
    efpConfig.spec_names[efpConfig.species])
print('      </td>')
print('   </tr>')
print('<tr><td align="middle">')
print('    <table>')
print('      <tr align = "center"><th>Data Source</th>')
print('      <th>Mode')
print(
    '<input type="checkbox" name="grey_stddev" title="In Absolute Mode, check to mask samples that exhibit a standard deviation of more than 50 percent of the signal value" ')
if grey_stddev == "on":
    print('checked')
print(' value="on" />')

print(
    '<input type="checkbox" name="grey_low" title="In Relative Mode, check to mask the use of low expression values in ratio calculations" ')
if grey_low == "on":
    print('checked')
print(' value="on" />')

if efpConfig.species == "CAMELINA":
    print('</th><th>Primary Gene ID<br>(C. sativa or Ath IDs)</th><th>Secondary Gene ID<br>(C. sativa or Ath IDs)</th>')
else:
    print('</th><th>Primary Gene ID</th><th>Secondary Gene ID</th>')

print(
    '      <th id="t1">Signal Threshold<input type="checkbox" name="useThreshold" title="Check to enable threshold" onclick="checkboxClicked(\'useThreshold\');" ')
if useThreshold == "on":
    print('checked')
print(' value="on" />')
print('</th><th></th></tr>')
print('      <tr><td>')

# Help Link
print(
    '      <img src="//bar.utoronto.ca/affydb/help.gif" border=0 align="top" alt="Click here for instructions in a new window" onClick="HelpWin = window.open(\'//bar.utoronto.ca/affydb/BAR_instructions.html#efp\', \'HelpWindow\', \'width=600,height=300,scrollbars,resizable=yes\'); HelpWin.focus();">&nbsp;')

# Build drop down list of Data Sources
if mode is None:
    print(
        '<select name="dataSource" onchange="location.href=\'efpWeb.cgi?dataSource=\' + this.options[this.selectedIndex].value ;">')
elif useThreshold is None:
    thresholdSwitch = ""
    print(
        '      <select name="dataSource" onchange="location.href=\'efpWeb.cgi?dataSource=\' + this.options[this.selectedIndex].value + \'&mode=%s&primaryGene=%s&secondaryGene=%s&useThreshold=%s&threshold=%s&grey_low=%s&grey_stddev=%s\' ;">' % (
            mode, primaryGene, secondaryGene, thresholdSwitch, threshold, grey_low, grey_stddev))
else:
    print(
        '      <select name="dataSource" onchange="location.href=\'efpWeb.cgi?dataSource=\' + this.options[this.selectedIndex].value + \'&mode=%s&primaryGene=%s&secondaryGene=%s&useThreshold=%s&threshold=%s&grey_low=%s&grey_stddev=%s\' ;">' % (
            mode, primaryGene, secondaryGene, useThreshold, threshold, grey_low, grey_stddev))

xml = efpBase.find_xml(efpConfig.dataDir)
for x in sorted(xml):
    print('    <option value="%s"' % x)
    # To preserve modes between form submits
    if dataSource == x:
        print('selected')
    x_text = x.replace('_', ' ')
    print('>%s</option>' % x_text)
print('      </select></td>')

# Build drop down list of modes
if mode is None:
    print('      <td><select selected="Absolute" name="mode" onchange="changeMode(\'mode\');">')
else:
    print(
        '		 <td><select selected="Absolute" name="mode" onchange="location.href=\'efpWeb.cgi?dataSource=%s&mode=\' + this.options[this.selectedIndex].text + \'&primaryGene=%s&secondaryGene=%s&grey_low=%s&grey_stddev=%s\' ">' % (
            dataSource, primaryGene, secondaryGene, grey_low, grey_stddev))

# Preserve mode between form submits. If the user selected 'Compare' as his/her
# mode, when the page reloads, the list should still have 'Compare' selected.
if mode == 'Relative':
    print('    <option>Absolute</option>')
    print('    <option selected>Relative</option>')
    # Skip Compare mode on eFP Maize
    if dataSource not in ["maize_leaf_gradient", "maize_rice_comparison", "rice_leaf_gradient", "rice_maize_comparison"]:
        print('    <option>Compare</option>')
elif mode == 'Compare':
    print('    <option>Absolute</option>')
    print('    <option>Relative</option>')
    # Skip Compare mode on eFP Maize
    if dataSource not in ["maize_leaf_gradient", "maize_rice_comparison", "rice_leaf_gradient", "rice_maize_comparison"]:
        print('    <option selected>Compare</option>')
else:  # Default (Absolute)
    print('    <option selected>Absolute</option>')
    print('    <option>Relative</option>')
    # Skip Compare mode on eFP Maize
    if dataSource not in ["maize_leaf_gradient", "maize_rice_comparison", "rice_leaf_gradient", "rice_maize_comparison"]:
        print('    <option>Compare</option>')

print('      </select></td><td>')
print('      <input type="text" id="g1" name="primaryGene" value="%s" size=17/></td><td>' % primaryGene)
print('      <input type="text" id="g2" name="secondaryGene" size=17 value="%s" ' % secondaryGene)

if mode != 'Compare' or dataSource == "rice_leaf_gradient" and dataSource == "rice_maize_comparison":
    # if mode != 'Compare':
    print('disabled')
print('      /></td><td>')

print('      <input type="text" id="t0" name="threshold" value="%s" ' % threshold)
if useThreshold is None:
    print('disabled')
print('      /></td>')
print('      <td><input type="submit" value="Go"/></td></tr>')
print('    </table>')
print('</td></tr>')

print('<tr><td>')
if error:
    print('    <ul>')
    for row in error_strings:
        print('<li class="error">%s</li>' % row)
    print('    </ul>')

if len(test_alert_strings) > 0:
    print('    <ul>')
    for row in test_alert_strings:
        print('<li>%s</li>' % row)
    print('    </ul>')

# print additional header text if configured for selected data source
if dataSource in efpConfig.datasourceHeader:
    print('%s' % efpConfig.datasourceHeader[dataSource])
elif 'default' in efpConfig.datasourceHeader:
    print('%s' % efpConfig.datasourceHeader['default'])

if len(alert_strings) > 0:
    print('    <ul>')
    for row in alert_strings:
        print('<li>%s</li>' % row)

    # eFP Wheat stuff
    if efpConfig.species == "WHEAT":
        homoeologues = gene1.get_homoeologues()

        if len(homoeologues) > 0:
            print('<li>Please see homoeologues of primary gene: ')

            for homoeolog in homoeologues:
                if homoeolog == gene1.gene_id:
                    pass
                else:
                    print("<a href='efpWeb.cgi?dataSource=%s&mode=%s&primaryGene=%s&secondaryGene=%s&override=%s&threshold=%s&modeMask_low=%s'>%s</a> " % (
                            dataSource, mode, homoeolog, secondaryGene, useThreshold, threshold, grey_low, homoeolog))

            print('</li>')

    print('    </ul>')
print('</td></tr>')

if mode and error == 0:
    # check external services
    # Serialize services data from XML file into a Info object
    info = efpService.Info()
    if info.load("%s/efp_info.xml" % efpConfig.dataDir) is None:
        print('<tr><td>')
        print('<table style="margin-left:auto;margin-right:auto"><tr>')
        for name in (info.get_services()):
            service = info.get_service(name)
            external = service.get_external()

            # eFP Maize stuff
            if name == "Metabolite" or name == "Enzyme":
                MaizeConvert3 = re.match(r"GRMZM(2|5)G[0-9]{6}_T[0-9]{1,2}", webservice_gene1.webservice_gene)
                if MaizeConvert3 is not None:
                    webservice_gene1.webservice_gene = re.sub(r"_T[0-9]{1,2}", "", webservice_gene1.webservice_gene)

            highlight1 = service.check_service(webservice_gene1.webservice_gene)
            highlight2 = None

            if mode == 'Compare':
                highlight2 = service.check_service(webservice_gene2.webservice_gene)
            if highlight1 == 'error' or highlight2 == 'error':
                print(
                    '<td><img title="connection error for service %s" width="50" height="50" alt="connection error" src="%s/error.png"></td>' % (
                        name, efpConfig.dataDirWeb))
                continue
            elif highlight1:
                link = service.get_link(webservice_gene1.webservice_gene)
                gene = webservice_gene1.webservice_gene
            elif highlight2:
                link = service.get_link(webservice_gene2.webservice_gene)
                gene = webservice_gene1.webservice_gene
            else:
                print(
                    '<td><img title="No %s data found" width="50" height="50" alt="No %s data found" style="opacity:0.30;filter:alpha(opacity=30);" src="%s/%s"></td>' % (
                        name, name, efpConfig.dataDirWeb, service.icon))
                continue
            if link:
                if external == "true":
                    print(
                        '<td><a target="_blank" title="%s gene %s" href="%s"><img width="50" height="50" alt="%s gene %s" src="%s/%s"></a></td>' % (
                            name, gene, link, name, gene, efpConfig.dataDirWeb, service.icon))
                else:
                    print(
                        '<td><a target="_blank" title="%s for gene %s" href="%s"><img width="50" height="50" alt="%s for gene %s" src="%s/%s"></a></td>' % (
                            name, gene, link, name, gene, efpConfig.dataDirWeb, service.icon))
            else:
                print(
                    '<td><img target="_blank" title="%s found for gene %s" width="50" height="50" alt="%s found for %s" src="%s/%s"></td>' % (
                        name, gene, name, gene, efpConfig.dataDirWeb, service.icon))
        print('</tr></table>')
        print('</td></tr>')

    # Tabular Navigation
    if nav_bar and len(max_dict.items()) != 0:
        colour_step = int(255 / len(max_dict.items()))
        next_colour = 255
        repeated_data = {}
        # This is a sorted list of tuples (source, max)
        data_source_max = sorted(max_dict.items(), key=operator.itemgetter(1))

        if useThreshold is None:
            useThreshold = ""

        # Determine font colour and hyperlink
        for i in range(len(data_source_max)):
            # Background color
            if (((255 * 299) + (next_colour * 587)) / 1000) > 160:
                font = 'black'
            else:
                font = 'white'

            # Hyperlink to load the new view
            # data_source_max = (source, max, color, font, link)
            href_link = '<a style="text-decoration:none;color:%s" href="efpWeb.cgi?dataSource=%s&mode=%s&primaryGene=%s&secondaryGene=%s&override=%s&threshold=%s&modeMask_low=%s&modeMask_stddev=%s">' % (
                font, data_source_max[i][0], mode, gene1.gene_id, gene2.gene_id, useThreshold, threshold, grey_low, grey_stddev)

            # Foreground color
            if next_colour < 17:
                data_source_max[i] = data_source_max[i] + ("0" + str(hex(next_colour)[-1:]), font, href_link)
            else:
                data_source_max[i] = data_source_max[i] + (str(hex(next_colour)[-2:]), font, href_link)

            # Decrease color by color step
            if next_colour > colour_step:
                next_colour -= colour_step
            else:
                next_colour = 0

        # Some eFPs, a source name is repeated like Root, Root II.
        # Make it Root, II
        try:
            for counter in range(len(data_source_max)):
                if data_source_max[counter][0] in efpConfig.repeatSources:
                    output = list(data_source_max[counter][:])
                    output[0] = output[0] + '_II'
                    output[4] = '<a style="text-decoration:none;color:%s" href="efpWeb.cgi?dataSource=%s&mode=%s&primaryGene=%s&secondaryGene=%s&override=%s&threshold=%s&modeMask_low=%s&modeMask_stddev=%s">' % (
                        output[3], output[0], mode, gene1.gene_id, gene2.gene_id, useThreshold, threshold, grey_low, grey_stddev)
                    repeated_data[data_source_max[counter][0]] = tuple(output)

            # Now add II right after I. In most cases, they are the same database.
            for item in efpConfig.repeatSources:
                for counter in range(len(data_source_max)):
                    if data_source_max[counter][0] == item:
                        data_source_max.insert(counter + 1, repeated_data[item])
                        break

            # Remove name and keep only II
            for item in efpConfig.repeatSources:
                for counter in range(len(data_source_max)):
                    if data_source_max[counter][0] == item + '_II':
                        new_item = list(data_source_max[counter])
                        new_item[0] = "II"
                        data_source_max[counter] = tuple(new_item)

        except AttributeError:
            pass

        current_source = 0
        while current_source < len(data_source_max):
            if data_source_max[current_source][0] == dataSource:
                break
            else:
                current_source = current_source + 1

        print('<tr><td><table style="width:100;white-space:nowrap;text-align:center;margin-left:auto;margin-right:auto">')

        # Append display name:
        # data_source_max = (source, max, color, font, link, name)
        for i in range(len(data_source_max)):
            try:
                if data_source_max[i][0] in efpConfig.shortNames.keys():
                    data_source_max[i] = list(data_source_max[i])
                    data_source_max[i].append(data_source_max[i][0].replace('_', ' '))
                    data_source_max[i][0] = efpConfig.shortNames[data_source_max[i][0]]
                    data_source_max[i] = tuple(data_source_max[i])
                else:
                    data_source_max[i] = list(data_source_max[i])
                    data_source_max[i][0] = data_source_max[i][0].replace('_', ' ')
                    data_source_max[i].append(data_source_max[i][0])
                    data_source_max[i] = tuple(data_source_max[i])
            except AttributeError:
                data_source_max[i] = list(data_source_max[i])
                data_source_max[i][0] = data_source_max[i][0].replace('_', ' ')
                data_source_max[i].append(data_source_max[i][0])
                data_source_max[i] = tuple(data_source_max[i])

        # Append (max) to last item
        max_key = list(data_source_max[-1])
        max_key[0] = max_key[0] + " (max)"
        data_source_max[-1] = tuple(max_key)

        # Print the current source row
        print('    <tr>')
        counter = 0
        while counter < len(data_source_max):
            if counter == current_source:
                print('        <td bgcolor="#%s%s%s" rowspan=2><table style="border:1px solid black"><tr><td><span title="%s: %s"><font color="%s" style="font-weight:bold">%s<b><u>%s</u></b></a></font></span></td></tr></table></td>' % (
                    data_source_max[counter][2], data_source_max[counter][2], data_source_max[counter][2], data_source_max[counter][5],
                    round(data_source_max[counter][1][0], 2), data_source_max[counter][3], data_source_max[counter][4], data_source_max[counter][0]))
            else:
                print('<td></td>')
            counter = counter + 1
        print('    </tr>')

        # print rest of selections
        print('    <tr>')
        counter = 0
        while counter < len(data_source_max):
            if counter == current_source:
                print('     ')
            elif data_source_max[counter][1]:
                # Asher 2018: Fails if one of the views doesn't have data.
                print('        <td bgcolor="#%s%s%s"><span title="%s: %s"><font color="%s">%s<u>%s</u></a></font></span></td>' % (
                    data_source_max[counter][2], data_source_max[counter][2], data_source_max[counter][2], data_source_max[counter][5],
                    round(data_source_max[counter][1][0], 2), data_source_max[counter][3], data_source_max[counter][4],
                    data_source_max[counter][0].replace('_', ' ')))
            counter = counter + 1
        print('    </tr>')
        print('</table></td></tr>')

    # print the image
    view_no = 1
    for view_name in views:
        print('<tr align="center"><td>')
        imgFile = img_filename[view_name]
        temp_imgPath = img_filename[view_name].split("/")
        last_element = temp_imgPath[-1]
        match = re.search(r"^efp", last_element)
        if match is not None:
            imgFile = efpConfig.OUTPUT_FILES + '/%s' % last_element
        print('  <img src="%s/%s" border="0" ' % (CONF["OUTPUT_FILES_WEB"], os.path.basename(imgFile)))
        if view_name in img_map:
            print('usemap="#imgmap_%s">' % view_name)
            print('%s' % img_map[view_name])
        else:
            print('>')
        print('</td></tr>')
        # Creates Button and Link to Page for Table of Expression Values
        print('<tr align="center"><td><br>')
        temp_tablePath = table_file[view_name][1].split("/")
        table_file_name = CONF["OUTPUT_FILES_WEB"] + '/%s' % (temp_tablePath[-1])
        print(
            '<input type="button" name="expressiontable" value="Click Here for Table of Expression Values" onclick="resizeIframe(\'ifr%d\', ifr%d);popup(\'table%d\', \'fadein\', \'center\', 0, 1)">&nbsp;&nbsp;' % (
                view_no, view_no, view_no))
        tableChart_name = table_file_name.replace(".html", ".png")
        print(
            '<input type="button" name="expressionchart" value="Click Here for Chart of Expression Values" onclick="popup(\'chart%d\', \'fadein\', \'center\', 0, 1)">' % (
                view_no))
        print('<script type="text/javascript">')
        popup_content = '<span style="color:#000000;"><b>For table download right click <a href="%s">here</a> and select "Save Link As ..."</b></span>' % table_file_name

        # Add Genes
        if primaryGene is not None and primaryGene != '':
            popup_content += '<br><br><span style="color:#000000;"><b>' + primaryGene + '</b></span><br>'
        if mode == "Compare" and secondaryGene is not None and secondaryGene != '':
            popup_content += '<span style="color:#000000;"><b>' + secondaryGene + '</b></span><br>'

        popup_content += '<div class="closewin_text">'
        popup_content += '<a href="" onclick="popdown(\\\'table%d\\\');return false;">' % (view_no)
        popup_content += '<span style="color:#000000">[Close]</span></a><br><br>'
        popup_content += '<a href="" onclick="switchPopup(\\\'table%d\\\', \\\'chart%d\\\');return false;">' % (
            view_no, view_no)
        popup_content += '<span style="color:#000000">[Switch to<br> Chart]</span></a></div>'
        popup_content += '<div class="chart"><iframe id="ifr%d" name="ifr%d" width=900 frameborder=0 src="%s">' % (
            view_no, view_no, table_file_name)
        popup_content += 'Your browser doesn\\\'t support iframes. Please use link above to open expression table</iframe></div>'
        popup_width = '1000'
        bg_color = '#FFFFFF'
        print("loadPopup(\'table%d\',\'%s\',\'%s\',%s);" % (view_no, popup_content, bg_color, popup_width))
        popup_content = '<div class="closewin_text">'
        popup_content += '<a href="" onclick="popdown(\\\'chart%d\\\');return false;">' % view_no
        popup_content += '<span style="color:#000000">[Close]</span></a><br><br>'
        popup_content += '<a href="" onclick="resizeIframe(\\\'ifr%d\\\', ifr%d);switchPopup(\\\'chart%d\\\', \\\'table%d\\\');return false;">' % (
            view_no, view_no, view_no, view_no)
        popup_content += '<span style="color:#000000">[Switch to<br>Table]</span></a><br><br>'
        popup_content += '<a href="" onclick="zoomElement(\\\'image%d\\\', 0.1);return false;">' % view_no
        popup_content += '<span style="color:#000000">[Zoom +]</span></a><br>'
        popup_content += '<a href="" onclick="zoomElement(\\\'image%d\\\', -0.1);return false;">' % view_no
        popup_content += '<span style="color:#000000">[Zoom -]</span></a><br>'
        popup_content += '<a href="" onclick="zoomElement(\\\'image%d\\\', 0);return false;">' % view_no
        popup_content += '<span style="color:#000000">[Reset<br>zoom]</span></a></div>'
        popup_content += '<div class="chart">'

        # Add Genes
        if primaryGene is not None and primaryGene != '':
            popup_content += '<br><br><span style="color:#000000;"><b>' + primaryGene + '</b></span><br>'
        if mode == "Compare" and secondaryGene is not None and secondaryGene != '':
            popup_content += '<span style="color:#000000;"><b>' + secondaryGene + '</b></span><br>'

        popup_content += '<img id="image%d" height="580px" src="%s"><br></div>' % (
            view_no, tableChart_name)
        print("loadPopup(\'chart%d\',\'%s\',\'%s\',%s);" % (view_no, popup_content, bg_color, popup_width))
        print("</script>")
        print('<br></td></tr>')
        view_no = view_no + 1
    print('  <tr><td><br><ul>')
    if efpConfig.species == "CAMELINA":
        print('  <li>%s is the Arabidopsis thaliana best hit for your Camelina sativa gene, %s </li>' %
              (gene1.get_gene_id(), gene1.get_probeset_id()))
    else:
        print('  <li>%s was used as the probe set identifier for your primary gene, %s (%s)</li>' %
              (gene1.get_probeset_id(), gene1.get_gene_id(), gene1.get_annotation()))
    if mode == 'Compare':
        if efpConfig.species == "CAMELINA":
            print('  <li>%s is the Arabidopsis thaliana best hit for your Camelina sativa gene, %s </li>' %
                  (gene2.get_gene_id(), gene2.get_probeset_id()))
        else:
            print('  <li>%s was used as the probe set identifier for the secondary gene, %s (%s)</li>' %
                  (gene2.get_probeset_id(), gene2.get_gene_id(), gene2.get_annotation()))
    print('  </ul>')
    if dataSource in efpConfig.datasourceFooter:
        print(efpConfig.datasourceFooter[dataSource])
    else:
        print(efpConfig.datasourceFooter['default'])
    print('</td></tr>')
    print('<tr><td></td></tr>')
else:
    print('<tr align="center"><td>')
    print('  <img src="%s" border="0">' % default_img_filename)
    print('</td></tr>')
print('</table>')
print('</form>')
print('</body>')
print('</html>')
