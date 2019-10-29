import shutil
import os
import sys
import cgi

import log_read

# escape_html from http://stackoverflow.com/questions/1061697/whats-the-easiest-way-to-escape-html-in-python
def escape_html(text):
    """escape strings for display in HTML"""
    return cgi.escape(text, quote=True).\
           replace(u'\n', u'<br />').\
           replace(u'\t', u'&emsp;').\
           replace(u'  ', u' &nbsp;')
        
class Html_report():
    def __init__(self, wd, data_in):
        self._wd = wd
        self._data_in = data_in

        self.filename = "adpssp_report.html"
        self.html = ""

        try:
            os.makedirs(self._wd)
        except:
            pass

        self.header()

    def header(self):
        self.html = """\
<html>
  <title>ADPSSP RESULTS</title>
  <head>
    <style type=text/css>
      <!--
        .monospace  {font-family:monospace;}
        table.stats {width: 800px;
                     border-collapse: separate;
                     border-spacing: 0px;
                     border-top: 1px solid #ccc;
                     border-left: 1px solid #ccc;
        }
        table.stats th {width: 45%%;
                        padding: 4px;
                        text-align: left;
                        vertical-align: top;
                        color: #444;
                        background-color: #feedf3;
                        border-left: 3px double #999;
                        border-top: 1px solid #fff;
                        border-right: 1px solid #ccc;
                        border-bottom: 1px solid #ccc;
        }
        table.stats td {width:55%%;
                        padding:4px;
                        background-color: #fafafa;
                        border-right: 1px solid #ccc;
                        border-bottom: 1px solid #ccc;
        }
        div.waku  {width: 800px;
                   border: 1px solid #badefe; 
                   padding:3px; 
                   margin:3px; 
                   display:none;
        }
      -->
    </style>
    <script language="JavaScript">
      <!--
        function log_open(geID, geID1, geID2){
            document.getElementById(geID).style.display = "block";
            document.getElementById(geID1).style.display = "none";
            document.getElementById(geID2).style.display = "inline";
        }
        function log_close(geID, geID1, geID2){
            document.getElementById(geID).style.display = "none";
            document.getElementById(geID1).style.display = "inline";
            document.getElementById(geID2).style.display = "none";
        }
      -->
    </script>
  </head>
  <body>
    <H1>NSLS-II ADPSSP RESULTS</H1>
    ADPSSP: Automatic Data Processing and Structure Solution Pipeline
    <p>Input file: %(data)s</p>
""" % dict(data=self._data_in)

        self.write_html()

    def write_html(self):
        text = self.html + """\
  </body>
</html>
"""
        open(os.path.join(self._wd, self.filename), "w").write(text)

    def diff_exp(self, wd):
        pass
        

    def fast_dp(self, logfile):
        self.fast_dp_stats = log_read.fast_dp(logfile)

        self.html += """
    <h2>Fast_dp</h2>
    <table class="stats">
      <tbody>
      <tr>
        <th>Point group</th>
        <td>%(point_group)s</td>
      </tr>
      <tr>
        <th>Unit cell</th>
        <td>%(unit_cell)s</td>
      </tr>
      <tr>
        <th>Resolution range</th>
        <td>%(low_reso)7.2f - %(high_reso)7.2f</td>
      </tr>
      <tr>
        <th>Overall Rmerge</th>
        <td>%(rmerge)6.3f</td>
      </tr>
      </tbody>
    </table>
    <a href="JavaScript:log_open('fast_dp_log','fast_dp_log1','fast_dp_log2');" id="fast_dp_log1">See log</a>
    <a href="JavaScript:log_close('fast_dp_log','fast_dp_log1','fast_dp_log2');" id="fast_dp_log2" style="display:none">Close</a>
    </div>
    <div id="fast_dp_log" class="waku">
      <span class="monospace">%(log)s</span>
    </div>
""" % dict(point_group=self.fast_dp_stats['point_group'],
           unit_cell=" ".join([str(x) for x in self.fast_dp_stats['cell']]),
           low_reso=self.fast_dp_stats['low_reso'][0],
           high_reso=self.fast_dp_stats['high_reso'][0],
           rmerge=self.fast_dp_stats['rmerge'][0],
           log=escape_html(open(logfile).read()))

        self.write_html()
        
    def fast_ep(self, logfile):
        self.fast_ep_stats = log_read.fast_ep(logfile)

        self.html += """
    <h2>Fast_ep</h2>
    <table class="stats">
      <tr>
        <th>Space group</hd>
        <td>%(spacegroup)s</td>
      </tr>
      <tr>
        <th>Number of sites</th>
        <td>%(nsite)s</td>
      </tr>
      <tr>
        <th>Best CFOM in SHELXD</th>
        <td>%(cfom)s</td>
      </tr>
      <tr>
        <th>Best solvent fraction in SHELXE</th>
        <td>%(solvent)4.2f</td>
      </tr>

    </table>
    <a href="JavaScript:log_open('fast_ep_log','fast_ep_log1','fast_ep_log2');" id="fast_ep_log1">See log</a>
    <a href="JavaScript:log_close('fast_ep_log','fast_ep_log1','fast_ep_log2');" id="fast_ep_log2" style="display:none">Close</a>
    </div>
    <div id="fast_ep_log" class="waku">
      <span class="monospace">%(log)s</span>
    </div>
""" % dict(spacegroup=self.fast_ep_stats['spacegroup'],
           nsite=self.fast_ep_stats['nsite'],
           cfom=self.fast_ep_stats['cfom'],
           solvent=self.fast_ep_stats['solvent'],
           log=escape_html(open(logfile).read()))

        self.write_html()
        
    def shelxe_autotrace(self, logfile):

        self.shelxe_autotrace_stats = log_read.shelxe_autotrace(logfile)

        self.html += """
    <h2>Shelxe auto tracing</h2>
    <table class="stats">
      <tr>
        <th>Number of chains built</th>
        <td>%(nchain)d</td>
      </tr>
      <tr>
        <th>Number of residues</th>
        <td>%(nresidue)d</td>
      </tr>
      <tr>
        <th>Correlation coefficient</th>
        <td>%(cc)6.2f</td>
      </tr>
    </table>
    <a href="JavaScript:log_open('shelxe_log','shelxe_log1','shelxe_log2');" id="shelxe_log1">See log</a>
    <a href="JavaScript:log_close('shelxe_log','shelxe_log1','shelxe_log2');" id="shelxe_log2" style="display:none">Close</a>
    </div>
    <div id="shelxe_log" class="waku">
      <span class="monospace">%(log)s</span>
    </div>
""" % dict(nchain=self.shelxe_autotrace_stats['nchain'],
           nresidue=self.shelxe_autotrace_stats['nresidue'],
           cc=self.shelxe_autotrace_stats['cc'],
           log=escape_html(open(logfile).read()))

        self.write_html()

    def dm(self, logfile):
        self.dm_stats = log_read.dm(logfile)

        self.html += """
    <h2>DM</h2>
    <a href="JavaScript:log_open('dm_log','dm_log1','dm_log2');" id="dm_log1">See log</a>
    <a href="JavaScript:log_close('dm_log','dm_log1','dm_log2');" id="dm_log2" style="display:none">Close</a>
    </div>
    <div id="dm_log" class="waku">
      <span class="monospace">%(log)s</span>
    </div>
""" % dict(log=escape_html(open(logfile).read()))

        self.write_html()
    # end of def dm

    def arpwarp(self, logfile):
        self.arpwarp_stats = log_read.arpwarp(logfile)

        self.html += """
    <h2>Arp/Warp</h2>
    <table class="stats">
      <tr>
        <th>Number of chains built</th>
        <td>%(nchain)d</td>
      </tr>
      <tr>
        <th>Number of residues</th>
        <td>%(nresidue)d</td>
      </tr>
      <tr>
        <th>Rmerge</th>
        <td>%(rmerge)6.2f</td>
      </tr>
    </table>
    <a href="JavaScript:log_open('dm_log','dm_log1','dm_log2');" id="dm_log1">See log</a>
    <a href="JavaScript:log_close('dm_log','dm_log1','dm_log2');" id="dm_log2" style="display:none">Close</a>
    </div>
    <div id="dm_log" class="waku">
      <span class="monospace">%(log)s</span>
    </div>
""" % dict(nchain=self.arpwarp_stats['nchain'],
           nresidue=self.arpwarp_stats['nresidue'],
           rmerge=self.arpwarp_stats['rmerge'],
           log=escape_html(open(logfile).read()))

        self.write_html()
    #end of def arpwarp
        
# end of class Html_report

