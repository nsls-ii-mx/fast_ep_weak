import os
import sys

def make_html(_wd, _best_results, target='.'):

    spacegroups = []
    nsites = []
    resos = []

    for (spacegroup, nsite, reso) in _best_results.keys():
        if not spacegroup in spacegroups:
            spacegroups.append(spacegroup)
        if not nsite in nsites:
            nsites.append(nsite)
        if not reso in resos:
            resos.append(reso)

    nsites.sort()
    resos.sort()

    output = '<!DOCTYPE html>\n'
    output += '<html>\n'
    output += '  <head>\n'
    output += '    <title>FAST_EP_WEAK RESULT</title>\n'
    output += '  </head>\n'
    output += '  <body>\n'
    for spacegroup in spacegroups:
        output += '    <p>SPACEGROUP: %s</p>\n' % spacegroup
        output += '    <table border="1">\n'
        output += '      <tr>\n'
        output += '        <td> </td>\n'

        for nsite in nsites:
            output += '        <td align="center">%d</td>\n' % nsite

        output += '      </tr>\n'

        for reso in resos:
            output += '      <tr>\n'
            output += '        <td>%.2f</td>\n' % reso

            for nsite in nsites:
                output += '        <td valign="top">\n'
                output += '          <img src="%s" width="300px"><br>\n' % (
                        os.path.join(spacegroup, str(nsite), str(reso), 'sad_fa.png'))
                if 'shelxe_fom' in _best_results[(spacegroup, nsite, reso)]:
                    shelxe_result(_wd, spacegroup, nsite, reso)
                    output += '          <a href="%s">SHELXE</a><br>\n' % (
                            os.path.join(spacegroup, str(nsite), str(reso), 'shelxe.html'))
                output += '        </td>\n'

            output += '      </tr>\n'

        output += '    </table>\n'
    output += '  </body>\n'
    output += '</html>\n'

    open(os.path.join(_wd, 'fast_ep_result.html'), 'w').write(output)

def shelxe_result(_wd, spacegroup, nsite, reso):
    shelxe_dir = os.path.join(_wd, spacegroup, str(nsite), str(reso))
    solvs = []

    output = """\
<!DOCTYPE html>
<html>
  <head>
    <title>SHELXE RESULT</title>
  </head>
  <body>
    <center>
    <table border="1">
      <tr>
        <td>Working directory</td><td>%(wd)s</td>
      </tr>
      <tr>
        <td>Spacegroup</td><td>%(sgname)s</td>
      </tr>
      <tr>
        <td>Number of site</td><td>%(nsite)d</td>
      </tr>
      <tr>
        <td>Resolution</td><td>%(reso)f</td>
      </tr>
    </table>

    <table border="1">
""" % dict(wd=_wd, sgname=spacegroup, nsite=nsite, reso=reso)

    for d in os.listdir(shelxe_dir):
        if os.path.isfile(os.path.join(shelxe_dir, d, 'sad.png')):
            solvs.append(d)

    for s in sorted(solvs):
        output +="""\
      <tr>
        <td>%(solv)s</td><td><img src="%(image)s" width="300px"></td>
      </tr>
""" % dict(solv=s, image=os.path.join(s,'sad.png'))

    output +="""\
    </table>
    </center>
  </body>
</html>
"""

    open(os.path.join(shelxe_dir, 'shelxe.html'), "w").write(output)

if __name__ == '__main__':
    dirbase = sys.argv[1]
    print dirbase
    best_conditions = {('P212121', 28, 2.6): {'shelxe_diff': (0.25, 0.11399999999999993), 'shelxd': (51.86, 15.26, 67.12, 36), 'shelxe_fom': (0.4, 0.584, 'original')}, ('P212121', 32, 2.6): {'shelxe_diff': (0.4, 0.10399999999999998), 'shelxd': (53.75, 16.88, 70.63, 41), 'shelxe_fom': (0.45, 0.574, 'original')}, ('P212121', 28, 2.5): {'shelxe_diff': (0.5, 0.12599999999999995), 'shelxd': (52.2, 16.37, 68.58, 38), 'shelxe_fom': (0.4, 0.63, 'inverted')}, ('P212121', 26, 2.6): {'shelxe_diff': (0.7, 0.12299999999999997), 'shelxd': (52.23, 13.15, 65.39, 36), 'shelxe_fom': (0.4, 0.57, 'original')}, ('P212121', 30, 2.4): {'shelxe_diff': (0.65, 0.16600000000000004), 'shelxd': (51.53, 18.11, 69.64, 40), 'shelxe_fom': (0.4, 0.62, 'inverted')}, ('P212121', 32, 2.5): {'shelxe_diff': (0.75, 0.14), 'shelxd': (55.97, 18.64, 74.61, 39), 'shelxe_fom': (0.35, 0.632, 'original')}, ('P212121', 30, 2.6): {'shelxe_diff': (0.7, 0.238), 'shelxd': (53.01, 17.02, 70.04, 42), 'shelxe_fom': (0.4, 0.567, 'original')}, ('P212121', 26, 2.4): {'shelxe_diff': (0.3, 0.10599999999999998), 'shelxd': (52.63, 14.18, 66.8, 33), 'shelxe_fom': (0.45, 0.621, 'original')}, ('P212121', 32, 2.4): {'shelxe_diff': (0.35, 0.07299999999999995), 'shelxd': (54.04, 19.44, 73.48, 41), 'shelxe_fom': (0.35, 0.602, 'inverted')}, ('P212121', 26, 2.5): {'shelxe_diff': (0.6000000000000001, 0.17000000000000004), 'shelxd': (51.6, 15.72, 67.32, 31), 'shelxe_fom': (0.35, 0.563, 'original')}, ('P212121', 28, 2.4): {'shelxe_diff': (0.65, 0.22499999999999998), 'shelxd': (53.25, 17.02, 70.27, 38), 'shelxe_fom': (0.45, 0.585, 'original')}, ('P212121', 30, 2.5): {'shelxe_diff': (0.45, 0.14500000000000007), 'shelxd': (53.59, 18.04, 71.63, 38), 'shelxe_fom': (0.5, 0.556, 'original')}}
    make_html(dirbase, best_conditions)

            
