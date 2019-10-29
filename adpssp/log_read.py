import os

from iotbx.pdb import hierarchy


def fast_dp(logfile):
    fast_dp_stats = {}

    for line in open(logfile).readlines():
        if line.startswith("Unit cell"):
            fast_dp_stats["cell"] = [float(x) for x in line.split()[2:8]]
        elif line.startswith("Merging point group"):
            fast_dp_stats["point_group"] = "".join(line.split()[3:7])
        elif line.startswith("      Low resolution"):
            fast_dp_stats["low_reso"] = [float(x) for x in line.split()[2:5]]
        elif line.startswith("     High resolution"):
            fast_dp_stats["high_reso"] = [float(x) for x in line.split()[2:5]]
        elif line.startswith("              Rmerge"):
            fast_dp_stats["rmerge"] = [float(x) for x in line.split()[1:4]]
        elif line.startswith("             I/sigma"):
            fast_dp_stats["isigma"] = [float(x) for x in line.split()[1:4]]
        elif line.startswith("        Completeness"):
            fast_dp_stats["completeness"] = [float(x) for x in line.split()[1:4]]
        elif line.startswith("        Multiplicity"):
            fast_dp_stats["multiplicity"] = [float(x) for x in line.split()[1:4]]
        elif line.startswith("              CC 1/2"):
            fast_dp_stats["cchalf"] = [float(x) for x in line.split()[2:5]]
        elif line.startswith("   Anom. Correlation"):
            fast_dp_stats["ccano"] = [float(x) for x in line.split()[2:5]]
        elif line.startswith("Mosaic spread"):
            fast_dp_stats["mos_min"] = float(line.split()[2])
            fast_dp_stats["mos_ave"] = float(line.split()[4])
            fast_dp_stats["mos_max"] = float(line.split()[6])

    wd = os.path.dirname(logfile)
    correct_lp = open(os.path.join(wd,'CORRECT.LP')).readlines()

    for i in range(len(correct_lp)):
        if 'SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION' in correct_lp[i]:
            sigano_low = correct_lp[i+4].split()[12]

    fast_dp_stats['sigano_low'] = sigano_low

    return fast_dp_stats

def fast_ep(logfile):
    fast_ep_stats = {}

    for line in open(logfile).readlines():
        if line.startswith("Best spacegroup"):
            fast_ep_stats['spacegroup'] = line.split()[2]
        elif line.startswith("Best solvent"):
            fast_ep_stats['solvent'] = float(line.split()[2])
        elif line.startswith("Best nsite"):
            fast_ep_stats['nsite'] = float(line.split()[2])
        elif line.startswith("Best CC / weak"):
            fast_ep_stats['ccall'] = float(line.split()[4])
            fast_ep_stats['ccweak'] = float(line.split()[6])
            fast_ep_stats['cfom'] = fast_ep_stats['ccall'] + fast_ep_stats['ccweak']

    return fast_ep_stats

def shelxe_autotrace(logfile):
    shelxe_autotrace_stats = {}

    for line in open(logfile).readlines():
        if line.startswith(" Best trace"):
            shelxe_autotrace_stats['cc'] = float(line.split()[6][:-2])
            pdb_name = line.split()[-1]

    logdir = os.path.dirname(logfile)

    pdb_in = hierarchy.input(file_name=os.path.join(logdir, pdb_name))

    shelxe_autotrace_stats['nchain'] = len(pdb_in.hierarchy.only_model().chains())

    nresidue = 0
    for i in pdb_in.hierarchy.only_model().chains():
        nresidue += len(i.residue_groups())

    shelxe_autotrace_stats['nresidue'] = nresidue

    return shelxe_autotrace_stats

def dm(logfile):
    dm_stats = {}

    for line in open(logfile).readlines():
        pass

    return dm_stats

def arpwarp(logfile):
    arpwarp_stats = {}

    for line in open(logfile).readlines():
        if "refmac" in line:
            arpwarp_stats['rmerge'] = float(line.split()[6])
        elif "Sequence coverage" in line:
            arpwarp_stats['seq_coverage'] = float(line.split()[3])
        elif "The working directory:" in line:
            arpwarp_dir = line.split()[3]

    pdb_in = hierarchy.input(file_name=os.path.join(arpwarp_dir, 'dm_warpNtrace.pdb'))
    arpwarp_stats['nchain'] = len(pdb_in.hierarchy.only_model().chains())

    nresidue = 0
    for i in pdb_in.hierarchy.only_model().chains():
        nresidue += len(i.residue_groups())

    arpwarp_stats['nresidue'] = nresidue

    return arpwarp_stats
                    
