import shutil
import os
import sys
import cgi
import time

from iotbx.phil import parse
from cctbx.sgtbx import space_group_symbols

# run_job module from Fast_ep
from run_job import run_job, run_job_cluster, is_cluster_job_finished

import log_read
import html_report


class logger:
    def __init__(self):
        self._fout = open('adpssp.log', 'w')
        self.start_time = {}
        return

    def __del__(self):
        self._fout.close()
        self._cout = None
        return

    def __call__(self, _line):
        sys.stdout.write('%s\n' % _line)
        self._fout.write('%s\n' % _line)
        return

    def start(self, _program):
        self.__call__("### %s: Started ###" % _program)
        self.start_time[_program] = time.time()

    def done(self, _program):
        self.__call__("### %s: Done (Elapsed time: %.2f) ###" % (_program,
                time.time() - self.start_time[_program]))
        return

class Adpssp_parameters():
    def __init__(self):
        self._phil = parse("""
adpssp {
    data = ""
        .type = str
    sequence = ""
        .type = str
    wd = ""
        .type = str
    jobs = "fast_dp,fast_ep,shelxe_autotrace,dm,arpwarp"
        .type = str
    fast_dp {
        sw_trigger = True
            .type = str
    }
}
""")
        argument_interpreter = self._phil.command_line_argument_interpreter(
            home_scope = 'adpssp')
        for argv in sys.argv[1:]:
            command_line_phil = argument_interpreter.process(arg = argv)
            self._phil = self._phil.fetch(command_line_phil)
        self._parameters = self._phil.extract()

        return

    def get_data(self):
        return self._parameters.adpssp.data

    def get_sequence(self):
        return self._parameters.adpssp.sequence

    def get_wd(self):
        return self._parameters.adpssp.wd

    def get_jobs(self):
        return self._parameters.adpssp.jobs.split(',')

    def get_fast_dp_sw_trigger(self):
        return self._parameters.adpssp.fast_dp.sw_trigger


class Adpssp():
    def __init__(self, _parameters):
        self._log = logger()
        self.start_time = time.time()
        self._wd = _parameters.get_wd()
        if self._wd == "":
            self._wd = os.getcwd()

        if _parameters.get_data() == "":
            self._log("No data was input!")
            self._log("Abort")
            sys.exit()
        self._parameters = _parameters
        self._data = os.path.join(self._wd, _parameters.get_data())
        self._sequence = _parameters.get_sequence()

        self.job_cnt = 0

        self.fast_dp_stats = {}
        self.fast_ep_stats = {}
        self.shelxe_autotrace_stats = {}
        self.dm_stats = {}
        self.parrot_stats = {}
        self.arpwarp_stats = {}

        self.html = html_report.Html_report(self._wd, self._data)

        self.jobs = self.check_input()
        self.next_job = self.jobs[self.job_cnt]
        self.mainloop()

    def check_input(self):
        jobs = self._parameters.get_jobs()

        suffix = self._data.split(".")[-1]
        if suffix in ("cbf" ,"h5"):
            self._log("The data provided is a set of diffraction images.")
            self.fast_dp_filein = self._data
        elif suffix in ("mtz", "hkl", "HKL"):
            self._log("The data provided is a .mtz file.")
            self.jobs = self.jobs[self.jobs.index('fast_dp')+1:]
            jobs = jobs[job.index('fast_dp')+1:]
            self.fast_ep_filein = self._data
        else:
            self._log("The file %s does not have a proper suffix." % self._data)
            self._log("Abort")
            sys.exit()

        return jobs

    def fast_dp(self):
        job_name = 'fast_dp'
        self._log.start(job_name)
        wd = self.make_dir(job_name)

        fast_dp_option = ""

        if self._parameters.get_fast_dp_sw_trigger():
            fast_dp_option += " --sw-trigger"

        output = run_job('fast_dp', ['%s' % fast_dp_option, '%s' % self.fast_dp_filein], [], wd)

        logfile = os.path.join(os.path.join(wd, 'fast_dp.log'))
        self.fast_dp_stats = log_read.fast_dp(logfile)
        self.html.fast_dp(logfile)

        self._log("fast_dp has been completed")
        self._log("Point group: %s" % self.fast_dp_stats["point_group"])
        self._log("Unit cell:   %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f" % tuple(self.fast_dp_stats["cell"]))
        self._log("Max. reso:   %5.2f" % self.fast_dp_stats["high_reso"][0])

        self._log.done(job_name)

        if self.fast_dp_stats["sigano_low"] > 1.5:
            self._log("There are strong anomalous signal")
            self._log("Phasing with Fast_ep will be tried")
            self.fast_ep_filein = os.path.join(wd, 'fast_dp.mtz')
            self.next_job = "fast_ep"
        elif self.fast_dp_stats["sigano_low"] > 1.0:
            self._log("There might be a weak anomalous signal")
            self._log("Phasing with Fast_ep_weak will be tried")
            self.fast_ep_filein = os.path.join(wd, 'fast_dp.mtz')
            self.next_job = "fast_ep_weak"
        else:
            self._log("There seems to be no anomalous signal")
            self._log("No phasing trial will be done")
            self.next_job = "end"

        self.previous_job = job_name
        self.previous_wd = wd
        self.current_mtz = os.path.join(wd, 'fast_dp.mtz')

        self.job_cnt += 1


    def fast_ep(self, weak=False):
        if weak:
            job_name = 'fast_ep_weak'
        else:
            job_name = 'fast_ep'
        self._log.start(job_name)
        wd = self.make_dir(job_name)

        job_output = run_job('fast_ep', ['sad=%s' % self.fast_ep_filein, 'weak=%s' % weak],[],wd)

        logfile = os.path.join(wd, 'fast_ep.log')

        self.fast_ep_stats = log_read.fast_ep(logfile)

        sg = space_group_symbols(self.fast_ep_stats['spacegroup'])

        job_output = run_job('convert2mtz', ['-hklin', 'sad.phs', '-mtzout', 'sad.mtz',
                                             '-colin', 'F FOM PHI SIGF',
                                             '-cell', '%f %f %f %f %f %f' % tuple(self.fast_dp_stats['cell']),
                                             '-spacegroup', '%d' % sg.number()],
                             [], wd)

        self.current_mtz = os.path.join(wd, 'sad.mtz')
        self.current_hat = os.path.join(wd, 'sad.hat')

        self.html.fast_ep(logfile)
        self._log.done(job_name)

        self.previous_job = job_name
        self.previous_wd = wd
        self.job_cnt += 1
        self.next_job = self.jobs[self.job_cnt]


    def shelxe_autotrace(self):
        job_name = 'shelxe_autotrace'
        self._log.start(job_name)
        wd = self.make_dir(job_name)

        solv = self.fast_ep_stats['solvent']

        for f in ("sad_fa.hkl", "sad.hkl"):
            shutil.copy(os.path.join(self.previous_wd, f),
                        os.path.join(wd, f))

        shutil.copy(os.path.join(self.previous_wd, 'sad.hat'),
                    os.path.join(wd, 'sad_fa.res'))
        
        job_output = run_job('shelxe', ['sad', 'sad_fa', '-s%.1f' % solv, '-a2'], [], wd)
            
        logfile = os.path.join(wd, 'sad.lst')
        self.shelxe_autotrace_stats = log_read.shelxe_autotrace(logfile)

        self.html.shelxe_autotrace(os.path.join(wd, 'sad.lst'))

        sg = space_group_symbols(self.fast_ep_stats['spacegroup'])

        job_output = run_job('convert2mtz', ['-hklin', 'sad.phs', '-mtzout', 'sad.mtz',
                                             '-colin', 'F FOM PHI SIGF',
                                             '-cell', '%f %f %f %f %f %f' % tuple(self.fast_dp_stats['cell']),
                                             '-spacegroup', '%d' % sg.number()],
                             [], wd)

        self.current_mtz = os.path.join(wd, 'sad.mtz')
        self.current_hat = os.path.join(wd, 'sad.hat')

        self._log.done(job_name)

        self.previous_job = 'shelxe_autotrace'
        self.previous_wd = wd
        self.job_cnt += 1
        self.next_job = self.jobs[self.job_cnt]

    def dm(self):
        job_name = 'dm'
        self._log.start(job_name)
        wd = self.make_dir(job_name)

        mtzin = self.current_mtz

        job_output = run_job('dm',
                             ['hklin', mtzin, 'hklout', 'dm.mtz'],
                             ['SOLC %f' % self.fast_ep_stats['solvent'],
                              'MODE SOLV HIST MULTI',
                              'COMBINE PERT',
                              'NCYCLE 10',
                              'LABIN FP=F SIGFP=SIGF PHIO=PHI FOMO=FOM FREE=FreeF_flag',
                              'LABO PHIDM=PHIDM FOMDM=FOMDM'],
                              wd)

        logfile = os.path.join(wd, 'dm.log')
        open(logfile, 'w').write("".join(job_output))

        self.html.dm(logfile)

        self._log.done(job_name)

        self.current_mtz = os.path.join(wd, 'dm.mtz')

        self.previous_job = 'dm'
        self.previous_d = wd

        self.job_cnt += 1
        self.next_job = self.jobs[self.job_cnt]

    def parrot(self):
        job_name = 'parrot'

        if not 'CLIBD' in os.environ:
            raise RuntimeError, 'CLIBD not set'

        ref_dir = os.path.join(os.environ['CLIBD'], 'reference_structures')

        self._log.start(job_name)
        wd = self.make_dir(job_name)
        
        mtzin = self.current_mtz
        mtzout = os.path.join(wd, 'parrot.mtz')

        if self.previous_job == 'fast_ep':
            colin_wrk_phifom = "PHI,FOM"
        elif self.previous_job == 'shelxe_autotrace':
            colin_wrk_phifom = "PHI,FOM"
        elif self.previous_job == 'shelxe_autotrace':
            colin_wrk_phifom = "PHIDM,FOMDM"

        cparrot_sh = """cparrot \
-pdbin-ref %(pdbin_ref)s \
-mtzin-ref %(mtzin_ref)s \
-colin-ref-fo '/*/*/[FP.F_sigF.F,FP.F_sigF.sigF]' \
-colin-ref-hl '/*/*/[FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D]' \
-mtzin-wrk %(mtzin_wrk)s \
-colin-wrk-fo '/*/*/[F,SIGF]' \
-colin-wrk-phifom '/*/*/[%(colin_wrk_phifom)s]' \
-colin-wrk-free '/*/*/[FreeF_flag]' \
-mtzout %(mtzout)s \
-colout parrot \
-solvent-flatten \
-histogram-match \
-ncs-average \
-anisotropy-correction \
-cycles 3 \
-resolution 1.0 \
-ncs-mask-filter-radius 6.0""" % dict(pdbin_ref=os.path.join(ref_dir, 'reference-1tqw.pdb'),
                                      mtzin_ref=os.path.join(ref_dir, 'reference-1tqw.mtz'),
                                      colin_wrk_phifom=colin_wrk_phifom,
                                      mtzin_wrk=mtzin,
                                      mtzout=mtzout)
        print cparrot_sh

        open(os.path.join(wd, 'parrot.sh'), "w").write(cparrot_sh)

        job_output = run_job('sh', ['parrot.sh',], [], wd)

        open(os.path.join(wd, '%s.log' % job_name), "w").write("".join(job_output))

        self._log.done(job_name)
        self.previous_job = 'parrot'
        self.previous_wd = wd

        self.job_cnt += 1
        self.next_job = self.jobs[self.job_cnt]

    def arpwarp(self):

        if self._sequence == "":
            self._log('The sequence file is not specified.')
            self._log('Arpwarp job was cancelled.')
            self.next_job = 'end'
            return

        job_name = 'arpwarp'
        self._log(job_name)
        wd = self.make_dir(job_name)

        shutil.copy(os.path.join(self._wd, self._sequence),
                    os.path.join(wd, 'sequence.txt'))

        data_file = self.current_mtz

        job_output = run_job('auto_tracing.sh', ['datafile', data_file,
                                                 'fp', 'F', 'sigfp', 'SIGF',
                                                 'freelabin', 'FreeF-flag',
                                                 'phibest', 'PHIDM', 'fom', 'FOM',
                                                 'modelin', '../shelxe/sad.pdb',
                                                 'seqin', 'sequence.txt',
                                                 'buildingcycles', '2',
                                                 'restrref', '1',
                                                 'rrcyc', '1'],
                             [], wd)


        logfile = os.path.join(wd, 'auto_tracing.log')
        open(os.path.join(wd, 'auto_tracing.log'), "w").write("".join(job_output))
        
        self.html.arpwarp(logfile)

        self._log(job_name)

        self.next_job = "end"

    def make_dir(self, job_name):
        wd = os.path.join(self._wd, "%02d_%s" % (self.job_cnt+1, job_name))
        try:
            os.makedirs(wd)
        except:
            pass

        return wd

    def mainloop(self):
        while 1:
            if self.next_job == "fast_dp":
                self.fast_dp()

            elif self.next_job == "fast_ep":
                self.fast_ep()

            elif self.next_job == "fast_ep_weak":
                self.fast_ep(weak=True)

            elif self.next_job == "shelxe_autotrace":
                self.shelxe_autotrace()

            elif self.next_job == "dm":
                self.dm()

            elif self.next_job == "parrot":
                self.parrot()

            elif self.next_job == "arpwarp":
                self.arpwarp()

            elif self.next_job == "end":
                self._log("Adpssp has been completed.")
                self.end_time = time.time()
                self._log("Total elapsed time: %.2f" % (self.end_time - self.start_time))
                break

            else:
                print "Internal error!"
                print "The variable next_job has a wrong keyword: %s " % self.next_job
                sys.exit(-1)

if __name__ == '__main__':
    wd = "/GPFS/CENTRAL/XF17ID1/yamada/work/fast_ep_weak/html_test"
    data_in = "/GPFS/CENTRAL/XF17ID1/yamada/work/fast_ep_weak/test/fast_dp.mtz"

    fin = "/GPFS/CENTRAL/XF17ID2/jjakoncic/Dec_2016_DO/projID/XtalSamp_21_15/3/XtalSamp_21_15_00004.cbf"
    fin = "/GPFS/CENTRAL/XF17ID2/jjakoncic/Aug_15_2016/projID/XtalSamp_6_3/10/XtalSamp_6_3_00004.cbf"
    adpssp = Adpssp(Adpssp_parameters())

#    fast_ep = Fast_ep_nsls2(Fast_ep_parameters())
#    fast_ep.shelxe_autotrace()
#    fast_ep.dm()
#    fast_ep.parrot()
#    fast_ep.arpwarp()
