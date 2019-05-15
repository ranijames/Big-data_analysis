import os

ALIGN = 'align'
CUFFLINKS = 'cufflinks'
CUFFLINKS_NOVEL = 'cufflinks_novel'
TRIM = 'trim'
SCRATCH_DIR = 'scratch'
TASKOUT_DIR = 'task_outputs'
QUALIMAP_EXON = 'qualimap_exon'
QUALIMAP_RNA = 'qualimap_rna'

__author__ = 'ARJ'


class WorkDirs():
    def __init__(self, exec_dir=None, sample=None, patient=None):
        self.exec_dir = exec_dir
        self.sample = sample
        self.patient = patient
        self.sample_subdir = None

        if patient is not None:
            self.sample_subdir = os.path.join(patient,sample)
        else:
            self.sample_subdir = self.sample

    def get(self, key, scratch=False):

        if scratch:
            scratch_or_final = SCRATCH_DIR
        else:
            scratch_or_final = TASKOUT_DIR

        if key == ALIGN:
            dir = os.path.join(self.exec_dir, scratch_or_final, self.sample_subdir, ALIGN)

        elif key == QUALIMAP_RNA:
            dir = os.path.join(self.exec_dir, scratch_or_final, self.sample_subdir, ALIGN, 'qc')

        else:
            dir = os.path.join(self.exec_dir, scratch_or_final, self.sample_subdir, key)

        if not os.path.exists(dir):
            os.makedirs(dir)

        return dir


class OutFiles():
    def __init__(self, work_dirs=None):
        self.exec_dir = work_dirs.exec_dir
        self.sample = work_dirs.sample
        self.patient = work_dirs.patient
        self.wd = work_dirs

    def generate(self, key, outfile):
        return os.path.join(self.wd.get(key), outfile)

    def get(self, key):
        if key == ALIGN:
            if self.sample is None:
                return "{}/{}.bam".format(self.wd.get(ALIGN), "*")
            else:
                return "{}/{}.bam".format(self.wd.get(ALIGN), self.sample)
        else:
            return None
