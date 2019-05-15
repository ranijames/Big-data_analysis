import logging
import sys
from allpy.exon.mutations.utils import Timer
from allpy.pipelineutils.filepaths import WorkDirs, OutFiles

__author__ = 'mpschr'


class SimpleTask():

    #logging:
    logging.basicConfig(stream=sys.stdout, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=logging.DEBUG)

    def __init__(self, task_name):
        self.task_name = task_name
        self.timer = None

    def run(self, **kwargs):
        logging.info('Running {}\n================================================='.format(self.task_name))
        self.timer = Timer(start_now=True)


class SampleTask(SimpleTask):

    def __init__(self, task_name):
        super().__init__(task_name)
        self.work_dirs = None
        self.output_files = None
        self.sample = None

    def run(self, task_input, exec_dir, cmd_sample_name=None, patient_id=None):
        super().run()

        self.sample = self.derive_sample_name(task_input, cmd_sample_name)
        self.patient = patient_id
        if len(self.sample) == 0:
            raise RuntimeError("Not valid sample name: {}".format(self.sample))
        else:
            logging.info("Sample: {}{}".format(self.sample, " / Patient: "+ self.patient if self.patient is not None else ""))

        self.work_dirs = WorkDirs(exec_dir, self.sample, patient=self.patient)
        self.output_files = OutFiles(work_dirs=self.work_dirs)

    def derive_sample_name(self, task_input, cmd_sample_name):
        raise NotImplementedError()


class PatientTask(SimpleTask):
    def __init__(self, task_name):
        super().__init__(task_name)
        self.patient = None

    def run(self, task_input, exec_dir, cmd_patient_name):
        super().run()

        self.patient = self.derive_patient_name(task_input, cmd_patient_name)
        if len(self.patient) == 0:
            raise RuntimeError("Not valid sample name: {}".format(self.patient))
        else:
            logging.info("Sample: {}".format(self.patient))

    def derive_patient_name(self, task_input, cmd_patient_name):
        raise NotImplementedError()