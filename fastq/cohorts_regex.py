__author__ = 'arj'


class FastqRegexCohort(object):
    """
    The patterns defined in this class serve to parse the fastq files and extract the
    strings for sample (patient and time_point) information.
    """
    def __init__(self):
        self.cohort_name = "OVERWRITE"
        self._SAMPLE_FASTQ_PATTERN = "OVERWRITE"
        self._TIMEPOINT_EXTRACTOR = "OVERWRITE"
        self._PATIENT_EXTRACTOR = "OVERWRITE"

    def get_fastq_pattern(self):
        fastqpattern = self._SAMPLE_FASTQ_PATTERN
        return self.fail_or_return(fastqpattern)


    def get_time_point_extractor(self):
        time_point = self._TIMEPOINT_EXTRACTOR
        return self.fail_or_return(time_point)

    def get_patient_extractor(self):
        patient = self._PATIENT_EXTRACTOR
        return self.fail_or_return(patient)

    def fail_or_return(self, var):
        if var == "OVERWRITE":
            raise RuntimeError("Not specified value in RegexCohort")
        return var

class VerbundsprojektDNASample(FastqRegexCohort):

    def __init__(self):
        super().__init__()
        self.cohort_name = "verbund_exon"
        self.PATIENT_PATTERN = "(PL|PE|AL|AE)[0-9]{2}"
        self.TP_PATTERN = "(REL|ID|CR)"

        self._SAMPLE_FASTQ_PATTERN = "^{}.+_L.+.fastq.gz$".format(self.PATIENT_PATTERN)
        self._TIMEPOINT_EXTRACTOR = '[0-9][0-9]_({}).+'.format(self.TP_PATTERN)
        self._PATIENT_EXTRACTOR = '({}).+'.format(self.PATIENT_PATTERN)


class VerbundsprojektRNASample(FastqRegexCohort):

    def __init__(self):
        super().__init__()
        self.cohort_name = "verbund_rna"

        self.PATIENT_PATTERN = "(PL|PE|AL|AE)[0-9]{2}"
        self.TP_PATTERN = "(rel|id)"
        self._SAMPLE_FASTQ_PATTERN = "^{}.+_L.+.fastq.gz$".format(self.PATIENT_PATTERN)
        self._TIMEPOINT_EXTRACTOR = '[0-9][0-9]({}).+'.format(self.TP_PATTERN)
        self._PATIENT_EXTRACTOR = '({}).+'.format(self.PATIENT_PATTERN)


FASTQ_COHORTS_REGEX = {
                        'verbund_rna': VerbundsprojektRNASample(),
                        'verbund_exon': VerbundsprojektDNASample()
                       }
