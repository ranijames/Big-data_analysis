import argparse
from logging import basicConfig, DEBUG, info, warning
import os
import itertools
import yaml
from copy import deepcopy
from socket import gethostname
import datetime

#logging:
import sys

basicConfig(stream=sys.stdout, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=DEBUG)


def listdir_abspath(d):
    d = os.path.abspath(d)
    files = []
    for sampledir in os.listdir(d):
        for fastq in os.listdir(os.path.join(d,sampledir)):
            files.append( os.path.join(d,sampledir,fastq))
            #files = files + [os.path.join(sampledir, f) for f in os.listdir(sampledir)]
    return files


class ConfigCreator():

    def __init__(self, yaml_file, patient_data_folder, config_dir, test_run=False, force_write=False, time_points=[]):
        self.force_write = True
        self.config_dir = config_dir
        self.test_run = test_run
        self.patient_data_folder = patient_data_folder
        self.yaml_file = yaml_file
        self.required_time_points = time_points
        print(self.config_dir)

    def run(self):

        info('CREATING CONFIG yaml\n=========================')

        sample_files = listdir_abspath(self.patient_data_folder)
        sample_files.sort()
        print(sample_files)

        patient = os.path.basename(self.patient_data_folder)

        time_points = {}
        for k, v in itertools.groupby(sample_files, key=lambda x: x.split('_')[-2]):
            print('---> {}'.format(k))
            if len(self.required_time_points) == 0 or k in self.required_time_points:
                time_points[k] = list(v)
                info([k,"-->", ":".join(time_points[k])])

        with open(self.yaml_file, 'r') as f:

            template = yaml.load(f)
            print(template)
            details = deepcopy(template['details'][0])
            #algorithm = deepcopy(template['details']['algorithm'])
            template['details'] = []

            for t in time_points.keys():
                this_details = deepcopy(details)
                this_details['description'] = patient + '_' + t
                this_details['files'] = time_points[t]
                metadata = {}
                metadata['phenotype'] = 'tumor' if t != 'CR' else 'normal'
                metadata['batch'] = self.get_batches(t, time_points.keys(), patient)
                this_details['metadata'] = metadata
                template['details'].append(this_details)

            template['fc_date'] = patient
            template['fc_name'] = gethostname() + str(datetime.date.today()).replace('-','')
            template['upload']['dir'] = os.path.join(template['upload']['dir'],patient,'bcbio')

            out_yaml = os.path.join(self.config_dir, patient + '.yaml')
            if not os.path.exists(out_yaml):
                yaml.dump(template, open(out_yaml, 'w'), Dumper=yaml.SafeDumper)
                info('Config file written for {}: {}'.format(patient, out_yaml))
            elif self.force_write:
                os.rename(out_yaml, out_yaml + '.bkp-' + str(datetime.datetime.now()).split('.')[0].replace('-','').replace(' ', '-').replace(':',''))
                yaml.dump(template, open(out_yaml, 'w'), Dumper=yaml.SafeDumper)
                warning('Old config file backed up {}'.format(patient))
                info('Config file written for {}: {}'.format(patient, out_yaml))
            else:
                raise Exception("YAML file already exists: {}".format(out_yaml))

    def get_batches(self, this_tt, all_time_points, patient):
        if this_tt == 'ID' or this_tt == 'REL' or this_tt == 'REL2':
            return 'b' + patient + '_' + this_tt
        elif this_tt == 'CR':
            return ['b' + patient + '_' + t for t in all_time_points if not t == 'CR']
        else:
            raise Exception('Wrong time point for sample: {}'.format(this_tt))



def cmdline():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='patient_data_folder', required=True)
    parser.add_argument('-y', dest='yaml_file', required=True)
    parser.add_argument('-t', '--test-run', action='store_true', dest='test',  default=False)
    parser.add_argument('-c', '--config-dir', dest='config_dir')
    options = parser.parse_args()
    task = ConfigCreator(yaml_file=os.path.abspath(options.yaml_file),
                         patient_data_folder=os.path.abspath(options.patient_data_folder),
                         test_run=options.test,
                         config_dir=os.path.abspath(options.config_dir))
    task.run()


if __name__ == "__main__":  # detects if called from system
    cmdline()
