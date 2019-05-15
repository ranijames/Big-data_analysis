import logging
import os

__author__ = 'ARJ'

import configobj

class Config(object):

    def __init__(self):
        self._config_object = self.load()
        self.initialized = True



    def load(self):
        try:
            conffile = os.path.expanduser('~/.agbaldus/agbaldus-seq.conf')
            if not os.path.exists(conffile):
                raise RuntimeError('Not existing confile: {}'.format(conffile))
            logging.debug("reading conf {}".format(conffile))
            return configobj.ConfigObj(conffile)
            #return parser.read('../config/pipelines.conf.template')
        except OSError:
            raise RuntimeError("Please create the configuration file at ~/.agbaldus/agbaldus-seq.conf")

    def get_config(self, keys):
        if type(keys) is str:
            keys = keys.split('.')

        c = self._config_object
        for k in keys:
            if k in c:
                c = c[k]
            else:
                raise RuntimeError("ERROR retrieving config from {}".format('.'.join(keys)))
        return c

    def get_config_dict(self, key_list):
        c = {}
        for k in key_list:
            c[k] = self.get_config(k)
            logging.debug("CONFIG loaded: {} = {}".format(k, c[k]))
        return c
