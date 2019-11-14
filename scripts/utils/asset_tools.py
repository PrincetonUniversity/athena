#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 ,-*
(_) Created on <Mo Okt 21 2019> @ 17:24:40

@author: Boris Daszuta
@function: Tools for controlling Athena++ and interacting with data
"""
###############################################################################
# BD: TODO
# this needs refactoring...
# should only deal with json and replacements.
# compile flags should also be converted appropriately
###############################################################################
# library and pathing
import os, sys
import re as _r
import collections as _c
import itertools as _it
import functools as _f
import operator as _op
import hashlib as _hl
import pathlib as _pl
import subprocess as _sub
import json as _j
import pickle as _p


# infer directory structure for Athena++ etc
_FN_SCRIPT = os.path.realpath(sys.argv[0])
_DIR_SCRIPTS = os.path.dirname(_FN_SCRIPT)
_DIR_ATHENA = os.path.split(_DIR_SCRIPTS)[0]
_DIR_REPO = os.path.split(_DIR_ATHENA)[0]
_DIR_PROBLEM_INP = os.path.join(_DIR_SCRIPTS, "problems")
_FN_SCRIPT_BASENAME = _FN_SCRIPT[len(_DIR_SCRIPTS)+1:-3]

sys.path.append(os.path.join(_DIR_ATHENA, 'vis', 'python'))

import athena_read as _ar   # for hdf5

###############################################################################
# constant strings
PROBLEM_TEMPLATE = '''
    <job>
    problem_id      = {problem_id}
    comment         = {comment}
    compile_flags   = {compile_flags}
    compile_hash    = {compile_hash}
    base_name       = {base_name}

    <time>
    cfl_number      = {cfl_number}
    tlim            = {tlim}
    integrator      = {integrator}

    <mesh>
    num_threads     = {num_threads}

    refinement      = {refinement}

    nx1             = {nx1}
    x1min           = {x1min}
    x1max           = {x1max}
    ix1_bc          = {ix1_bc}
    ox1_bc          = {ox1_bc}

    nx2             = {nx2}
    x2min           = {x2min}
    x2max           = {x2max}
    ix2_bc          = {ix2_bc}
    ox2_bc          = {ox2_bc}

    nx3             = {nx3}
    x3min           = {x3min}
    x3max           = {x3max}
    ix3_bc          = {ix3_bc}
    ox3_bc          = {ox3_bc}

    <meshblock>
    nx1             = {nx1}
    nx2             = {nx2}
    nx3             = {nx3}
'''

###############################################################################
# functions
def _sanitize_script_like(spec):
    # clean up an input string describing an Athena++ problem

    # strip comments
    tmp = _r.sub("#(.+?)\n", "\n", spec + "\n")
    # strip extraneous spaces
    tmp = _r.sub(" {2,}", " ", tmp)
    # pack new-lines
    tmp = tmp.replace(" \n", "\n").replace("\n ", "\n")
    # strip multiple new-lines
    tmp = _r.sub("\n{2,}", "\n", tmp)
    # get rid of useless new-lines
    tmp = tmp.split("\n")
    tmp = [el for el in tmp if (len(el) > 0) and (el != " ")]
    return "\n".join(tmp)

def _extract_from_dodn(inp, key_seq):
    # attempt extraction of value from arbitrarily nested dictionary
    try:
        # treat as final element
        return inp[key_seq]
    except:
        pass

    try:
        return _extract_from_dodn(inp[key_seq[0]], key_seq[1:])
    except:
        return inp[key_seq[0]]

def _dump(obj, filename):
    with open(os.path.join(_DIR_REPO, filename), 'wb') as f:
        _p.dump(obj, f)

###############################################################################
# classes
class _format_dict_filler(dict):
    def __missing__(self, key):
        return key.join("{}")

class problem_specification(object):
    def __init__(self, *args, template=None, **kwargs):
        # template parsing
        self._template_sections = None
        self._hash = None

        if template is not None:
            self.template = template
        else:
            self.template = PROBLEM_TEMPLATE


        self._specification_sections = {}

    @property
    def template(self):
        return self._template

    @template.setter
    def template(self, val):
        # sanitize
        self._template = _sanitize_script_like(val)

        # parse to sections
        s = self._template_section_extract
        l = list(self._template_section_labels) # conv. req. below

        # store as sorted and ordered
        ts = _c.OrderedDict()
        for k, v in sorted(zip(l, s)):
            ts[k] = v
        self._template_sections = ts

        # if sorting required regenerate
        if l != sorted(l):
            self._template = self._sections_to_template(
                populate_template=False)

        enc = self._template.encode("utf-8")
        self._hash = _hl.sha256(enc).hexdigest()

    @property
    def template_as_dod(self):
        # convert to dict of dicts
        od_dod = _c.OrderedDict()
        for sec, settings in self._template_sections.items():
            if len(settings) > 0:
                tmp = [setting.split(' = ', maxsplit=1) for setting in
                       settings.split('\n')]

                d_tmp = _c.OrderedDict()
                for el in tmp:
                    d_tmp[el[0]] = el[1]
                od_dod[sec] = d_tmp
            else:
                od_dod[sec] = None

        return od_dod

    @property
    def sections(self):
        return tuple(self._template_sections.keys())

    def get_section(self, section=None):
        try:
            return self._template_sections[section]
        except:
            raise AttributeError(
                "Section {section} ".format(section=section) +
                "is unknown.")

    @property
    def _template_section_labels(self):
        return _r.findall('<(.+?)>', self.template)

    @property
    def _template_section_extract(self):
        t = self.template.split("<")
        # strip all elements not containing ">"
        t = [_sanitize_script_like(el[el.find(">")+2:])
             for el in t if ">" in el]
        return t

    def _sections_to_template(self, populate_template=True):
        # convert dict repr to template string

        lbls = self._template_sections.keys()
        tstr = ""
        for l in lbls:
            tstr += "\n<" + l + ">\n" + self._template_sections[l]

        if populate_template:
            self.template = tstr
        else:
            return tstr

    def add_section(self, section=None, **kwargs):
        '''
        Add a new section to the template.
        '''
        section = str(section).replace("<", "").replace(">", "")
        if section in self.sections:
            raise AttributeError(
                "Section {section} ".format(section=section) +
                "already exists.")
        str_sec = "<{section}>\n".format(section=section)
        for k, v in kwargs.items():
            str_sec += str(k) + " = " + str(v) + "\n"
        str_sec = _sanitize_script_like(str_sec)

        # put and sanitize
        self.template = self.template + "\n" + str_sec

    def set_value(self, section=None, **kwargs):
        '''
        Set variables in a section.

        If section is 'None' then replacement is attempted over the entire
        template.
        '''
        section = str(section).replace("<", "").replace(">", "")
        fkwargs = _format_dict_filler(kwargs)

        if section == "None":
            self.template = self.template.format_map(fkwargs)
            return
        elif section not in self.sections:
            raise AttributeError(
                "Section {section} ".format(section=section) +
                "is unknown.")

        # modify specified section
        v = self._template_sections[section].format_map(fkwargs)
        self._template_sections[section] = v
        # update template
        self._sections_to_template(populate_template=True)


    @property
    def hash(self):
        # reproducible and update on template mutation
        return self._hash


class problem_compile_flags(object):
    def __init__(self, *args, compile_flags=None, **kwargs):
        self._hash = None
        self._compile_flags = None
        self.compile_flags = compile_flags

    @property
    def compile_flags(self):
        return self._compile_flags

    @compile_flags.setter
    def compile_flags(self, val):
        val = _sanitize_script_like(val)
        cmp_list = [s for s in val.replace('\n', ' ').split(' ')
                    if len(s)>0]
        cmp_list = sorted(cmp_list)

        for ix, el in enumerate(cmp_list):
            if '$' in el:
                ix_el = el.find('$')
                try:
                    cmp_list[ix] = el.replace(el[ix_el:],
                                              os.environ[el[ix_el+1:]])
                except:
                    pass
        self._compile_flags = ' '.join(cmp_list)

        enc = self._compile_flags.encode("utf-8")
        self._hash = _hl.sha256(enc).hexdigest()

    @property
    def compile_flags_as_dod(self):
        # convert to dict of dicts
        od_dod = _c.OrderedDict()
        cf = self.compile_flags.split(' ')
        for fl in cf:
            fl_s = fl.split('=')
            try:
                od_dod[fl_s[0]] = fl_s[1]
            except:
                od_dod[fl_s[0]] = None

        return od_dod

    @property
    def hash(self):
        return self._hash


class problem_IOC_handler(object):
    # provides directory structures and compiler / execution functionality

    _PSPEC = problem_specification
    _PCF = problem_compile_flags
    def __init__(self, *args,
                 problem_specification=None,
                 problem_compile_flags=None,
                 dir_athena=_DIR_ATHENA,
                 dir_outputs=None,
                 dir_outputs_nested_structure=None,
                 fn_base_name=_FN_SCRIPT_BASENAME,
                 verbose=False,
                 make_threads=1,
                 **kwargs):

        if not isinstance(problem_specification, self._PSPEC):
            raise ValueError("")
        else:
            self._problem_specification = problem_specification

        if not isinstance(problem_compile_flags, self._PCF):
            raise ValueError("")
        else:
            self._problem_compile_flags = problem_compile_flags

        self.make_threads = make_threads

        # generate a hash based on inputs
        self._hash = None
        self._mk_hash()

        self.dir_athena = dir_athena
        if dir_outputs is None:
            dir_outputs = os.path.join(self.dir_athena, "outputs")
        self.dir_outputs = dir_outputs
        if dir_outputs_nested_structure is None:
            raise ValueError(
                "dir_outputs_nested_structure must be provided")
        if isinstance(dir_outputs_nested_structure, str):
            dir_outputs_nested_structure = dir_outputs_nested_structure.split(
                os.path.sep)
        t = list(dir_outputs_nested_structure)
        t.extend([self.hash])
        self.dir_outputs_nested_structure = tuple(t)

        self.dir_exe = os.path.join(
            self.dir_outputs,
            *self.dir_outputs_nested_structure[:-1])

        self.fn_base_name = fn_base_name
        self.fn_exec_name = self.fn_base_name + "_" + \
                            self._problem_compile_flags.hash + ".x"

        tmp = self.fn_base_name + "_" + \
              self._problem_compile_flags.hash + ".cmp"
        self.fn_abs_cmp_info = os.path.join(
            self.dir_outputs, *_it.chain(
                self.dir_outputs_nested_structure[:-1], (tmp,)))

        self.fn_abs_problem_input = os.path.join(
            self.dir_outputs, *_it.chain(self.dir_outputs_nested_structure,
                                         ('problem.inp',)))

        self.fn_abs_exec_name = os.path.join(self.dir_exe, self.fn_exec_name)

        # hide random messages
        self._verbose = verbose
        if self._verbose:
            self._quiet_cmd = ""
        else:
            self._quiet_cmd = " > /dev/null 2>&1"

        # add compile hash
        self._problem_specification.set_value(
            section='job',
            compile_hash=self._problem_compile_flags.hash)

        self._problem_specification.set_value(
            section='job',
            base_name=self.fn_base_name)


    def _mk_hash(self):
        # hash here is based on problem_specification and problem_compile_flags
        hps = self._problem_specification.hash
        hcf = self._problem_compile_flags.hash

        enc = (hps + hcf).encode('utf-8')
        self._hash = _hl.sha256(enc).hexdigest()

    def ensure_dir_structure(self):
        '''
        Provide requisite directory structure if it doesn't exist.
        '''

        # record current dir and jump to working directory
        cwd = os.getcwd()
        data_path = os.path.join(
            self.dir_outputs, *self.dir_outputs_nested_structure)

        tar_dir = list(_pl.Path(data_path).parts)

        # attempt to construct subdir structure
        for i in range(1, len(tar_dir)):
            tmp = os.path.join(*tar_dir[:i+1])
            try:
                os.mkdir(tmp)
            except:
                pass

    def clean_athena(self):
        # remove extant objects / binaries
        os.system("rm -rf " + os.path.join(self.dir_athena, "obj") +
                  self._quiet_cmd)
        os.system("rm -rf " + os.path.join(self.dir_athena, "bin") +
                  self._quiet_cmd)

    @property
    def _default_fn_abs_athena(self):
        return os.path.join(self.dir_athena, 'bin', 'athena')

    @property
    def exists_fn_exec_name(self):
        return os.path.exists(self.fn_abs_exec_name)

    @property
    def hash(self):
        return self._hash

    def athena_configure(self):
        '''
        Apply specified configuration to Athena++
        '''
        print("Configuring Athena++")
        # prepare athena configuration for compilation
        cwd = os.getcwd()
        os.chdir(self.dir_athena)
        compile_flags = self._problem_compile_flags.compile_flags
        os.system("python configure.py " + compile_flags +
                  self._quiet_cmd)

        os.chdir(cwd)
        print("Done!")

    def athena_compile(self):
        '''
        Compile Athena++ and move binaries to appropriate locations.
        '''
        print("Compiling Athena++... make_threads={make_threads}".format(
            make_threads=self.make_threads))
        print("flags:")
        print(self._problem_compile_flags.compile_flags)
        # prepare athena configuration for compilation
        cwd = os.getcwd()
        os.chdir(self.dir_athena)

        os.system("make -j" + str(self.make_threads) +
                  self._quiet_cmd)

        # move to the outputs / run directory
        os.system("mv {loc} {tar}".format(
            loc=self._default_fn_abs_athena,
            tar=self.fn_abs_exec_name))

        # cleanup
        self.clean_athena()
        os.chdir(cwd)
        print("Done!")

    def provide_problem_input(self, write_json=True):
        # ensure that we have dumped the content of the problem template
        # to the specific directory
        if not os.path.exists(self.fn_abs_problem_input):
            print("Preparing problem input file")
            with open(self.fn_abs_problem_input, 'w') as f:
                f.write(self._problem_specification.template)

        if not write_json:
            print("Done!")
            return

        fn = os.path.join(os.path.dirname(self.fn_abs_problem_input),
                          "problem.json")
        if not os.path.exists(fn):
            print("Preparing problem input file [json]")
            with open(fn, 'w') as f:
                _j.dump(self._problem_specification.template_as_dod, f)

            print("Done!")

    def provide_compile_flags(self, write_json=True):
        # dump compile flags to file
        if not os.path.exists(self.fn_abs_cmp_info):
            print("Dumping compile_flags to file")
            with open(self.fn_abs_cmp_info, 'w') as f:
                f.write(self._problem_compile_flags.compile_flags)
        if not write_json:
            print("Done!")
            return

        fn = self.fn_abs_cmp_info + ".json"
        if not os.path.exists(fn):
            print("Dumping compile_flags to file [json]")

            with open(fn, 'w') as f:
                _j.dump(self._problem_compile_flags.compile_flags_as_dod, f)

            print("Done!")

    def prepare_athena(self):
        # create directory structure and compile athena as required
        # also prepare problem input together with compilation info

        self.ensure_dir_structure()
        if not self.exists_fn_exec_name:
            self.athena_configure()
            self.athena_compile()

        # dump problem input file
        self.provide_problem_input()

        # dump compile flags
        self.provide_compile_flags()

    def get_data_files(self, dir_problem_hash=None):
        # inspect any extant data to determine if execution is required
        if dir_problem_hash is not None:
            dons = _it.chain(self.dir_outputs_nested_structure[:-1],
                             [dir_problem_hash])
        else:
            dons = self.dir_outputs_nested_structure

        data_dir = os.path.join(self.dir_outputs, *dons)

        with os.scandir(path=data_dir) as sdir:
            problem_files = tuple(entry.name for entry in sdir
                                  if entry.is_file() and
                                  len(entry.name) > 5 and
                                  entry.name[-5:] == "athdf")

        return tuple(os.path.join(data_dir, f) for f in problem_files)

    @property
    def get_sim_time(self):
        # target sim time and known sim time
        str_tlim = _r.findall(
            "tlim = (.+?)\n",
            self._problem_specification.template)

        tar_tlim = float(str_tlim[0])

        df = self.get_data_files()
        if len(df) == 0:
            hav_tlim = 0
        else:
            hav_tlim = _ar.athdf(df[-1])['Time']

        return tar_tlim, hav_tlim

    @property
    def is_exec_required(self):
        # check sim_time within tolerance
        tar_tlim, hav_tlim = self.get_sim_time

        rel_diff_tlim = 1 - hav_tlim / tar_tlim

        if rel_diff_tlim > 1e-12:
            return True

        return False

    def execute_local(self, return_only_command=False):
        '''
        Execute problem locally if required or restart
        '''
        # check if execution is required
        exec_required = self.is_exec_required

        if exec_required:
            exec_dir = self.dir_exe
            exec_name = self.fn_exec_name

            data_dir = os.path.join(
                self.dir_outputs, *self.dir_outputs_nested_structure)

            cmds = []
            cwd = os.getcwd()
            cmd_exec = './' + exec_name + \
                      ' -i ' + \
                      os.path.join(data_dir, "problem.inp") + \
                      ' -d ' + data_dir
            cmds.append("cd " + exec_dir)
            cmds.append(cmd_exec)
            cmds.append("cd " + cwd)

            if return_only_command:
                return cmds

            os.chdir(exec_dir)
            os.system(cmd_exec)
            os.chdir(cwd)


    def submit_cluster(self, *args, **kwargs):
        '''
        Submit specified problem for cluster execution
        '''
        raise NotImplementedError("")

    ###########################################################################
    @property
    def known_datasets(self):
        '''
        For the current 'dirs_outputs_nested_structure' locate all known
        datasets.
        '''
        ddir = self.dir_exe
        dirs = []
        with os.scandir(path=ddir) as sdir:
            dirs = [entry.name for entry in sdir if entry.is_dir()]
        files = {d: self.get_data_files(dir_problem_hash=d) for d in dirs}
        return files


    def dataset_parser(self, *args,
                       filter_problem_input=None,
                       filter_dataset=None,
                       key_format=None,
                       **kwargs):
        '''
        Parse known datasets based on filters [these are treated as unwanted
        elements].

        The filter variables must be specified as an iterable of lambdas that
        evaluate to true when operating on the salient target.
        '''
        data_dir = os.path.join(self.dir_outputs,
                                *self.dir_outputs_nested_structure[:-1])

        known_datasets = self.known_datasets

        relevant_datasets = []

        # for problem input as dict
        json_store = {}

        # filter based on the problem input specification
        if filter_problem_input is not None:
            kkd = tuple(known_datasets.keys())  # required as we pop
            for dir in kkd:
                fn = os.path.join(data_dir, dir, 'problem.json')
                with open(fn, 'r') as f:
                    template_dod = _j.load(f)

                # now apply problem_input filters
                for fil in filter_problem_input:
                    if fil(template_dod):
                        # current set satisfies a filter so throw
                        known_datasets.pop(dir)
                        break
                else: # finally
                    # only store if no pop
                    json_store[dir] = template_dod

        # filter based on dataset specification
        if filter_dataset is not None:
            kkd = tuple(known_datasets.keys())  # required as we pop
            for dir in kkd:

                ffns = tuple(fn for fn in known_datasets[dir]
                             if not all(filt(_ar.athdf(fn, quantities=[]))
                                        for filt in filter_dataset))
                relevant_datasets.extend(ffns)
        else:
            for key, fns in known_datasets.items():
                relevant_datasets.extend(fns)

        # accumulate data satisfying conditions based on input key

        # for each file, merge problem spec (based on json) and athena spec
        # unfortunately, this has a potential to clobber some variables

        if key_format is None:
            return tuple(relevant_datasets)

        relevant_data_keyed = {}

        for ds in relevant_datasets:
            # ensure we have a json
            prob_dir = os.path.dirname(ds)
            base_dir = os.path.split(prob_dir)
            ddir = base_dir[-1]
            if ddir not in json_store:
                fn = os.path.join(prob_dir, 'problem.json')
                with open(fn, 'r') as f:
                    json_store[ddir] = _j.load(f)

            # assemble full meta for each set
            data = _ar.athdf(ds, quantities=[])

            full_metadata = {k: v for k, v in json_store[ddir].items()}
            full_metadata.update(data)

            # get compile_flags
            chash = json_store[ddir]['job']['compile_hash']
            bname = json_store[ddir]['job']['base_name']
            cf_json = os.path.join(base_dir[0],
                                   bname + '_' + chash + '.cmp.json')
            with open(cf_json, 'r') as f:
                full_metadata.update({'compile_flags': _j.load(f)})


            # assemble key:
            key = tuple(str(_extract_from_dodn(full_metadata, kf))
             for kf in key_format)


            relevant_data_keyed[key] = ds

        return relevant_data_keyed


    def data_load(self, data_dict, quantities=None, **athdf_kwargs):
        '''
        Load data into dictionary
        '''
        return {key: _ar.athdf(fn, quantities=quantities, **athdf_kwargs)
                for key, fn in data_dict.items()}

#
# :D
#
