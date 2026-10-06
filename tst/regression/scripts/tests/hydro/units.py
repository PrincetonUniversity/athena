# Regression test for the units module
#
# Runs the unit_test problem generator with each preset unit system and several
# custom_basis/custom systems, parses the code and basis units printed at setup, and
# compares them to values computed independently here. Also checks that invalid
# <units> input is rejected.

# Modules
import logging
import os
import re
import subprocess
import scripts.utils.athena as athena

logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

# c.g.s. constants, as in src/units/units.hpp
pc = 3.08567758e18
kpc = 3.08567758e21
au = 1.495978707e13
yr = 3.15576e7
Myr = 3.15576e13
Msun = 1.9884099e33
mH = 1.6733e-24
c = 2.99792458e10

rtol = 1.0e-5
_input = os.path.join(athena.athena_rel_path, 'inputs', 'hydro', 'athinput.unit_test')
_run_args = ['mesh/nx1=4', 'mesh/nx2=4', 'mesh/nx3=4', 'time/nlim=0',
             'output1/file_type=vtk']


def expected(L, T, M, mu=1.4, lunit=None, tunit=None, munit=None, vunit=None,
             nunit=None):
    """Code units in c.g.s. and, optionally, basis values in the given units."""
    out = {'code Length': L, 'code Time': T, 'code Mass': M,
           'code density': M / L**3, 'code velocity': L / T,
           'code pressure': M / L / T**2}
    if lunit:
        out['basis length'] = L / lunit
    if tunit:
        out['basis time'] = T / tunit
    if munit:
        out['basis mass'] = M / munit
    if vunit:
        out['basis velocity'] = L / T / vunit
    if nunit:
        out['basis ndensity'] = M / (mu * mH * L**3) / nunit
    return out


# (name, <units> block, expected values)
_cases = [
    ('ism', 'unit_system = ism',
     expected(pc, pc / 1.0e5, 1.4 * mH * pc**3, tunit=Myr, munit=Msun, nunit=1.0)),
    ('galaxy', 'unit_system = galaxy',
     expected(kpc, Myr, 1.4 * mH * kpc**3, vunit=1.0e5, munit=Msun)),
    ('galaxypc', 'unit_system = galaxypc',
     expected(pc, Myr, 1.4 * mH * pc**3, vunit=1.0e5, munit=Msun)),
    ('ism_SR', 'unit_system = ism_SR',
     expected(pc, pc / c, 1.4 * mH * pc**3, tunit=Myr, vunit=1.0e5)),
    ('cgs', 'unit_system = cgs',
     expected(1.0, 1.0, 1.0, vunit=1.0, nunit=1.0)),
    ('SI', 'unit_system = SI',
     expected(1.0e2, 1.0, 1.0e3, vunit=1.0e2, nunit=1.0e-6)),
    ('custom_basis_km',
     'unit_system = custom_basis\nlength = 2.0\nlength_unit = km\n'
     'time = 3.0\ntime_unit = s\nmass = 5.0\nmass_unit = kg',
     expected(2.0e5, 3.0, 5.0e3, lunit=1.0e5, tunit=1.0, munit=1.0e3, vunit=1.0e5,
              nunit=1.0)),
    ('custom_basis_au',
     'unit_system = custom_basis\nmass_per_hydrogen = 1.0\nlength = 1.0\n'
     'length_unit = au\nvelocity = 30.0\nndensity = 1.0e6\nndensity_unit = n/m^3',
     expected(au, au / 3.0e6, mH * au**3, mu=1.0, tunit=Myr, munit=Msun, nunit=1.0e-6)),
    ('custom_kpc',
     'unit_system = custom\nmass_cgs = 1.0e40\nlength_cgs = 3.08567758e21\n'
     'time_cgs = 3.15576e13',
     expected(kpc, Myr, 1.0e40, lunit=pc, tunit=Myr, munit=Msun, vunit=1.0e5,
              nunit=1.0)),
]

# (name, <units> block) that must make the Units constructor raise a fatal error
_bad_cases = [
    ('bad_system', 'unit_system = not_a_system'),
    ('bad_unit', 'unit_system = custom_basis\nlength = 1.0\nlength_unit = lightyear\n'
     'time = 1.0\nndensity = 1.0'),
    ('time_and_velocity', 'unit_system = custom_basis\nlength = 1.0\ntime = 1.0\n'
     'velocity = 1.0\nndensity = 1.0'),
    ('no_mass_or_ndensity', 'unit_system = custom_basis\nlength = 1.0\ntime = 1.0'),
]

_outputs = {}
_bad_outputs = {}


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(prob='unit_test', **kwargs)
    athena.make()


def write_input(name, units_block):
    with open(_input, 'r') as f:
        text = f.read()
    text = text[:text.index('<units>')] + '<units>\n' + units_block + '\n'
    filename = 'athinput.units_' + name
    with open(os.path.join('bin', filename), 'w') as f:
        f.write(text)
    return filename


def run_athena(filename):
    cmd = ['./athena', '-i', filename] + _run_args + athena.global_run_args
    logger.debug('Executing: ' + ' '.join(cmd))
    return subprocess.run(cmd, cwd='bin', stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, universal_newlines=True)


# Run Athena++
def run(**kwargs):
    for name, units_block, _ in _cases:
        result = run_athena(write_input(name, units_block))
        if result.returncode != 0:
            raise athena.AthenaError('unit system ' + name + ' failed:\n'
                                     + result.stdout)
        _outputs[name] = result.stdout
    for name, units_block in _bad_cases:
        _bad_outputs[name] = run_athena(write_input(name, units_block)).stdout


def parse(output):
    values = {}
    for line in output.splitlines():
        m = re.match(r'^(code \w+|basis \w+) = (\S+)', line)
        if m:
            values[m.group(1)] = float(m.group(2))
    return values


# Analyze outputs
def analyze():
    analyze_status = True
    for name, _, ref in _cases:
        values = parse(_outputs[name])
        for key, val in ref.items():
            if key not in values:
                logger.warning('%s: "%s" not found in output', name, key)
                analyze_status = False
            elif abs(values[key] - val) > rtol * abs(val):
                logger.warning('%s: %s = %g, expected %g', name, key, values[key], val)
                analyze_status = False
    # main() catches the exception and still returns 0, so check the message instead
    for name, output in _bad_outputs.items():
        if '### FATAL ERROR in Units' not in output:
            logger.warning('%s: invalid <units> input was accepted', name)
            analyze_status = False
    return analyze_status
