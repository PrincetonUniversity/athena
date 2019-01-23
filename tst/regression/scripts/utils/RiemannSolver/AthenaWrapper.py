from scripts.utils.EquationOfState.eos import Ideal, SimpleHydrogen
from .riemann import riemann_problem
import os
import sys
sys.path.insert(0, '../../../../../vis/python')
import athena_read  # noqa


def athinput2riemann(athinput, eos=None):
    """Uses parameters from "athinput" to define a Riemann problem with a
    specified (optional) eos."""
    if os.path.isfile(athinput):
        athinput = athena_read.athinput(athinput)
    if eos is None:
        # assume ideal EOS with gamma from athinput
        eos = Ideal(athinput['hydro']['gamma'])
    if eos == 'H':
        eos = SimpleHydrogen()
    pin = athinput['problem']
    # only look at d[lr], p[lr], u[lr], T[lr] (state) parameters
    states = {i: pin[i]
              for i in pin if (len(i) == 2) and (i[0] in 'dpuT') and (i[1] in 'lr')}
    return riemann_problem(states, eos)
