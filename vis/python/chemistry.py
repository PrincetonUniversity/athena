"""
Utilities for reading chemistry output data files
"""


# ========================================================================================

def get_gow17_fields(data, xC_tot=1.6e-4, xO_tot=3.2e-4):
    """Calculate the species abundances from conservation and gas temperature
       in the gow17 network."""
    mh = 1.673534e-24       # hydrogen mass CGS
    km = 1.0e5              # km in cm
    kb = 1.380658e-16       # Boltzmann's const CGS
    gm1 = 5./3. - 1.
    unit_density = 1.4 * mh
    unit_velocity = km
    unit_energy_density = unit_density * unit_velocity**2
    # species
    data["species"] = ["He+", "OHx", "CHx", "CO", "C+", "HCO+", "H2", "H+",
                       "H3+", "H2+", "O+", "Si+"]
    # add ghost species from elemental conservation
    data["species_all"] = data["species"] + ["e", "H", "C", "O"]
    data["re"] = (data["rHe+"] + data["rC+"] + data["rHCO+"]
                  + data["rH+"] + data["rH3+"] + data["rH2+"])
    data["rH"] = 1.0 - (data["rH2"]*2 + 3*data["rH3+"] + 2*data["rH2+"]
                        + data["rH+"] + data["rHCO+"] + data["rCHx"] + data["rOHx"])
    data["rC"] = xC_tot - data["rCHx"] - data["rCO"] - data["rC+"] - data["rHCO+"]
    data["rO"] = xO_tot - data["rOHx"] - data["rCO"] - data["rHCO+"]
    # temperature in Kelvin
    if "press" in data.keys():
        Cv = 1.5 * kb * (1.1 + data["re"] - data["rH2"])
        data["T"] = data["press"]/data["rho"]/gm1*unit_energy_density/Cv
    return
