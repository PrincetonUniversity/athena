'''
extract from io piped
'''
import re
import numpy as np

match_S = "{{{"
match_E = "}}}"


f = open("out", "r").read()


str_match = re.findall('\{(.*?\})', f, re.DOTALL)

