# Functions for comparing two datasets

# Modules
import numpy as np
import math


# Function for computing L1 norm of a quantity
def l1_norm(faces, vals):
    return math.fsum(abs(vals) * np.ediff1d(faces))


# Function for computing L1 difference between datasets
def l1_diff(faces_1, vals_1, faces_2, vals_2):
    faces_refined = np.union1d(faces_1, faces_2)

    def fill_to_refined(faces, vals, faces_refined):
        vals_refined = np.empty(len(faces_refined)-1)
        for left_face, right_face, val in zip(faces[:-1], faces[1:], vals):
            indices = np.where(
                (faces_refined >= left_face) & (
                    faces_refined < right_face))[0]
            vals_refined[indices] = val
        return vals_refined
    vals_refined_1 = fill_to_refined(faces_1, vals_1, faces_refined)
    vals_refined_2 = fill_to_refined(faces_2, vals_2, faces_refined)
    return math.fsum(abs(vals_refined_1 - vals_refined_2) * np.ediff1d(faces_refined))
