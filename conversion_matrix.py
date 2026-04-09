import pandas as pd
import numpy as np
from itertools import product

# def reverse(string):
#     return(string[::-1])
#
# def converted(string):
#     change_to = {"A":"T", "C":"G", "T": "A", "G": "C"}
#     return("".join([change_to[string[i]] for i in range(len(string))]))
#
# NA = ["A", "T", "G", "C"]
# triplet = ["".join(p) for p in product(NA, repeat = 3)]
#
# conversion_matrix = pd.DataFrame(np.zeros((64,64)),index=triplet, columns=triplet)
# for i in range(64):
#         conversion_matrix.loc[triplet[i], converted(reverse(triplet[i]))] = 1
#
#
#
#
# conversion_matrix.to_csv("conversion_matrix.csv")
#

conversion_matrix = pd.read_csv("conversion_matrix.csv", index_col=0)

def get_triplet_number(triplet):
    return(conversion_matrix.index.get_loc(triplet))

def get_category_name(triplet):
    reverse_converted = conversion_matrix.loc[conversion_matrix[triplet] == 1].index.to_list()[0]
    row_number_q = conversion_matrix.index.get_loc(triplet)
    col_number_rc = conversion_matrix.columns.get_loc(reverse_converted)

    row_number_rc = conversion_matrix.index.get_loc(reverse_converted)
    col_number_q = conversion_matrix.columns.get_loc(triplet)
    if row_number_q <= col_number_rc:
        row = row_number_q
        col = col_number_rc
    else:
        row = row_number_rc
        col = col_number_q
    return reverse_converted, str(row) + ":" + str(col)

def get_category_number (triplet1, triplet2):
    reverse_converted1 = conversion_matrix.loc[conversion_matrix[triplet1] == 1].index.to_list()[0]
    reverse_converted2 = conversion_matrix.loc[conversion_matrix[triplet2] == 1].index.to_list()[0]

    id1 = int(str(get_triplet_number(triplet1)) + str(get_triplet_number(triplet2)))
    id2 = int(str(get_triplet_number(reverse_converted1)) + str(get_triplet_number(reverse_converted2)))
    if id1 < id2:
        id = str(id1) + str(id2)
    else:
        id = str(id2) + str(id1)

    return reverse_converted1 + "_" +  reverse_converted2, id

#
# print(get_category_number("AAT", "AGT"))
# print(get_category_number("ATT", "ACT"))
#



