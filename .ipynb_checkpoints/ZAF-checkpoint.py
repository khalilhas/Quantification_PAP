import math
import pandas as pd
import numpy as np
import sqlite3
def get_mac(elements, file_mac="Henke 1993.txt"):
    # import MAC
    mac_df = pd.read_csv(file_mac)
    mac_df_sample = mac_df.query('zAbs in @elements and zMes in @elements')
    return mac_df_sample

def get_element(elements):
    con = sqlite3.connect('xraydb.sqlite')
    cur = con.cursor()
    query = "SELECT * FROM elements WHERE "
    for idx, element in enumerate(elements):
        if idx == len(elements)-1:
            query += "element='"+ str(element) + "'"
        else:
            query += "element='"+ str(element) + "' or "
    df = pd.read_sql_query(query, con)
    con.close()
    return df

# def ZAF_absorption(elements, C, Ec, A, Z):
#     """
#     Attributes
#     ----------
#     elements: list
#         the list of elements present in sample
#     C: dict
#         dictionary of the concentration of elements
#     Ec: dict
#         critical energies for the compounds present in the sample
#         ### Check if you have to add all the lines and what to put for H
#     A: dict
#         the atomic weight of the elements present in the sample
#     Z: dict
#         atomic numbers of the elements present in the sample
#     """
#     get_mac(elements)

# xray database
# con = sqlite3.connect('xraydb.sqlite')
# cur = con.cursor()
# query = """
# SELECT * FROM xray_transitions WHERE """

# ## Experimental conditions
# take_off_angle = 40
# E0 = 5 # keV

# # for the emission of 
# emission = {'element':'Li', 'line':'Ka'}

# # calculating X for each element
# print("[+] The program is considering the emission of ", emission['element'], emission['line'])
# # Extracting the needed MAC from the Henke database
# sample_mac = pd.DataFrame()
# for idx, value in enumerate(elements):
#     new_row = pd.DataFrame({'abs': [value], 'mac': [mac_df[(mac_df['zAbs'] == value) & (mac_df['zMes'] == emission['element'])][emission['line']].values[0]] })
#     print(new_row)
#     sample_mac = sample_mac.append(new_row, ignore_index=True)
# print(sample_mac)
# print(mac_df)