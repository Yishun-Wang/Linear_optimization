import random
import copy
import math
from gurobipy import *
import random
import itertools
import copy
import math
from shutil import copyfile
import numpy as np
import matplotlib.pyplot as plt
import timeit
import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
from email.MIMEBase import MIMEBase
from email import encoders
from operator import itemgetter
random.seed(1)
# import matplotlib.pyplot as plt


def topology(networkfieldsize, V_E, Number_of_EHnodes, Tx, Data_Links_to_schedule, gamma, N_0, Number_of_Links, Gn_dic, distant_dictionary, P_max):

    set1 = set()
    while len(set1) < len(Tx):
        x = random.uniform(0, networkfieldsize)
        y = random.uniform(0, networkfieldsize)
        coordinate = tuple((x, y))
        set1.add(coordinate)
    # print('set1')

    l1 = list(set1)
    T_dic = {}
    for i in range(len(Tx)):
        T_dic[Tx[i]] = l1[i]

    R_dic = {}
    set2 = set()
    for k, v in T_dic.items():
        Tx_copy = copy.deepcopy(Tx)
        Tx_copy.remove(k)
        while True:
            xr = random.uniform(0, networkfieldsize)
            yr = random.uniform(0, networkfieldsize)
            rloc = tuple((xr, yr))
            r_rest_list = list()
            t_r_dist = math.hypot(xr - v[0], yr - v[1])
            for i in Tx_copy:
                t_loc_rest = T_dic[i]
                r_rest_dist = math.hypot(
                    xr - t_loc_rest[0], yr - t_loc_rest[1])
                r_rest_list.append(r_rest_dist)
            if rloc not in set1 and rloc not in set2 and t_r_dist >= 30.0 and t_r_dist <= 50.0 and all(b > 1.0 for b in r_rest_list):
                set2.add(rloc)
                R_dic[k] = rloc
                break
    # print('set2')
    R_fin_dic = {}
    for k, v in R_dic.items():
        left = len(k) - 3

        if left <= 1:
            R_fin_dic["Rx_" + str(k[3])] = v
        elif left < 3 and left >= 2:
            R_fin_dic["Rx_" + str(k[3]) + str(k[4])] = v
        elif left >= 3:
            R_fin_dic["Rx_" + str(k[3]) + str(k[4]) + str(k[5])] = v

    n_set = set()

    for i in V_E:
        while True:

            xn = random.uniform(0, networkfieldsize)
            yn = random.uniform(0, networkfieldsize)
            dist_list = list()
            for i in Tx:
                t_location = T_dic[i]
                t_n_dist = math.hypot(
                    t_location[0] - xn, t_location[1] - yn)
                dist_list.append(t_n_dist)
            nloc = tuple((xn, yn))
            if nloc not in set1 and nloc not in set2 and nloc not in n_set and all(di > 1.0 for di in dist_list) and any(ddd <= 2.0 for ddd in dist_list):
                n_set.add(nloc)
                break
    # print('n_set')
    l4 = list(n_set)
    N_dic = {}
    for i in range(Number_of_EHnodes):
        N_dic[V_E[i]] = l4[i]
    R_N_dic = dict(N_dic.items() + R_fin_dic.items())
    # print(T_dic)
    # print(R_N_dic)
    distant_dictionary = {}
    for tk, tv in T_dic.items():
        for rnk, rnv in R_N_dic.items():
            distant = math.hypot(tv[0] - rnv[0], tv[1] - rnv[1])
            distant_dictionary[(tk, rnk)] = distant
            Gn_dic[(tk, rnk)] = float(
                1.0 / distant**(2.5))  # generate gain

    G_min = min(Gn_dic.values())

    G_max = max(Gn_dic.values())

    SNR_power = {}
    for i in Data_Links_to_schedule:
        SNR_power[i] = gamma * N_0 / Gn_dic[i]

    Phi = gamma * (N_0 + G_max * (Number_of_Links + len(Tx) - 1)
                   * P_max)
    return Gn_dic, G_max, G_min, SNR_power, Phi, distant_dictionary


def MILP(Time_Links, Time, Data_Links_to_schedule, Gn_dic, gamma, N_0, Phi, Links,
         V_E, alpha, energy_t, P_max, P_min, Time_Data_link, Tx, ave_list, Time_Energy_link,
         f_counter, T, number_co_use_slot, number_pure_data_link_slot,
         number_pure_energy_link_slot, J_e_list, J_d_list, J_co_list, i_co_n, i_d_n, i_e_n,
         i_eff_ie_n, i_eff_id_n, i_eff_co_n):

    m = Model("MILP_debug_model")  # build the model
    m.Params.LogToConsole = 1
    m.Params.OutputFlag = 0  # do not show the log
    m.Params.timeLimit = 500.0
    m.Params.NumericFocus = 3
    # Define variables
    X_t_ij = m.addVars(Time_Links, vtype=GRB.BINARY, name="X_t_ij")
    P_t_ij = m.addVars(Time_Links,
                       vtype=GRB.CONTINUOUS, name="P_t_ij")
    Z_t = m.addVars(Time, vtype=GRB.BINARY, name="Z_t")
    # print("varibales are done ")
    # Define constraints
    c_1 = {}  # data link has to active at least once
    for i in Data_Links_to_schedule:
        c_1[i] = m.addConstr(quicksum(
            [v for k, v in X_t_ij.iteritems() if k[1:] == i]), GRB.GREATER_EQUAL, 1)

    # SINR constraint
    c_2 = {}
    for t in Time:
        for dl in Data_Links_to_schedule:
            FST = Gn_dic[dl] * P_t_ij[(t, dl[0], dl[1])]
            SCD = gamma * quicksum([P_t_ij[(t, odl[0], odl[1])] *
                                    Gn_dic[(odl[0], dl[1])]for odl in Links if odl != dl])

            THD = gamma * N_0
            FTH = Phi * (X_t_ij[(t, dl[0], dl[1])] - 1)
            c_2[(t, dl[0], dl[1])] = m.addConstr(
                (FST - SCD - THD), GRB.GREATER_EQUAL, FTH)

    # Energy constraint
    c_3 = {}
    energy_from_each_link = []
    for eh_node in V_E:
        for l in Links:
            one_link_energy_in_whole_frame = alpha * Gn_dic[(l[0], eh_node)] *\
                quicksum([P_t_ij[(t, l[0], l[1])] for t in Time])
            energy_from_each_link.append(one_link_energy_in_whole_frame)
        c_3[eh_node] = m.addConstr(
            quicksum(energy_from_each_link), GRB.GREATER_EQUAL, energy_t * 1000.0)
        energy_from_each_link = []

   # Activation constratin
    c_4 = {}
    for t in Time:
        for l in Links:
            c_4[(t, l[0], l[1])] = m.addConstr(
                Z_t[t], GRB.GREATER_EQUAL, X_t_ij[t, l[0], l[1]])

   # Two power constraints

    c_5 = {}
    for tl in Time_Links:
        c_5[tl] = m.addConstr(
            P_t_ij[tl], GRB.LESS_EQUAL, P_max * X_t_ij[tl])

    c_5_scd = {}
    for tl in Time_Data_link:
        lower_bond = gamma * N_0 / Gn_dic[(tl[1], tl[2])]
        c_5_scd[tl] = m.addConstr(
            # SNR_power[(tl[1],tl[2])]*X_t_ij[tl], GRB.LESS_EQUAL, P_t_ij[tl])
            P_min * X_t_ij[tl], GRB.LESS_EQUAL, P_t_ij[tl])

    # Half Duplex constraint

    c_6 = {}
    for t in Time:
        for transmitter in Tx:
            fromx = sum(X_t_ij.select(t, transmitter, "*"))
            c_6[(t, transmitter)] = m.addConstr(fromx, GRB.LESS_EQUAL, 1)

    # There is no c7 only c8

    # Number of time slots constraint
    c_8 = []
    c_8.append(m.addConstr(quicksum(Z_t[t]
                                    for t in Time), GRB.LESS_EQUAL, T))
    # Set Objective
    obj = []
    m.setObjective(quicksum(Z_t[t] for t in Time), GRB.MINIMIZE)
    sum_z = sum(Z_t[t] for t in Time)
    obj.append(sum_z)

    # Update and process
    m.update()
    m.write('main_1.lp')
    m.optimize()
    print(m.status)
    # Results generation
    # 1. Different type of time slots
    if m.status == 2:
        ave_list.append(m.objVal)
        X_t_ij_copy = {}
        schedule = {}
        P_dic = {}
        X_dic = {}
        P_e_dic = {}
        X_e_dic = {}
        pure_data_link_slot = []
        pure_energy_link_slot = []
        co_use_slot = []
        J_energy = []
        J_data = []
        J_co = []
        for i in Time_Links:
            temp = round(X_t_ij[i].X)
            X_t_ij_copy[i] = temp

        for t in Time_Data_link:
            if X_t_ij_copy[t] == 1.0 and P_t_ij[t].X > 0.0:
                X_dic[t] = X_t_ij_copy[t]
                P_dic[t] = P_t_ij[t].X
        for l in Time_Energy_link:
            if X_t_ij_copy[l] == 1.0 and P_t_ij[l].X > 0.0:
                P_e_dic[l] = P_t_ij[l].X
                X_e_dic[l] = 1.0
        for tl in Time:
            if round(Z_t[tl].X) > 0.0:  # here  it has to be int because the results is int
                schedule[tl] = []

        for tl in Time_Links:
            if X_t_ij_copy[tl] == 1.0 and P_t_ij[tl].X > 0.0:
                schedule[tl[0]].append(
                    tuple((tl[0], tl[1], tl[2], P_t_ij[tl].X)))

        for slot, values in schedule.items():
            if len(values) == 0.0:
                del schedule[slot]
        multiple_link_dic = {}
        for k, v in schedule.items():
            if len(v) > 1:
                multiple_link_dic[k] = v
        schedule_copy = copy.deepcopy(schedule)
        for i in multiple_link_dic:
            del schedule_copy[i]
        d_link_set = set(P_dic.keys())
        e_link_set = set(P_e_dic.keys())
        for vl in schedule_copy.values():
            v = vl[0]
            if tuple((v[0], v[1], v[2])) in d_link_set:
                pure_data_link_slot.append(tuple((v[0], v[1], v[2])))

            elif tuple((v[0], v[1], v[2])) in e_link_set:
                pure_energy_link_slot.append(tuple((v[0], v[1], v[2])))

        cc1 = 0
        cc2 = 0
        cc3 = 0
        for i in multiple_link_dic.values():
            third_position = []
            for j in i:
                third_position.append(j[2])
            if all(three == 'refpoint' for three in third_position):
                cc1 += 1
                for j in i:

                    pure_energy_link_slot.append(tuple((j[0], j[1], j[2])))
            elif all(three != 'refpoint' for three in third_position):
                cc2 += 1
                for j in i:
                    pure_data_link_slot.append(tuple((j[0], j[1], j[2])))
            else:
                cc3 += 1
                for j in i:
                    co_use_slot.append(tuple((j[0], j[1], j[2])))
        pure_energy_link_set = set()
        pure_data_link_set = set()
        co_use_set = set()
        for i in pure_energy_link_slot:
            pure_energy_link_set.add(i[0])
        for i in pure_data_link_slot:
            pure_data_link_set.add(i[0])
        for i in co_use_slot:
            co_use_set.add(i[0])
        number_pure_energy_link_slot.append(len(pure_energy_link_set))
        number_pure_data_link_slot.append(len(pure_data_link_set))
        number_co_use_slot.append(len(co_use_set))
        # Energy
        for i in pure_energy_link_slot:
            J_energy.append(P_e_dic[i])
        for i in co_use_slot:
            J_co.append(P_t_ij[i].X)
        for i in pure_data_link_slot:
            J_data.append(P_dic[i])
        J_e_list.append(sum(J_energy))
        J_d_list.append(sum(J_data))
        J_co_list.append(sum(J_co))
        i_e = []
        i_d = []
        i_co = []
        i_eff_ie = []
        i_eff_id = []
        i_eff_co = []
        for i in V_E:
            n_i_e = []
            n_i_d = []
            n_i_co = []

            for lk in pure_energy_link_slot:

                i_e.append(P_e_dic[lk] * Gn_dic[(lk[1], i)] * alpha)
                n_i_e.append(P_e_dic[lk] * Gn_dic[(lk[1], i)] * alpha)
            for lk in pure_data_link_slot:
                i_d.append(P_dic[lk] * Gn_dic[(lk[1], i)] * alpha)
                n_i_d.append(P_dic[lk] * Gn_dic[(lk[1], i)] * alpha)
            for lk in co_use_slot:
                n_i_co.append(P_t_ij[lk].X * Gn_dic[(lk[1], i)] * alpha)
                i_co.append(P_t_ij[lk].X * Gn_dic[(lk[1], i)] * alpha)
            i_eff_ie.append(
                sum(n_i_e) / (sum(n_i_e) + sum(n_i_d) + sum(n_i_co)))
            i_eff_id.append(
                sum(n_i_d) / (sum(n_i_e) + sum(n_i_d) + sum(n_i_co)))
            i_eff_co.append(
                sum(n_i_co) / (sum(n_i_e) + sum(n_i_d) + sum(n_i_co)))

        if len(V_E) != 0.0:
            i_co_n.append(float(sum(i_co) / len(V_E)))
            i_e_n.append(float(sum(i_e) / len(V_E)))
            i_d_n.append(float(sum(i_d) / len(V_E)))
            i_eff_ie_n.append(sum(i_eff_ie) / len(V_E))
            i_eff_id_n.append(sum(i_eff_id) / len(V_E))
            i_eff_co_n.append(sum(i_eff_co) / len(V_E))
        pure_data_link_slot = []
        pure_energy_link_slot = []
        co_use_slot = []
        J_energy = []
        J_data = []
        J_co = []

    else:
        f_counter += 1

    if m.status == 2:
        return ave_list, X_dic, P_dic, X_e_dic, P_e_dic, schedule, J_co_list, J_d_list,\
            J_e_list, number_co_use_slot, number_pure_data_link_slot, number_pure_energy_link_slot,\
            f_counter, i_co_n, i_d_n, i_e_n, i_eff_ie_n, i_eff_id_n, i_eff_co_n

    else:
        f_counter += f_counter
        # ave_list = []
        # X_dic = {}
        # P_dic = {}
        # X_e_dic = {}
        # P_e_dic = {}
        # schedule = {}
        # J_co_list = []
        # J_d_list = []
        # J_e_list = []
        # number_co_use_slot = []
        # number_pure_data_link_slot = []
        # number_pure_energy_link_slot = []
        # i_co_n = []
        # i_d_n = [] 
        # i_e_n = []
        # i_eff_ie_n = []
        # i_eff_id_n = []
        # i_eff_co_n = []  
        return ave_list, None, None, None, None, None, J_co_list, J_d_list,\
            J_e_list, number_co_use_slot, number_pure_data_link_slot, number_pure_energy_link_slot,\
            f_counter, i_co_n, i_d_n, i_e_n, i_eff_ie_n, i_eff_id_n, i_eff_co_n



def CLLCTMP(l, number_pure_energy_link_slot, number_pure_data_link_slot,
            number_co_use_slot, J_e_list, J_d_list, J_co_list, i_e_n, i_d_n, i_co_n,
            i_eff_co_n, i_eff_id_n, i_eff_ie_n, max_dic, max_co_dic, max_d_dic, max_e_dic,
            min_dic, min_co_dic, min_d_dic, min_e_dic, link_level_J_co, link_level_J_data,
            link_level_J_energy, link_level_schedule_length, link_level_co_length,
            link_level_data_length, link_level_energy_length,
            link_level_n_co, link_level_n_d, link_level_n_e, link_level_n_propo_co,
            link_level_n_propo_d, link_level_n_propo_e, ave_list, samples, link_level_n_harvested):

    if len(ave_list) != 0:
        link_level_schedule_length[l] = sum(
            ave_list) / float(len(ave_list))
        print(link_level_schedule_length)
        max_dic[l] = max(ave_list)
        min_dic[l] = min(ave_list)

    if len(number_pure_data_link_slot) != 0.0:
        link_level_data_length[l] = sum(
            number_pure_data_link_slot) / float(len(number_pure_data_link_slot))
        max_d_dic[l] = max(number_pure_data_link_slot)
        min_d_dic[l] = min(number_pure_data_link_slot)

    else:
        link_level_data_length[l] = 0.0

    if len(number_pure_energy_link_slot) != 0.0:
        link_level_energy_length[l] = sum(
            number_pure_energy_link_slot) / float(len(number_pure_energy_link_slot))
        max_e_dic[l] = max(number_pure_energy_link_slot)
        min_e_dic[l] = min(number_pure_energy_link_slot)

    else:
        link_level_energy_length[l] = 0.0

    if len(number_co_use_slot) != 0.0:
        link_level_co_length[l] = sum(
            number_co_use_slot) / float(len(number_co_use_slot))
        max_co_dic[l] = max(number_co_use_slot)
        min_co_dic[l] = min(number_co_use_slot)

    else:
        link_level_co_length[l] = 0.0

    if len(J_e_list) != 0.0:
        link_level_J_energy[l] = sum(
            J_e_list) / float(len(J_e_list))
    else:
        link_level_J_energy[l] = 0.0
    if len(J_d_list) != 0.0:
        link_level_J_data[l] = sum(
            J_d_list) / float(len(J_d_list))
    else:
        link_level_J_data[l] = 0.0

    if len(J_co_list) != 0:
        link_level_J_co[l] = sum(
            J_co_list) / float(len(J_co_list))
    else:
        link_level_J_co[l] = 0.0

    link_level_n_e[l] = float(sum(i_e_n) / samples)
    link_level_n_d[l] = float(sum(i_d_n) / samples)
    link_level_n_co[l] = float(sum(i_co_n) / samples)
    link_level_n_harvested[l] = float(sum(i_e_n)+sum(i_d_n)+sum(i_co_n))/ samples

    link_level_n_propo_d[l] = float(
        sum(i_eff_id_n) / samples)
    link_level_n_propo_e[l] = float(
        sum(i_eff_ie_n) / samples)
    link_level_n_propo_co[l] = float(
        sum(i_eff_co_n) / samples)

    return link_level_schedule_length, link_level_co_length, link_level_data_length, link_level_J_co,\
        link_level_J_data, link_level_J_energy, link_level_n_co, link_level_n_d, link_level_n_e,\
        link_level_n_propo_co, link_level_n_propo_d, link_level_n_propo_e, max_dic, min_dic, max_co_dic,\
        max_d_dic, max_e_dic, min_co_dic, min_d_dic, min_e_dic


def LP1(P_snr, gamma, N_0, Gn_dic, potential_list, P_min, P_max, transmission_set, T_counter,
        Data_link_power_dic, V_E, E_dic, alpha, energy_t):
    while len(P_snr) > 0:
        # print("1")
        author__ = 'artemr'
        min_psnr = min(P_snr, key=P_snr.get)
        potential_list.append(min_psnr)
        # print("before po = " + str(potential_list))
        # Heuristic_log_file.write("LP1_potentiallist_befor_model=" + str(potential_list))
        # Heuristic_log_file.write('LP1 starts there')
        LP1 = Model("Heur_LP1")
        LP1.Params.timeLimit = 7200.0
        # LP1.Params.NumericFocus = 3
        LP1.Params.LogToConsole = 0
        LP1.Params.OutputFlag = 0
        # LP1.Params.OutputFlag = 0

        D_P_ij = LP1.addVars(
            potential_list, vtype=GRB.CONTINUOUS, name="D_P_ij")

        LP1C1 = {}
        for i in potential_list:
            FST = Gn_dic[i] * D_P_ij[(i[0], i[1])]
            SCD = gamma * quicksum([D_P_ij[(oi[0], oi[1])] * Gn_dic[(oi[0], i[1])]
                                    for oi in potential_list if oi != i])
            THD = gamma * N_0

            LP1C1[(i[0], i[1])] = LP1.addConstr(
                (FST - SCD - THD), GRB.GREATER_EQUAL, 0)

        LP1C2 = {}
        for i in potential_list:
            LP1C2[i] = LP1.addConstr(P_min, GRB.LESS_EQUAL, D_P_ij[i])

        LP1C3 = {}
        for i in potential_list:
            LP1C3[i] = LP1.addConstr(D_P_ij[i], GRB.LESS_EQUAL, P_max)

        LP1.setObjective(quicksum(D_P_ij[i]
                                  for i in potential_list), GRB.MINIMIZE)
        LP1OBJ = quicksum(D_P_ij[i] for i in potential_list)

        LP1.update()
        LP1.write('testlp1.lp')
        LP1.optimize()

        if LP1.status == 3:

            # print('failed_po' + str(potential_list))

            transmission_set[T_counter] = potential_list[0:len(
                potential_list) - 1]
            for i in potential_list[0:len(potential_list) - 1]:
                if i in P_snr.keys():
                    del P_snr[i]

            potential_list = []
            T_counter += 1

        elif LP1.status == 2:

            # print("success_pl" + str(potential_list))

            for i in potential_list:
                if i in P_snr.keys():
                    del P_snr[i]
            for i in potential_list:

                Data_link_power_dic[(T_counter, i[0], i[1])] = D_P_ij[i].X

            if len(P_snr) == 0:
                transmission_set[T_counter] = potential_list
                break

    TLP = len(transmission_set)
    for n in V_E:
        E_dic[n] = energy_t * 1000
    E_dic_copy = copy.deepcopy(E_dic)
    # print('before = ' + str(E_dic))

    for i in Data_link_power_dic:
        for n in V_E:
            E_dic_copy[n] -= Gn_dic[(i[1], n)] * \
                Data_link_power_dic[i] * alpha

    return transmission_set, TLP, T_counter, Data_link_power_dic, E_dic, E_dic_copy


def LP2(transmission_set, Gn_dic, gamma, alpha,
        N_0, Data_link_power_dic,
        Energy_Links, P_max, P_min, LP2schedule, T_counter, E_dic,
        power_schedule, energy_t, V_E, Number_of_Links, energy_transmitter_set):
    for each_slot in transmission_set:
        LP2 = Model("Heur_LP2")
        LP2.Params.timeLimit = 7200.0
        # LP2.Params.NumericFocus = 3
        LP2.Params.LogToConsole = 0
        LP2.Params.OutputFlag = 0
        # 1. take the data link to make as decicions variable
        D_P_xy = LP2.addVars(
            transmission_set[each_slot], vtype=GRB.CONTINUOUS, name="D_P_xy")

        # 2. take the potntail energy link
        selected_energy_list = []
        data_transmit_set = set()
        for link in transmission_set[each_slot]:
            data_transmit_set.add(link[0])
        potential_set = energy_transmitter_set - data_transmit_set
        for i in Energy_Links:
            if i[0] in potential_set:
                selected_energy_list.append(i)
        E_P_xy = LP2.addVars(selected_energy_list,
                             vtype=GRB.CONTINUOUS, name="E_P_xy")
        LP2.update()
        SINR_check_list = selected_energy_list + \
            transmission_set[each_slot]

        LP2C1 = {}  # SINR constraint
        for i in transmission_set[each_slot]:
            FST = Gn_dic[i] * D_P_xy[i]
            SCD1 = gamma * quicksum([D_P_xy[(oi[0], oi[1])] * Gn_dic[(oi[0], i[1])]
                                     for oi in transmission_set[each_slot] if oi != i])
            SCD2 = gamma * quicksum([E_P_xy[ei] * Gn_dic[(ei[0], i[1])]
                                     for ei in selected_energy_list])
            THD = gamma * N_0
            LP2C1[i] = LP2.addConstr(
                (FST - SCD1 - SCD2 - THD), GRB.GREATER_EQUAL, 0)
        LP2C2 = {}  # power constraint for data links
        LP2C3 = {}
        for i in transmission_set[each_slot]:
            LP2C2[i] = LP2.addConstr(Data_link_power_dic[(
                each_slot, i[0], i[1])], GRB.LESS_EQUAL, D_P_xy[i])
            LP2C3[i] = LP2.addConstr(
                D_P_xy[i], GRB.LESS_EQUAL, P_max)
        LP2C4 = {}  # power constratin for energy link
        LP2C5 = {}
        for i in selected_energy_list:
            LP2C4[i] = LP2.addConstr(0, GRB.LESS_EQUAL, E_P_xy[i])
            LP2C5[i] = LP2.addConstr(
                E_P_xy[i], GRB.LESS_EQUAL, P_max)

        Data_power = quicksum([D_P_xy[i]
                               for i in transmission_set[each_slot]])
        Energy_power = quicksum([E_P_xy[ei]
                                 for ei in selected_energy_list])
        LP2.setObjective((Data_power + Energy_power), GRB.MAXIMIZE)
        LP2.update()
        LP2.write("LP2.lp")
        LP2.optimize()
        # print("LP2 status="+str(LP2.status))
        if LP2.status == 3:
            print('failed')
        elif LP2.status == 2:

            chosen_energy_list = []
            for i in selected_energy_list:
                if E_P_xy[i].X > 0:
                    chosen_energy_list.append(i)

            LP2schedule[each_slot] = list(
                D_P_xy.keys() + chosen_energy_list)

            for i in selected_energy_list:
                if E_P_xy[i].X > 0:

                    power_schedule[(each_slot, i[0], i[1])
                                   ] = E_P_xy[i].X
            for i in D_P_xy:
                power_schedule[(each_slot, i[0], i[1])
                               ] = D_P_xy[i].X
        # print('\n')

    T_DATA_COUNTER = T_counter

    for lk in power_schedule:
        for cn in V_E:
            E_dic[cn] -= Gn_dic[(lk[1], cn)] * \
                power_schedule[lk] * alpha
    while any(ele > 0 for ele in E_dic.values()):
        T_counter += 1
        LP2schedule[T_counter] = ' all energy links '
        power_schedule[(T_counter, 'all_energy_link')
                       ] = Number_of_Links * P_max
        for n in V_E:
            for elk in Energy_Links:
                E_dic[n] -= Gn_dic[(elk[0], n)] * P_max * alpha
    fn_schedule = LP2schedule
    fn_power = power_schedule

    return fn_schedule, E_dic, T_DATA_COUNTER, T_counter, fn_power


def CLLCTLA(T_DATA_COUNTER, fn_schedule, fn_power, co_slot_number,
            data_slot_number, energy_slot_number, fn_slot_number, T_counter, Energy_Links,
            E_dic, ave_n_harvested_list, ave_n_co_list, ave_n_data_list, ave_n_energy_list,
            Gn_dic, alpha, Number_of_EHnodes, P_max, Tx):
    co_counter = 0
    d_counter = 0
    # e_counter =0
    co_name_list = {}

    co_slot = []

    data_slot = []

    energy_slot = []

    # print(fn_schedule)
    for i in range(1, T_DATA_COUNTER + 1):
        # print(i)
        links_in_slot = fn_schedule[i]
        # print('links_in_slot=' + str(links_in_slot))
        second_part = []
        for lk in links_in_slot:
            second_part.append(lk[1])
        # print('second_part=' + str(second_part))
        if any(second == 'refpoint' for second in second_part):
            # print('DOES IF WORKS?')
            # co_counter += 1
            # co_slot_number.append(co_counter)
            co_name_list[i] = fn_schedule[i]
            for lkc in links_in_slot:
                co_slot.append((i, lkc[0], lkc[1]))
        elif all(second != 'refpoint' for second in second_part):
            # d_counter += 1
            # data_slot_number.append(d_counter)

            for lkd in links_in_slot:
                data_slot.append((i, lkd[0], lkd[1]))

    co_counter = len(co_name_list)
    co_slot_number.append(co_counter)
    d_counter = len(data_slot)
    data_slot_number.append(d_counter)

    left_e = T_counter - T_DATA_COUNTER
    if left_e > 0:
        e_counter = left_e
        energy_slot_number.append(e_counter)
        for i in range(T_DATA_COUNTER + 1, T_counter + 1):
            energy_slot.append((i, 'all energylinks'))
    else:
        e_counter = 0
        energy_slot_number.append(e_counter)

    fn_slot_number.append(len(fn_schedule))

######################################################### initial schedule collection stops here ########

    harvested_energy = {}
    from_co = {}
    from_energy = {}
    from_data = {}
    total_slot_record = {}
    for i in range(T_DATA_COUNTER + 1, T_counter + 1):
        for elink in Energy_Links:
            total_slot_record[(i, elink[0], elink[1])] = P_max
    non_energy_list = data_slot + co_slot

    for i in non_energy_list:
        total_slot_record[i] = fn_power[i]

    # Heuristic_log_file.write(
    #     'total_slot_record=' + str(total_slot_record))
    for n in E_dic:
        from_co[n] = sum([fn_power[coslot] * alpha *
                          Gn_dic[(coslot[1], n)] for coslot in co_slot])
        from_data[n] = sum([fn_power[dslot] * alpha *
                            Gn_dic[(dslot[1], n)] for dslot in data_slot])
        from_energy[n] = left_e * P_max * alpha * \
            (sum([Gn_dic[(transmitters, n)] for transmitters in Tx]))
        harvested_energy[n] = sum(
            [total_slot_record[lkslot] * alpha * Gn_dic[(lkslot[1], n)] for lkslot in total_slot_record])

    average_n_harvested = sum(
        harvested_energy.values()) / float(Number_of_EHnodes)
    average_n_co = sum(from_co.values()) / float(Number_of_EHnodes)
    average_n_data = sum(from_data.values()) / float(Number_of_EHnodes)
    average_n_energy = sum(from_energy.values()) / \
        float(Number_of_EHnodes)
    ave_n_harvested_list.append(average_n_harvested)
    ave_n_co_list.append(average_n_co)
    ave_n_data_list.append(average_n_data)
    ave_n_energy_list.append(average_n_energy)
    return fn_slot_number, co_slot_number, energy_slot_number, data_slot_number,\
        ave_n_harvested_list, ave_n_co_list, ave_n_data_list,\
        ave_n_energy_list


def CLLCTLA2(fn_slot_number, samples, data_slot_number, energy_slot_number,
             co_slot_number, ave_n_energy_list, ave_n_data_list, ave_n_co_list, ave_n_harvested_list,
             ave_schedule, ave_data_time_slot, ave_co_time_slot, ave_energy_time_slot,
             n_energy, n_data, n_co, n_harvested, Number_of_Links,fn_max_dic, fn_min_dic):

    ave_schedule[Number_of_Links] = sum(fn_slot_number) / float(samples)
    ave_data_time_slot[Number_of_Links] = sum(
        data_slot_number) / float(samples)
    ave_co_time_slot[Number_of_Links] = sum(
        co_slot_number) / float(samples)
    ave_energy_time_slot[Number_of_Links] = sum(
        energy_slot_number) / float(samples)
    n_energy[Number_of_Links] = sum(ave_n_energy_list) / float(samples)
    n_data[Number_of_Links] = sum(ave_n_data_list) / float(samples)
    n_co[Number_of_Links] = sum(ave_n_co_list) / float(samples)
    n_harvested[Number_of_Links] = sum(
        ave_n_harvested_list) / float(samples)
    fn_max_dic[Number_of_Links] = max(fn_slot_number)
    fn_min_dic[Number_of_Links] = min(fn_slot_number)
    return ave_schedule, ave_data_time_slot, ave_energy_time_slot, ave_co_time_slot,\
        n_energy, n_co, n_data, n_harvested, fn_max_dic, fn_min_dic
