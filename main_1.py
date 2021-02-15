
from gurobipy import *
import random
import itertools
import copy
import math
import matplotlib.pyplot as plt
from decimal import *
from shutil import copyfile
import numpy as np
import matplotlib.pyplot as plt
import smtplib
# from email.MIMEMultipart import MIMEMultipart
# from email.MIMEText import MIMEText
# from email.MIMEBase import MIMEBase
# from email import encoders
import timeit
import numpy as np
from functions import topology
random.seed(1)


# 1. Enironment paramter set up stage

T = int(15)  # 1.Given time slot value
# alpha = 0.85  # 2. Harvesting efficiency cited from other paper, do not change
gamma = float(100.0)  # 3. the SINR threshold, equal to 20 dB
# 4. Thermal noise. As we do not have a specific platlfrom hence averagely
# is -90 dBm which is 10 **-9 mW
N_0 = float(1 * 10 ** -9)
# / # 5. The energy requirment of EH nodes
P_max = float(5000.0)  # 6. Maximum power in milliwatt
P_min = float(750.0)
# e_start = 1.0
# e_end = 1.1
# e_step = 0.1
energy_t = 0.1
# The following are network specifc parameters
loop_time = 1
Number_of_EHnodes = 1
networkfieldsize = 80.0  # meters This set the topology size
start_number_of_links = 5
end_number_of_links = 5
alpha = 1.00
c_counter = 0
step = 1
# start_time = timeit.default_timer()  # start to record time

g_loop = [100]
#  empty dictaionries and list variables which are saved for later use

max_dic = dict()  # save maximum time slot
min_dic = dict()  # save minimum time slot
max_d_dic = {}
min_d_dic = {}
max_e_dic = {}
min_e_dic = {}
max_co_dic = {}
min_co_dic = {}
ave_dic = dict()  # save average time slot
ave_list = list()  # store m.objVal
pure_energy_link_slot = []
number_pure_energy_link_slot = []
pure_data_link_slot = []
number_pure_data_link_slot = []
co_use_slot = []
number_co_use_slot = []
J_energy = []
J_data = []
J_co = []
J_e_list = []
J_d_list = []
J_co_list = []
i_e_n = []
i_d_n = []
i_co_n = []
schedule_ratio_e = {}
schedule_ratio_d = {}
schedule_ratio_co = {}
energy_ratio_e = {}
energy_ratio_d = {}
energy_ratio_co = {}
Gn_dic = {}  # channge gain dictionary
distant_dictionary = {}  # store distanve between two nodes
link_level_schedule_length = {}
link_level_energy_length = {}
link_level_data_length = {}
link_level_co_length = {}
link_level_J_energy = {}
link_level_J_data = {}
link_level_J_co = {}
link_level_n_e = {}
link_level_n_d = {}
link_level_n_co = {}
link_level_n_propo_e = {}
link_level_n_propo_d = {}
link_level_n_propo_co = {}
i_eff_ie_n = []
i_eff_id_n = []
i_eff_co_n = []
'''
Note, during the loop, some of the precreate dict or set has to be
cleaned
Some of them are just for temperoray use, they can be defined on the way
'''


# 2. NETWORK SET UP
for gl in g_loop:
    gamma = gl
    f_counter = 0
    start_time = timeit.default_timer()
    for Given in range(start_number_of_links, end_number_of_links, step):

        Number_of_Links = Given
        start_time = timeit.default_timer()
        V_E = ['n_' + str(EHnodes)for EHnodes in range(Number_of_EHnodes)]
        Tx = ['Tx_' + str(nodes) for nodes in range(Number_of_Links)]
        Rx = ['Rx_' + str(nodes) for nodes in range(Number_of_Links)]
        Recevie_nodes = Rx + V_E
        Number_of_nodes = Number_of_EHnodes + \
            Number_of_Links * 2
        The_total_nodes_in_network = Tx + Rx + V_E
        R_N = Rx + V_E
        All_links = tuplelist(itertools.product(Tx, R_N))
        for i in All_links:
            Gn_dic[i] = None

        # Links in the network

        Data_Links_to_schedule = [('Tx_' + str(Lnkindex), 'Rx_' + str(Lnkindex))
                                  for Lnkindex in range(Number_of_Links)]
        Energy_Links = [(i, "refpoint") for i in Tx if Number_of_EHnodes != 0]
        Links = Data_Links_to_schedule + Energy_Links

        # assign time slot to links

        Time = tuplelist(t + 1 for t in range(T))
        Time_Data_link = tuplelist([(t, dlink[0], dlink[1])
                                    for t in Time for dlink in Data_Links_to_schedule])

        Time_Energy_link = tuplelist([(t, elink[0], elink[1])
                                      for t in Time for elink in Energy_Links])

        Time_Links = tuplelist([(t, l[0], l[1])
                                for t in Time for l in Links])
        print("set up generated")

        for loop in range(loop_time):
            # print('alpha='+str(alpha))
            print('numberoflink='+str(Number_of_Links))
            print('loop='+str(loop))
            print("in innner loop + '\n'+current loop number = " + str(loop))
            set1 = set()
            while len(set1) < len(Tx):
                x = random.uniform(0, networkfieldsize)
                y = random.uniform(0, networkfieldsize)
                coordinate = tuple((x, y))
                set1.add(coordinate)

            l1 = list(set1)
            T_dic = {}
            for i in range(len(Tx)):
                T_dic[Tx[i]] = l1[i]

            R_dic = {}
            set2 = set()
            l3 = []

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

            l4 = list(n_set)
            N_dic = {}
            for i in range(Number_of_EHnodes):
                N_dic[V_E[i]] = l4[i]
            R_N_dic = dict(N_dic.items() + R_fin_dic.items())

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
                           * P_max)  # Phi is the only parameter\

            # only for plot purpose
            topology_f = plt.figure(3)
            plt.clf()
            plt.axis([0, networkfieldsize, 0, networkfieldsize])

            T_x_axis = []
            T_y_axis = []
            for i in T_dic:
                coordi = T_dic[i]
                x_value = coordi[0]
                T_x_axis.append(x_value)
                y_value = coordi[1]
                T_y_axis.append(y_value)
                plt.annotate(i, (x_value, y_value), size=15)
                plt.plot(T_x_axis, T_y_axis, '^', ms=10, c='r')

            R_x_axis = []
            R_y_axis = []
            for i in R_fin_dic:
                co = R_fin_dic[i]
                x_r_value = co[0]
                R_x_axis.append(x_r_value)
                y_r_value = co[1]
                R_y_axis.append(y_r_value)
                plt.annotate(i, (x_r_value, y_r_value), size=15)
                plt.plot(R_x_axis, R_y_axis, 'o', ms=10, c='blue')

            n_x_axis = []
            n_y_axis = []
            for i in N_dic:
                c = N_dic[i]
                x_n_value = c[0]
                n_x_axis.append(x_n_value)
                y_n_value = c[1]
                n_y_axis.append(y_n_value)
                plt.annotate(i, (x_n_value, y_n_value), size=15)
                plt.plot(n_x_axis, n_y_axis, 's', ms=10, c='green')

            for i in Tx:
                t = T_dic[i]

                r = R_dic[i]

                plt.plot((t[0], r[0]), (t[1], r[1]),
                         color='black', linestyle='solid')
            plt.show()
            Gn_dic, G_max, G_min, SNR_power, Phi, distant_dictionary = topology(networkfieldsize, V_E, Number_of_EHnodes,
                                                                                Tx, Data_Links_to_schedule, gamma, N_0, Number_of_Links,
                                                                                Gn_dic, distant_dictionary, P_max)

            # print("start to generate optimization model ")

            # plot end

            # END OF NETWORK TcOPOLOGY

            # 3. OPTIMIZATION MODEL START NOW
            m = Model("MILP_debug_model")  # build the model
            m.Params.LogToConsole = 1
            m.Params.OutputFlag = 0  # do not show the log
            m.Params.timeLimit = 100.0
            # m.Params.NumericFocus = 3
            X_t_ij = m.addVars(Time_Links, vtype=GRB.BINARY, name="X_t_ij")
            P_t_ij = m.addVars(Time_Links,
                               vtype=GRB.CONTINUOUS, name="P_t_ij")
            Z_t = m.addVars(Time, vtype=GRB.BINARY, name="Z_t")
            # print("varibales are done ")
            # constraintslen9
            c_1 = {}  # data link has to active at least once
            for i in Data_Links_to_schedule:
                c_1[i] = m.addConstr(quicksum(
                    [v for k, v in X_t_ij.iteritems() if k[1:] == i]), GRB.GREATER_EQUAL, 1)
            # print("c_1 is done ")
            # ------------------SINR constraint --------------
            ''''
l1=[]
				for k ,v in X_t_ij:
					if k == 1 :
						l1.appned(v)
				quicksum(l1)


'''

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
            # print("c_2 is done ")
            # energy requirement constarint

            c_3 = {}
            energy_from_each_link = []
            for eh_node in V_E:
                for l in Links:
                    one_link_energy_in_whole_frame = alpha * Gn_dic[(l[0], eh_node)] *\
                        quicksum([P_t_ij[(t, l[0], l[1])] for t in Time])
                    energy_from_each_link.append(
                        one_link_energy_in_whole_frame)
                c_3[eh_node] = m.addConstr(
                    quicksum(energy_from_each_link), GRB.GREATER_EQUAL, energy_t * 1000.0)
                energy_from_each_link = []
            # print("c_3 is done ")
            #   one link activated, there is a time slot
            c_4 = {}
            for t in Time:
                for l in Links:
                    c_4[(t, l[0], l[1])] = m.addConstr(
                        Z_t[t], GRB.GREATER_EQUAL, X_t_ij[t, l[0], l[1]])
            # print("c_4 is done ")
            #  link power constraint for P_max
            c_5 = {}
            for tl in Time_Links:
                c_5[tl] = m.addConstr(
                    P_t_ij[tl], GRB.LESS_EQUAL, P_max * X_t_ij[tl])
            # print("c_5 is done ")
            # lower bond
            c_5_scd = {}
            for tl in Time_Data_link:
                lower_bond = gamma * N_0 / Gn_dic[(tl[1], tl[2])]
                c_5_scd[tl] = m.addConstr(
                    # SNR_power[(tl[1],tl[2])]*X_t_ij[tl], GRB.LESS_EQUAL, P_t_ij[tl])
                    P_min * X_t_ij[tl], GRB.LESS_EQUAL, P_t_ij[tl])
            # print("c_5_scd is done ")
            # ________each transmiter can only transmit data link or energy link in one time slot_
            c_6 = {}
            for t in Time:
                for transmitter in Tx:
                    fromx = sum(X_t_ij.select(t, transmitter, "*"))
                    c_6[(t, transmitter)] = m.addConstr(
                        fromx, GRB.LESS_EQUAL, 1)
            # print("c_6 is done")
            # ____________maximum energy time slot can given

            c_8 = []
            c_8.append(m.addConstr(quicksum(Z_t[t]
                                            for t in Time), GRB.LESS_EQUAL, T))
            # print("c_8 is done")
            # _______________ set the objective__________
            obj = []
            m.setObjective(quicksum(Z_t[t] for t in Time), GRB.MINIMIZE)
            sum_z = sum(Z_t[t] for t in Time)
            obj.append(sum_z)
          #  print("objective is done ")

            # model adjustment
            m.update()  # update when calculating
           # print("model updated ")

            #  write the optimization file
            m.write('main_1.lp')

           # print("strat_ optimization ")
            m.optimize()
            if m.status == 2:
                ave_list.append(m.objVal)
                X_t_ij_copy = {}
                schedule = {}
                P_dic = {}
                X_dic = {}
                P_e_dic = {}
                X_e_dic = {}
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
                    # here  it has to be int because the results is int
                    if round(Z_t[tl].X) > 0.0:
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

                            pure_energy_link_slot.append(
                                tuple((j[0], j[1], j[2])))
                    elif all(three != 'refpoint' for three in third_position):
                        cc2 += 1
                        for j in i:
                            pure_data_link_slot.append(
                                tuple((j[0], j[1], j[2])))
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

    # calculate energy for each type of slots
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
                        n_i_co.append(
                            P_t_ij[lk].X * Gn_dic[(lk[1], i)] * alpha)
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

            elapsedtime = timeit.default_timer() - start_time
        txt1 = 'MILPtime'+str(Number_of_Links)
        document1 = open(str(txt1), 'w')
        document1.write('time='+str(elapsedtime))
        document1.write('number of links='+str(Number_of_Links))
        document1.close()
        if len(ave_list) != 0:
            link_level_schedule_length[Given] = sum(
                ave_list) / float(len(ave_list))
            max_dic[Given] = max(ave_list)
            min_dic[Given] = min(ave_list)
        ave_list = []

        if len(number_pure_data_link_slot) != 0.0:
            link_level_data_length[Given] = sum(
                number_pure_data_link_slot) / float(len(number_pure_data_link_slot))
            max_d_dic[Given] = max(number_pure_data_link_slot)
            min_d_dic[Given] = min(number_pure_data_link_slot)

        else:
            link_level_data_length[Given] = 0.0
        number_pure_data_link_slot = []

        if len(number_pure_energy_link_slot) != 0.0:
            # file.write('energy_slot'+str(number_pure_energy_link_slot))
            link_level_energy_length[Given] = sum(
                number_pure_energy_link_slot) / float(len(number_pure_energy_link_slot))
            max_e_dic[Given] = max(number_pure_energy_link_slot)
            min_e_dic[Given] = min(number_pure_energy_link_slot)
        else:
            link_level_energy_length[Given] = 0.0
        number_pure_energy_link_slot = []

        if len(number_co_use_slot) != 0.0:
            link_level_co_length[Given] = sum(
                number_co_use_slot) / float(len(number_co_use_slot))
            max_co_dic[Given] = max(number_co_use_slot)
            min_co_dic[Given] = min(number_co_use_slot)
        else:
            link_level_co_length[Given] = 0.0

        number_co_use_slot = []

        if len(J_e_list) != 0.0:
            link_level_J_energy[Given] = sum(J_e_list) / float(len(J_e_list))
        else:
            link_level_J_energy[Given] = 0.0
        if len(J_d_list) != 0.0:
            link_level_J_data[Given] = sum(J_d_list) / float(len(J_d_list))
        else:
            link_level_J_data[Given] = 0.0

        if len(J_co_list) != 0:
            link_level_J_co[Given] = sum(J_co_list) / float(len(J_co_list))
        else:
            link_level_J_co[Given] = 0.0

        link_level_n_e[Given] = float(sum(i_e_n) / loop_time)
        link_level_n_d[Given] = float(sum(i_d_n) / loop_time)
        link_level_n_co[Given] = float(sum(i_co_n) / loop_time)
        link_level_n_propo_d[Given] = float(sum(i_eff_id_n) / loop_time)
        link_level_n_propo_e[Given] = float(sum(i_eff_ie_n) / loop_time)
        link_level_n_propo_co[Given] = float(sum(i_eff_co_n) / loop_time)
        J_e_list = []
        J_co_list = []
        J_d_list = []
        i_e_n = []
        i_d_n = []
        i_co_n = []
        i_eff_ie_n = []
        i_eff_id_n = []
        i_eff_co_n = []

    efficienty_ratio_e = {}
    efficienty_ratio_d = {}
    efficienty_ratio_co = {}

    for i in link_level_schedule_length:
        if link_level_J_energy[i] != 0:
            efficienty_ratio_e[i] = float(
                link_level_n_e[i] / link_level_J_energy[i])
        else:
            efficienty_ratio_e[i] = 0
        if link_level_J_co[i] != 0:
            efficienty_ratio_co[i] = float(
                link_level_n_co[i] / link_level_J_co[i])
        else:
            efficienty_ratio_co[i] = 0
        if link_level_J_data[i] != 0:
            efficienty_ratio_d[i] = float(
                link_level_n_d[i] / link_level_J_data[i])
        else:
            efficienty_ratio_d[i] = 0

    for i in link_level_schedule_length:
        schedule_ratio_e[i] = link_level_energy_length[i] / \
            float(link_level_schedule_length[i])
        schedule_ratio_d[i] = link_level_data_length[i] / \
            float(link_level_schedule_length[i])
        schedule_ratio_co[i] = link_level_co_length[i] / \
            float(link_level_schedule_length[i])

    denomitor = {}
    for i in link_level_J_data:
        denomitor[i] = link_level_J_data[i] + \
            link_level_J_co[i] + link_level_J_energy[i]
    for i in link_level_J_data:
        if denomitor[i] != 0:
            energy_ratio_e[i] = link_level_J_energy[i] / denomitor[i]
            energy_ratio_d[i] = link_level_J_data[i] / denomitor[i]
            energy_ratio_co[i] = link_level_J_co[i] / denomitor[i]

    elapsed = timeit.default_timer() - start_time
    txt = 'MRnew'+'SNR-1000_noenergy_requried_EHNODE'+'SNRpower' + 'gamma_'+str(gamma)+'J'+'_'+str(loop_time)+'_eh_'+str(Number_of_EHnodes)+'_size='+str(networkfieldsize)+'_link='\
          + str(start_number_of_links)+'_'+str(end_number_of_links)+'.txt'

    document = open(str(txt), 'w')
    document.write('MRnew'+'750one_noenergy_requried_EHNODE'+'SNRpower'+'N_0 = '+str(N_0)+'mW'+'\n'+'gamma = '+str(gamma)+'\n'+'Number_of_EHnodes = '+str(Number_of_EHnodes) + '\n'+'networkfieldsize = ' + str(networkfieldsize)+' m^2'+'\n' +
                   'P_max = '+str(P_max)+'mW'+'\n'+'P_min = '+str(P_min)+'mW'+'\n'+'E_min ='+str(energy_t)+'J'+'\n'+'samples = '+str(loop_time)+'\n')
    document.write('\n')
    document.write('link_level_schedule_length = ' +
                   str(link_level_schedule_length))
    document.write('\n')
    document.write('link_level_data_length = '+str(link_level_data_length))
    document.write('\n')
    document.write('link_level_co_length = '+str(link_level_co_length))
    document.write('\n')
    document.write('link_level_energy_length= '+str(link_level_energy_length))
    document.write('\n')
    document.write('max_dic = '+str(max_dic))
    document.write('\n')
    document.write('min_dic = '+str(min_dic))
    document.write('\n')
    document.write('max_d_dic = ' + str(max_d_dic))
    document.write('\n')
    document.write('min_d_dic = ' + str(min_d_dic))
    document.write('\n')
    document.write('max_e_dic = '+str(max_e_dic))
    document.write('\n')
    document.write('min_e_dic = '+str(min_e_dic))
    document.write('\n')
    document.write('max_co_dic = '+str(max_co_dic))
    document.write('\n')
    document.write('min_co_dic = '+str(min_co_dic))
    document.write('\n')
    document.write('link_level_J_energy = '+str(link_level_J_energy))
    document.write('\n')
    document.write('link_level_J_data = '+str(link_level_J_data))
    document.write('\n')
    document.write('link_level_J_co = '+str(link_level_J_co))
    document.write('\n')
    document.write('link_level_n_d = '+str(link_level_n_d))
    document.write('\n')
    document.write('link_level_n_e = '+str(link_level_n_e))
    document.write('\n')
    document.write('link_level_n_co ='+str(link_level_n_co))
    document.write('\n')
    document.write('link_level_n_propo_d = '+str(link_level_n_propo_d))
    document.write('\n')
    document.write('link_level_n_propo_e = '+str(link_level_n_propo_e))
    document.write('\n')
    document.write('link_level_n_propo_co ='+str(link_level_n_propo_co))
    document.write('\n')
    document.write('efficienty_ratio_d =' + str(efficienty_ratio_d))
    document.write('\n')
    document.write('efficienty_ratio_co = '+str(efficienty_ratio_co))
    document.write('\n')
    document.write('efficienty_ratio_e =' + str(efficienty_ratio_e))
    document.write('\n')
    document.write('schedule_ratio_co= '+str(schedule_ratio_co))
    document.write('\n')
    document.write('schedule_ratio_d ='+str(schedule_ratio_d))
    document.write('\n')
    document.write('schedule_ratio_e = '+str(schedule_ratio_e))
    document.write('\n')
    document.write("energy_ratio_d = "+str(energy_ratio_d))
    document.write('\n')
    document.write('energy_ratio_co ='+str(energy_ratio_co))
    document.write('\n')
    document.write('energy_ratio_e = '+str(energy_ratio_e))
    document.write('\n')
    document.write('time = '+str(elapsed)+' s')
    document.write('\n')
    document.write('f_counter=' + str(f_counter))
    document.write('\n')
    document.write('c_counter= '+str(c_counter))
    document.close()
    print('milp'+str(link_level_schedule_length))
    print('d='+str(link_level_data_length))
    print('e='+str(link_level_energy_length))
    print('co='+str(link_level_co_length))
    # print(schedule)
    # print(link_level_n_d[Given]+link_level_n_co[Given]+link_level_n_e[Given])
    # print(Gn_dic)
    # print(P_dic)
    # print(P_e_dic)
    # fromaddr = "yw955@uowmail.edu.au"
    # toaddr = "yw955@uowmail.edu.au"
    # msg = MIMEMultipart()
    # msg['From'] = fromaddr
    # msg['To'] = toaddr
    # msg['Subject'] = str(txt)
    # body = "from lab computer 2"
    # msg.attach(MIMEText(body, 'plain'))
    # filename = str(txt)
    # attachment = open(str(txt),'r')
    # part = MIMEBase('application', 'octet-stream')
    # part.set_payload((attachment).read())
    # encoders.encode_base64(part)
    # part.add_header('Content-Disposition', "attachment; filename= %s" % filename)
    # msg.attach(part)
    # server = smtplib.SMTP('pod51008.outlook.com', 587)
    # server.starttls()
    # server.login(fromaddr, "asdASD18")
    # text = msg.as_string()
    # server.sendmail(fromaddr, toaddr, text)
    # server.quit()
    # print('message'+str(energy_t)+'_sent')
