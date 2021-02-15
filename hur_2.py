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
import timeit
import smtplib
from functions import topology
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
from email.MIMEBase import MIMEBase
from email import encoders
from operator import itemgetter
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
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
from email.MIMEBase import MIMEBase
from email import encoders
import timeit
import numpy as np
from functions import topology
random.seed(1)




# 1. Heurisitc
# Parameter setting

  # start to record time
# alpha = 0.85  # 2. Harvesting efficiency cited from other paper, do not change
# gamma = float(100.0)  # 3. the SINR threshold, equal to 20 dB
# 4. Thermal noise. As we do not have a specific platlfrom hence averagely
# is -90 dBm which is 10 **-9 mW
N_0 = float(1 * 10**-9)
energy_t = float(10.0)  # 5. The energy requirment of EH nodes J
P_max = float(1000.0)  # 6. Maximum power in milliwatt
P_min = float(750.0)
# The following are network specifc parameters
loop_time = 1
Number_of_EHnodes = 1

networkfieldsize = 80.0  # meters This set the topology size
start_number_of_links = 15
end_number_of_links = 16
alpha = 1.0
step = 1
n_loop = [1]
alacounter =0
# start_time = timeit.default_timer()
# 2. NETWORK SET UP
Gn_dic = {}
distant_dictionary = {}
gamma = 10000
P_snr = {}
E_dic = {}
energy_transmitter_set = set()
ave_schedule = {}
ave_data_time_slot = {}
ave_co_time_slot = {}
ave_energy_time_slot = {}
n_harvested = {}
n_co = {}
n_data = {}
n_energy = {}
fn_max_dic = {}
fn_min_dic = {}
for nl in n_loop:
  	Number_of_EHnodes = nl

  	for Given in range(start_number_of_links, end_number_of_links, step):
	    
	    Number_of_Links = Given
	    start_time = timeit.default_timer() 
	    V_E = ['n_' + str(EHnodes)for EHnodes in range(Number_of_EHnodes)]
	    for n in V_E:
	      E_dic[n] = energy_t * 1000  # change the unit to mJ
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

	    for i in Energy_Links:
	      energy_transmitter_set.add(i[0])
	    Links = Data_Links_to_schedule + Energy_Links
	    ave_schedule_list = []  # this point is to clean the list cell for store schedule /
	    ave_data_slot_list = []
	    ave_energy_slot_list = []
	    ave_co_exist_list = []
	    ave_n_energy_list = []
	    ave_n_data_list = []
	    ave_n_harvested_list = []
	    ave_n_co_list = []
	    energy_slot_number = []
	    data_slot_number = []
	    co_slot_number = []
	    fn_slot_number = []
	  ####################################################################
	    for loop in range(loop_time):  # start the inner loop maning different position of links
	      schedule = {}  # create to record all the when energy links are added in
	      power_schedule = {}
	      Data_link_power_dic = {}
	      transmission_set = {}
	      # print('current_Number_of_Links=' + str(Given))
	      # print('loop=' + str(loop))
	      T_counter = 1  # slot counter for later use
	      Gn_dic, G_max, G_min, SNR_power, Phi, distant_dictionary = topology(networkfieldsize, V_E, Number_of_EHnodes,
	                                                                          Tx, Data_Links_to_schedule, gamma, N_0, Number_of_Links,
                                                                          Gn_dic, distant_dictionary, P_max)

	      




	      for i in Data_Links_to_schedule:
	        P_snr[i] = gamma * N_0 / Gn_dic[i]

	      potential_list = []
	      P_snr_copy = copy.deepcopy(P_snr)
	      # print('topology generated')
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
	        alacounter +=1

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

	      T = len(transmission_set)
	      for n in V_E:
	        E_dic[n] = energy_t * 1000
	      E_dic_copy = copy.deepcopy(E_dic)
	      # print('before = ' + str(E_dic))

	      for i in Data_link_power_dic:
	        for n in V_E:
	          E_dic_copy[n] -= Gn_dic[(i[1], n)] * \
	              Data_link_power_dic[i] * alpha

	      # print('E_dic_copy_before=' + str(E_dic_copy))

	      if any(v > 0 for v in E_dic_copy.values()):
	        # print('have to insert energy links ')

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
	            LP2C3[i] = LP2.addConstr(D_P_xy[i], GRB.LESS_EQUAL, P_max)
	          LP2C4 = {}  # power constratin for energy link
	          LP2C5 = {}
	          for i in selected_energy_list:
	            LP2C4[i] = LP2.addConstr(0, GRB.LESS_EQUAL, E_P_xy[i])
	            LP2C5[i] = LP2.addConstr(E_P_xy[i], GRB.LESS_EQUAL, P_max)

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

	            schedule[each_slot] = list(
	                D_P_xy.keys() + chosen_energy_list)

	            for i in selected_energy_list:
	              if E_P_xy[i].X > 0:

	                power_schedule[(each_slot, i[0], i[1])
	                               ] = E_P_xy[i].X
	            for i in D_P_xy:
	              power_schedule[(each_slot, i[0], i[1])] = D_P_xy[i].X
	          # print('\n')

	        T_DATA_COUNTER = T_counter

	        for lk in power_schedule:
	          for cn in V_E:
	            E_dic[cn] -= Gn_dic[(lk[1], cn)] * \
	                power_schedule[lk] * alpha
	        print(E_dic)
	        # print('after insert energy links =' + str(E_dic))
	        while any(ele > 0 for ele in E_dic.values()):
	          # print('have to insert E slot')
	          T_counter += 1
	          schedule[T_counter] = ' all energy links '
	          power_schedule[(T_counter, 'all_energy_link')
	                         ] = Number_of_Links * P_max
	          for n in V_E:
	            for elk in Energy_Links:
	              E_dic[n] -= Gn_dic[(elk[0], n)] * P_max * alpha

	        fn_schedule = schedule
	        fn_power = power_schedule

	      else:
	        # print('only data links are enough ')
	        T_DATA_COUNTER = T_counter
	        fn_schedule = transmission_set
	        fn_power = Data_link_power_dic
	      # print('after add e time slot =' + str(E_dic))
	    

	      co_counter = 0
	      d_counter = 0
	      # e_counter =0
	      co_name_list = {}

	      co_slot = []

	      data_slot = []

	      energy_slot = []

	      # print(fn_schedule)
	      for i in range(1, T_DATA_COUNTER + 1):
	        print(i)
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

	      # Noc_slot = []
	      # for cpl in 
	      co_counter = len(co_name_list)

	      co_slot_number.append(co_counter)
	      

	      # d_counter = len(data_slot)
	      Nod_slot = []

	  #     if data_slot !=0 :
			# for dlp in data_slot:	      	
			# 	Nod_slot.append(dlp[0])	      
	  # 	 	data_slot_number.append(max(Nod_slot))
   #        else:
			# data_slot_number.append(0)
	      data_slot_number.append(len(data_slot))
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

	      # fn_schedule = {}
	      # fn_power = {}
	     
       #    document1 = open(str(txt1),'w')
	      # document1.write('lpaatime='+str(elapsedt)+'NOL'+str(Number_of_Links))
	   	  
	      # document1.close()
	  	
	  	# elapsedt = timeit.default_timer() - start_time
	   #  txt1= 'lpaatime'+'number_of_links'+str(Number_of_Links)
	   #  documentl =open(str(txt1),'w')
	   #  documentl.write('lpaatime='+str(elapsedt)+'nol'+str(Number_of_Links))
	   #  documentl.close()
	   	# elapsedt = timeit.default_timer() - start_time
	    ave_schedule[Given] = sum(fn_slot_number) / float(loop_time)
	    elapsedt = timeit.default_timer() - start_time
	    txtlp='lpaatime'+str(Number_of_Links)+'.txt'
	    documentlp=open(str(txtlp),'w')
	    documentlp.write('time='+str(elapsedt)+'nol'+str(Number_of_Links))
	    documentlp.close()
	    fn_max_dic[Given] = max(fn_slot_number)
	    fn_min_dic[Given] = min(fn_slot_number)
	    ave_data_time_slot[Given] = sum(data_slot_number) / float(loop_time)
	    ave_co_time_slot[Given] = sum(co_slot_number) / float(loop_time)
	    ave_energy_time_slot[Given] = sum(
	        energy_slot_number) / float(loop_time)
	    n_energy[Given] = sum(ave_n_energy_list) / float(loop_time)
	    n_data[Given] = sum(ave_n_data_list) / float(loop_time)
	    n_co[Given] = sum(ave_n_co_list) / float(loop_time)
	    n_harvested[Given] = sum(ave_n_harvested_list) / float(loop_time)
	    elapsed = timeit.default_timer() - start_time
	  # print(fn_schedule)
	    print('LPAA=' + str(ave_schedule))
	    # print(fn_power)
	    # print(n_harvested)
	    print('hd=' + str(ave_data_time_slot))
	    print('he=' + str(ave_energy_time_slot))
	    print('hco=' + str(ave_co_time_slot))
	    txt = '100LPAAnpower_500-1000'+'energy_t'+str(energy_t)+'gamma'+str(gamma)+'_'+str(loop_time)+'_eh_'+str(Number_of_EHnodes)+'_size='+str(networkfieldsize)+'_link='\
	        +str(start_number_of_links)+'_'+str(end_number_of_links)+'.txt'

	    document = open(str(txt),'w')
	    document.write('LPAA'+'N_0 = '+str(N_0)+'mW'+'\n'+'gamma = '+str(gamma)+'\n'+'Number_of_EHnodes = '+str(Number_of_EHnodes) + '\n'+'networkfieldsize = '+ str(networkfieldsize)+' m^2'+'\n'+
	    'P_max = '+str(P_max)+'mW'+'\n'+'P_min = '+str(P_min)+'mW'+'\n'+'E_min ='+str(energy_t)+'J'+'\n'+'samples = '+str(loop_time)+'\n')
	    document.write('\n')
	    document.write('ave_schedule = '+str(ave_schedule))
	    document.write('\n')
	    document.write('fn_max_dic = '+str(fn_max_dic))
	    document.write('\n')
	    document.write('fn_min_dic = '+str(fn_min_dic))
	    document.write('\n')
	    document.write('ave_data_time_slot= '+str(ave_data_time_slot))
	    document.write('\n')
	    document.write('ave_co_time_slot = '+str(ave_co_time_slot))
	    document.write('\n')
	    document.write('ave_energy_time_slot'+str(ave_energy_time_slot))
	    document.write('\n')
	    document.write('n_energy =' +str(n_energy))
	    document.write('\n')
	    document.write('n_data = ' +str(n_data))
	    document.write('\n')
	    document.write('n_co = '+str(n_co))
	    document.write('\n')
	    document.write('n_harvested = '+str(n_harvested))
	    document.write('\n')
	    # document.write('max_co_dic = '+str(max_co_dic))
	    # document.write('\n')
	    # document.write('min_co_dic = '+str(min_co_dic))
	    # document.write('\n')
	    # document.write('link_level_J_energy = '+str(link_level_J_energy))
	    # document.write('\n')
	    # document.write('link_level_J_data = '+str(link_level_J_data))
	    # document.write('\n')
	    # document.write('link_level_J_co = '+str(link_level_J_co))
	    # document.write('\n')
	    # document.write('link_level_n_d = '+str(link_level_n_d))
	    # document.write('\n')
	    # document.write('link_level_n_e = '+str(link_level_n_e))
	    # document.write('\n')
	    # document.write('link_level_n_co ='+str(link_level_n_co))
	    # document.write('\n')
	    # document.write('link_level_n_propo_d = '+str(link_level_n_propo_d))
	    # document.write('\n')
    # document.write('link_level_n_propo_e = '+str(link_level_n_propo_e))
    # document.write('\n')
    # document.write('link_level_n_propo_co ='+str(link_level_n_propo_co))
    # document.write('\n')
    # document.write('efficienty_ratio_d =' +str(efficienty_ratio_d))
    # document.write('\n')
    # document.write('efficienty_ratio_co = '+str(efficienty_ratio_co))
    # document.write('\n')
    # document.write('efficienty_ratio_e =' +str(efficienty_ratio_e))
    # document.write('\n')
    # document.write('schedule_ratio_co= '+str(schedule_ratio_co))
    # document.write('\n')
    # document.write('schedule_ratio_d ='+str(schedule_ratio_d))
    # document.write('\n')
    # document.write('schedule_ratio_e = '+str(schedule_ratio_e))
    # document.write('\n')
    # document.write("energy_ratio_d = "+str(energy_ratio_d))
    # document.write('\n')
    # document.write('energy_ratio_co ='+str(energy_ratio_co))
    # document.write('\n')
    # document.write('energy_ratio_e = '+str(energy_ratio_e))
    # document.write('\n')
	    document.write('time = '+str(elapsed)+' s' )
	    document.write('\n')
	    # document.write('f_counter=' +str(f_counter))
	    # document.write('\n')
	    # document.write('c_counter= '+str(c_counter))
	    document.close()
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
