import random
import os
from itertools import product, combinations

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

from scipy import stats
from astropy.stats import kuiper, kuiper_two



###################
# JUNCTION ANGLES #
###################

def junction_angles(filepath_to_junction_angles = "./basins/PART2/all_junction_angles.csv", destination_directory = "./junction_angles"):
    """
    This function coordinates the analysis of the junction angles part of the analysis. 
    The result is a figure with CDF's of junction angles of each basin, CDF's of Area Ratios of each basin, and a boxplot with Kuiper's V-statistic for the groups VRP-VRP, AP-AP and AP-VRP
   
    11/01/2024
    """
    
    # 1. INITIAL STUFF
    # make directory
    if os.path.exists(destination_directory) == False:
        os.mkdir(destination_directory)
        
    # open junction angles csv
    ja_df = pd.read_csv(filepath_to_junction_angles)
    
    # 2. CALCULATE KUIPER'S V STATISTICS
    # Initiate threee lists to store Kuiper's V statistics: one for each angle of the junction
    # Each list gets three sublists: One for V statistics between basins in the VRP, one for V statistics between basins in the AP, and one for V stastistics between AP and VRP basins
    A_Vstats  = [[],[],[]]
    B1_Vstats = [[],[],[]]
    B2_Vstats = [[],[],[]]
    
    # get all new_basin_keys for the populations of VRP and AP (the new_basin_key is a unique identifier for each basin)
    VRP_mask = ja_df["Province"] == "VRP"
    VRP_basin_keys = ja_df.loc[VRP_mask, "new_basin_key"].unique().tolist()
    
    AP_mask = ja_df["Province"] == "AP"
    AP_basin_keys = ja_df.loc[AP_mask, "new_basin_key"].unique().tolist()
    print("The number of basins in the VRP and in the AP are:", len(VRP_basin_keys), len(AP_basin_keys))
    
    # find all unique combinations of basins
    VRP_VRP_combinations = list(combinations(VRP_basin_keys, 2))
    AP_AP_combinations = list(combinations(AP_basin_keys, 2))
    VRP_AP_combinations = list(product(VRP_basin_keys, AP_basin_keys))
    print("The number of unique combinations for VRP-VRP, AP-AP and VRP-AP are:", len(VRP_VRP_combinations), len(AP_AP_combinations), len(VRP_AP_combinations))
    
    # add all these combinations into one list of 'pairs' of basins for which the V statistics needs to be calculated
    basin_pairs = VRP_VRP_combinations + AP_AP_combinations + VRP_AP_combinations
    
    # now loop through all the basin pairs, extract the data from the junction angle dataframe, calculate the V statistic and store this in the right list (A_Vstats, B1_Vstat, B2_Vstats) 
    for basin_pair in basin_pairs:
        # get the new_basin_keys of the basin_pairs (X denotes one basin, Y the other)
        X_new_basin_key = basin_pair[0]
        Y_new_basin_key = basin_pair[1]
        
        # get the junction angles (in list form) for the respective basins
        X_A = ja_df[ja_df["new_basin_key"] == X_new_basin_key]["donors_junction_angle"].tolist()
        X_BENDING1 = ja_df[ja_df["new_basin_key"] == X_new_basin_key]["donor1_receiver_junction_angle"].tolist()
        X_BENDING2 = ja_df[ja_df["new_basin_key"] == X_new_basin_key]["donor2_receiver_junction_angle"].tolist()
        Y_A = ja_df[ja_df["new_basin_key"] == Y_new_basin_key]["donors_junction_angle"].tolist()
        Y_BENDING1 = ja_df[ja_df["new_basin_key"] == Y_new_basin_key]["donor1_receiver_junction_angle"].tolist()
        Y_BENDING2 = ja_df[ja_df["new_basin_key"] == Y_new_basin_key]["donor2_receiver_junction_angle"].tolist()
        
        # the bending angles BENDING1 and BENDING2 do not conform yet to (my) definiton of bending angles, where B1 is greater than B2.
        # the following couple lines ind B1 and B2 for both X and Y based on the Bending angles extracted from the dataframe
        X_B1 = [max(val1, val2) for val1, val2 in zip(X_BENDING1, X_BENDING2)]
        X_B2 = [min(val1, val2) for val1, val2 in zip(X_BENDING1, X_BENDING2)]
        
        Y_B1 = [max(val1, val2) for val1, val2 in zip(Y_BENDING1, Y_BENDING2)]
        Y_B2 = [min(val1, val2) for val1, val2 in zip(Y_BENDING1, Y_BENDING2)]
        
        """        
        # now ready to perform kuiper tests! Using kuiper_two function from Astropy
        A_kuiperresult = kuiper_two(X_A, Y_A)
        B1_kuiperresult = kuiper_two(X_B1, Y_B1)
        B2_kuiperresult = kuiper_two(X_B2, Y_B2)
        
        # get the V statistic from the kuiperresult tuple
        A_V = A_kuiperresult[0]
        B1_V = B1_kuiperresult[0]
        B2_V = B2_kuiperresult[0]
        """
    
        A_V = calculate_kuipers_V(X_A, Y_A)
        B1_V = calculate_kuipers_V(X_B1, Y_B1)
        B2_V = calculate_kuipers_V(X_B2, Y_B2)
        
        # append the result to the right list
        index = 0 if basin_pair in VRP_VRP_combinations else 1 if basin_pair in AP_AP_combinations else 2
        A_Vstats[index].append(A_V)
        B1_Vstats[index].append(B1_V)
        B2_Vstats[index].append(B2_V)

        
    # 3. GET THE JUNCTION ANGLES and AREA RAIOS OF ALL BASINS IN THE VRP AND OF ALL BASINS IN THE VRP
    # set up lists
    VRP_A = []
    VRP_B1 = []
    VRP_B2 = []
    AP_A = []
    AP_B1 = []
    AP_B2 = []
    
    VRP_AR = []
    AP_AR = []
    
    all_basin_keys = AP_basin_keys + VRP_basin_keys
    for basin_key in all_basin_keys:
        A = ja_df[ja_df["new_basin_key"] == basin_key]["donors_junction_angle"].tolist()
        BENDING1 = ja_df[ja_df["new_basin_key"] == basin_key]["donor1_receiver_junction_angle"].tolist()
        BENDING2 = ja_df[ja_df["new_basin_key"] == basin_key]["donor2_receiver_junction_angle"].tolist()
        donor1drainagearea = ja_df[ja_df["new_basin_key"] == basin_key]["donor1_drainage_area"].tolist()
        donor2drainagearea = ja_df[ja_df["new_basin_key"] == basin_key]["donor2_drainage_area"].tolist()
        
        B1 = [max(val1, val2) for val1, val2 in zip(BENDING1, BENDING2)]
        B2 = [min(val1, val2) for val1, val2 in zip(BENDING1, BENDING2)]
        AR = [min(val1, val2) / max(val1, val2) for val1, val2 in zip(donor1drainagearea, donor2drainagearea)]
        
        if basin_key in VRP_basin_keys:
            VRP_A.append(A)
            VRP_B1.append(B1)
            VRP_B2.append(B2)
            VRP_AR.append(AR)
        if basin_key in AP_basin_keys:
            AP_A.append(A)
            AP_B1.append(B1)
            AP_B2.append(B2)
            AP_AR.append(AR)
       
    # 4. PLOT THE FIGURE (contains BOXPLOTS, CDFs of JUNCTION ANGLES and CDFs OF AREA RATIOS)
    # set up figure with gridspec
    fig = plt.figure(figsize=(7,8))
    
    gs_outer = gridspec.GridSpec(2,1, height_ratios = [1,1])
    gs_upper = gridspec.GridSpecFromSubplotSpec(1,1, subplot_spec=gs_outer[0])
    gs_lower = gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec=gs_outer[1], width_ratios = [1,1], wspace = 0.3) # first gridspec position for the boxplot, second for AR cdf
    
    # create the boxplot (using nested gridspec)
    gs_bplot = gridspec.GridSpecFromSubplotSpec(1,3, subplot_spec=gs_lower[0], wspace=0)
    ax1 = plt.subplot(gs_bplot[0])
    ax2 = plt.subplot(gs_bplot[1])
    ax3 = plt.subplot(gs_bplot[2])
    
    bp1 = ax1.boxplot(A_Vstats, patch_artist=True, medianprops = dict(color = "black"))
    bp2 = ax2.boxplot(B2_Vstats, patch_artist=True, medianprops = dict(color = "black"))
    bp3 = ax3.boxplot(B1_Vstats, patch_artist=True, medianprops = dict(color = "black"))
    
  
    colours = ["darkorange", "limegreen", "limegreen"]
    for bp in (bp1, bp2, bp3):
        for i, (patch, colour) in enumerate(zip(bp["boxes"], colours)):
            patch.set_facecolor(colour)
            if i == 2: 
                patch.set_edgecolor("darkorange")
        
    
    ax2.set_yticks([])
    ax3.set_yticks([])
    ax2.set_ylabel("")
    ax3.set_ylabel("")
    ax1.spines["right"].set_visible(False)
    ax2.spines["left"].set_color("grey")
    ax2.spines["right"].set_visible(False)
    ax3.spines["left"].set_color("grey")
       
    ax1.set_xticks([])
    ax2.set_xticks([])  
    ax3.set_xticks([])  
        
    ax1.set_xlabel(r"$\alpha$")
    ax2.set_xlabel(r"$\beta_2$")
    ax3.set_xlabel(r"$\beta_1$")
    
    ax1.set_ylabel("V")
    
    limits_axes = (0, max(ax1.get_ylim()[1],ax2.get_ylim()[1],ax3.get_ylim()[1]))
    ax1.set_ylim(limits_axes)
    ax2.set_ylim(limits_axes)
    ax3.set_ylim(limits_axes)
                      
    # create subplot for junction angle CDFs 
    ax4 = plt.subplot(gs_upper[0])
    
    VRP_all_angles = [VRP_A, VRP_B1, VRP_B2]
    for list_of_angle_lists in VRP_all_angles:
        for angle_list in list_of_angle_lists:
            ax4.plot(np.sort(angle_list), np.linspace(0,1,len(angle_list), endpoint=True), color = "darkorange", linewidth = 0.2, alpha = 0.5)
  
    AP_all_angles = [AP_A, AP_B1, AP_B2]
    for list_of_angle_lists in AP_all_angles:
        for angle_list in list_of_angle_lists:
            ax4.plot(np.sort(angle_list), np.linspace(0,1,len(angle_list), endpoint=True), color = "limegreen", linewidth = 0.2, alpha = 0.5)
    
    for lst in [VRP_A, VRP_B1, VRP_B2]:
        angle_list = [value for sublist in lst for value in sublist] # flattens list
        ax4.plot(np.sort(angle_list), np.linspace(0,1,len(angle_list), endpoint=True), color = "#ff6200", linewidth = 1, linestyle = "--")
    for lst in [AP_A, AP_B1, AP_B2]:
        angle_list = [value for sublist in lst for value in sublist] # flattens list
        ax4.plot(np.sort(angle_list), np.linspace(0,1,len(angle_list), endpoint=True), color = "#00b300", linewidth = 1, linestyle = "--")
    
    ax4.set_yticks([0,0.25,0.5,0.75,1])
    ax4.set_ylim([0,1])
    ax4.set_ylabel("Probability")
    
    ax4.set_xlim([0,240])
    ax4.set_xticks([0,30,60,90,120,150,180,210,240], labels = ["0","30","60","90","120","150","180","210","240"])
    ax4.set_xticks([0,10,20,40,50,70,80,100,110,130,140,160,170,190,200,220,230], minor = True)
    ax4.set_xlabel("Angle")
    
    # create subplot for area ratios CDFs
    ax5 = plt.subplot(gs_lower[1])
    for list_of_AR in VRP_AR:
        ax5.plot(np.sort(list_of_AR), np.linspace(0,1,len(list_of_AR), endpoint=True), color = "darkorange", linewidth = 0.2, alpha = 0.5)
    for list_of_AR in AP_AR:
        ax5.plot(np.sort(list_of_AR), np.linspace(0,1,len(list_of_AR), endpoint=True), color = "limegreen", linewidth = 0.2, alpha = 0.5)

    list_of_AR = [value for sublist in VRP_AR for value in sublist] # flattens list
    ax5.plot(np.sort(list_of_AR), np.linspace(0,1,len(list_of_AR), endpoint=True), color = "#ff6200", linewidth = 1, linestyle = "--")
    
    list_of_AR = [value for sublist in AP_AR for value in sublist] # flattens list
    ax5.plot(np.sort(list_of_AR), np.linspace(0,1,len(list_of_AR), endpoint=True), color = "#00b300", linewidth = 1, linestyle = "--")

    ax5.set_yticks([0,0.25,0.5,0.75,1])
    ax5.set_ylim([0,1])
    ax5.set_ylabel("Probability")
    
    ax5.set_xlim([0,1])
    ax5.set_xlabel("$A_R$")
    # adjust plot layout
    plt.tight_layout()
        

    # save figure
    plt.savefig(destination_directory+"/junctionanglesfig.jpg", format = "jpg", dpi = 360)
    
######################################################
# SEGMENT ORIENTATIONS UNIFORMITY/BIMODALITY TESTING #  
######################################################
def segment_orientations(filepath_to_segment_bearings = "./basins/PART2/all_segment_bearings.csv", destination_directory = "./segment_orientations"):
    # 1. INITIAL STUFF
    # make directory
    if os.path.exists(destination_directory) == False:
        os.mkdir(destination_directory)
        
    # open junction angles csv
    sb_df = pd.read_csv(filepath_to_segment_bearings)
    
    # remove directionality from bearings, i.e. make orientations column on [0,180)
    sb_df["segment_orientation"] = remove_directionality(sb_df["segment_bearing_whole_segment"].tolist())
    
    # get all new_basin_keys for the populations of VRP and AP (the new_basin_key is a unique identifier for each basin)
    VRP_mask = sb_df["Province"] == "VRP"
    VRP_basin_keys = sb_df.loc[VRP_mask, "new_basin_key"].unique().tolist()
    
    AP_mask = sb_df["Province"] == "AP"
    AP_basin_keys = sb_df.loc[AP_mask, "new_basin_key"].unique().tolist()
    N_basins_VRP = len(VRP_basin_keys)
    N_basins_AP = len(AP_basin_keys)
    print("The number of basins in the VRP and in the AP are:", N_basins_VRP, N_basins_AP)
    
    # 2. GENERATE DISTRIBUTIONS
    # make three lists, one for unimodal, one for bimodal distributions and one for skewed bimodal distributions
    # each list will contain 11 lists, one for every 0.1 interval of p, i.e. 0.0, 0.1, ... , 0.9, 1.0
    # these sublists will contain either 180 (unimodal) or 90 (bimodal) lists with distributions, one for each 1 degree interval of location of mode
    unimodal_distributions = []
    bimodal_distributions = []
    skewed_bimodal_distributions = []
    
    # now generate the distributions
    p_values = np.linspace(0,1,11)
    unimodal_peaks = np.arange(0,180,1)
    bimodal_peaks = np.arange(0,90,1)
    skewed_bimodal_peaks = np.arange(0,90,1)
    
    for p_value in p_values:
        unimodal_distributions_this_p = []
        bimodal_distributions_this_p = []
        skewed_bimodal_distributions_this_p = []
        for peak in unimodal_peaks:
            unimodal_distributions_this_p.append(generate_distribution(p_value, peak, "unimodal", s = 30))
        for peak in bimodal_peaks:
            bimodal_distributions_this_p.append(generate_distribution(p_value, peak, "bimodal", s = 30))
        for peak in skewed_bimodal_peaks:
            skewed_bimodal_distributions_this_p.append(generate_distribution(p_value, peak, "skewed_bimodal", s = 30))
            
        unimodal_distributions.append(unimodal_distributions_this_p)
        bimodal_distributions.append(bimodal_distributions_this_p)
        skewed_bimodal_distributions.append(skewed_bimodal_distributions_this_p)
    
    # 3. FOR VRP AND ARP BASINS, FIND LOWEST VALUE OF V FOR EACH P_VALUE FOR THE UNIMODAL AND BIMODAL DISTRIBUTIONS
    # make four lists, one for unimodal V values and one for bimodal V values, for the VRP and for the AP
    # these lists will then contain 11 lists, one for every 0.1 interval of p
    # these sublists will contain for each basin the lowest V value found in that interval
    VRP_UNIMODAL_V = []
    AP_UNIMODAL_V = []
    VRP_BIMODAL_V = []
    AP_BIMODAL_V = []
    
    """
    # first deal with unimodal distributions
    for distributions_this_p in unimodal_distributions:
        # first do VRP basins
        # initiate V values list for this p
        VRP_thisp_V = []
        for basin_key in VRP_basin_keys:
            #get segment_orientations
            segment_orientations = sb_df[sb_df["new_basin_key"]==basin_key]["segment_orientation"].tolist()
            # find minimum V value for this p for this basin
            Vvals = []
            for distribution in distributions_this_p:
                # perform kuipers test
                #Vvals.append(kuiper_two(segment_orientations, distribution)[0])
                Vvals.append(calculate_kuipers_V(segment_orientations, distribution))
            VRP_thisp_V.append(min(Vvals))
        VRP_UNIMODAL_V.append(VRP_thisp_V)
        # do the same thing for the AP
        AP_thisp_V = []
        for basin_key in AP_basin_keys:
            #get segment_orientations
            segment_orientations = sb_df[sb_df["new_basin_key"]==basin_key]["segment_orientation"].tolist()
            # find minimum V value for this p for this basin
            Vvals = []
            for distribution in distributions_this_p:
                # perform kuipers test
                #Vvals.append(kuiper_two(segment_orientations, distribution)[0])
                Vvals.append(calculate_kuipers_V(segment_orientations, distribution))
            AP_thisp_V.append(min(Vvals))
        AP_UNIMODAL_V.append(AP_thisp_V)
    """
        
    # now deal with bimodal distributions 
    for distributions_this_p in bimodal_distributions:
        # first do VRP basins
        # initiate V values list for this p
        VRP_thisp_V = []
        for basin_key in VRP_basin_keys:
            #get segment_orientations
            segment_orientations = sb_df[sb_df["new_basin_key"]==basin_key]["segment_orientation"].tolist()
            # find minimum V value for this p for this basin
            Vvals = []
            for distribution in distributions_this_p:
                # perform kuipers test
                # Vvals.append(kuiper_two(segment_orientations, distribution)[0])
                Vvals.append(calculate_kuipers_V(segment_orientations, distribution))
            VRP_thisp_V.append(min(Vvals))
        VRP_BIMODAL_V.append(VRP_thisp_V)
        # do the same thing for the AP
        AP_thisp_V = []
        for basin_key in AP_basin_keys:
            #get segment_orientations
            segment_orientations = sb_df[sb_df["new_basin_key"]==basin_key]["segment_orientation"].tolist()
            # find minimum V value for this p for this basin
            Vvals = []
            for distribution in distributions_this_p:
                # perform kuipers test
                # Vvals.append(kuiper_two(segment_orientations, distribution)[0])
                Vvals.append(calculate_kuipers_V(segment_orientations, distribution))
            AP_thisp_V.append(min(Vvals))
        AP_BIMODAL_V.append(AP_thisp_V)
    
    # 3. FOR VRP BASINS, FIND - PER STREAM ORDER - LOWEST VALUE OF V FOR EACH P VALUE OF SKEWED BIMODAL DISTRIBUTIONS
    # make a list for the VRP
    # the list will contain 7 sublists, one for each stream order
    # these sublists will then contain 11 subsublists, one for every 0.1 interval of p
    # these subsublists will contain for each basin the lowest V value found in that interval
    VRP_SKEWED_BIMODAL_V_ALL_ORDERS = []
    # keep track of how many basins n per order
    VRP_n_basins = []
    
    stream_orders = np.linspace(1,6,6) # note there are technically basins with stream order 7 but none of them have more than 10 individual segment bearings in them so don't want to consider these low samples
    
    for stream_order in stream_orders:  
        print("Now doing skewed bimodal test for stream order:", int(stream_order))
        VRP_SKEWED_BIMODAL_V = []
        for distributions_this_p in skewed_bimodal_distributions:
            # first do VRP basins
            # initiate V values list for this p
            VRP_thisp_V = []
            for basin_key in VRP_basin_keys:
                #get segment_orientations
                segment_orientations = sb_df[(sb_df["new_basin_key"]==basin_key) & (sb_df["junction_stream_order"] == int(stream_order))]["segment_orientation"].tolist()
                # only inlcude basin if it has 10 or more segment bearings for this stream order
                if len(segment_orientations) >= 10:
                    # find minimum V value for this p for this basin
                    Vvals = []
                    for distribution in distributions_this_p:
                        # perform kuipers test
                        # Vvals.append(kuiper_two(segment_orientations, distribution)[0])
                        Vvals.append(calculate_kuipers_V(segment_orientations, distribution))
                    VRP_thisp_V.append(min(Vvals))
            VRP_SKEWED_BIMODAL_V.append(VRP_thisp_V)
        VRP_n_basins.append(len(VRP_SKEWED_BIMODAL_V[0]))
        VRP_SKEWED_BIMODAL_V_ALL_ORDERS.append(VRP_SKEWED_BIMODAL_V)
        
    # 5. PLOT THE FIGURES
    # Firt plot
    # intitiate figure
    fig = plt.figure(figsize=(5,5))
    
    # set up gridspecs
    gs = gridspec.GridSpec(1,2, wspace = 0) 
    
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    
    bp1 = ax1.boxplot(VRP_BIMODAL_V, labels = p_values, patch_artist=True, medianprops = dict(color = "black"))
    bp2 = ax2.boxplot(AP_BIMODAL_V, labels = p_values, patch_artist=True, medianprops = dict(color = "black"))
    
    for box in bp1['boxes']:
        box.set(facecolor="darkorange")
    for box in bp2['boxes']:
        box.set(facecolor="limegreen")
    
    ax1.set_xlabel("$P_1$")
    
    ax1.set_xticklabels(["0.0", "", "0.2", "", "0.4", "", "0.6", "", "0.8", "", "1.0"], rotation="vertical")
    ax2.set_xticklabels(["0.0", "", "0.2", "", "0.4", "", "0.6", "", "0.8", "", "1.0"], rotation="vertical")
    
    ax1.set_ylabel("V")               
    limits_axes_bimodal = (0, max(ax1.get_ylim()[1],ax2.get_ylim()[1]))
    ax1.set_ylim(limits_axes_bimodal)
    ax2.set_ylim(limits_axes_bimodal)
    
    ax1.text(1.0, 0.005, "N = " + str(N_basins_VRP), color = "darkorange")
    ax2.text(1.0, 0.005, "N = " + str(N_basins_AP), color = "limegreen")
    
    ax2.set_yticks([])

    ax2.spines["left"].set_color("grey") 
    
    plt.savefig(destination_directory+"/segmentorientationsfig1.jpg", format = "jpg", dpi = 360)
    
    
    # Second figure
    plt.clf()
    # intitiate figure
    fig = plt.figure(figsize=(8.5,5))
    gs_skewedbimodal = gridspec.GridSpec(1,6, wspace = 0.1)
    
    for stream_order in stream_orders:
        ax5 = plt.subplot(gs_skewedbimodal[int(stream_order)-1])
        bp5 = ax5.boxplot(VRP_SKEWED_BIMODAL_V_ALL_ORDERS[int(stream_order)-1], labels = p_values , patch_artist=True, medianprops = dict(color = "black"))
        for box in bp5['boxes']:
            box.set(facecolor="darkorange")
    
        ax5.set_xticklabels(["0.0", "", "0.2", "", "0.4", "", "0.6", "", "0.8", "", "1.0"], rotation="vertical")
        
        ax5.set_ylabel("V")
        ax5.set_title("Order: " + str(int(stream_order)))
        ax5.text(0.5, 0.015, "N = " + str(VRP_n_basins[int(stream_order)-1]), color = "darkorange")
        ax5.set_ylim(0,0.7)
        
        if stream_order == 1:
            ax5.spines["right"].set_visible(False)
        if stream_order > 1:
            ax5.set_yticks([])
            ax5.set_ylabel("")
            ax5.spines["left"].set_visible(False) 
        if stream_order == 3: 
            ax5.set_xlabel("$P_2$")
    
    plt.savefig(destination_directory+"/segmentorientationsfig2.jpg", format = "jpg", dpi = 360)
    
#####################################################
# HALF-DIAL PLOT OF SEGMENT ORIENTATION AND DVALUES #
#####################################################
def halfdials(filepath_to_segmentbearings_bedddingdata = "./basins/PART1/all_segment_bearings_bedding_stats.csv", filepath_to_AP_segmentbearingsdata = "./basins/PART2/AP_MID_M/AP_MID_M_COP30_UTM_segment_bearings_filtered_nolakes.csv", filepath_to_beddingdata = "./basins/PART1/all_unique_bedding_measurements.csv", destination_directory = "./halfdials", number_of_bins = 12, include_means_d_strike_stream_orders = [0,2,3,4,5,6], include_means_d_dip_stream_orders = [0,1,2,3], figsize = (7.5,11), hspace = 0):
    # 1. INITIAL STUFF
    # make directory
    if os.path.exists(destination_directory) == False:
        os.mkdir(destination_directory)
        
    # open segment bearings bedding data csv, AP segment bearings and unique bedding data csv 
    VRP_sb_df = pd.read_csv(filepath_to_segmentbearings_bedddingdata)
    AP_sb_df = pd.read_csv(filepath_to_AP_segmentbearingsdata)
    unique_beddingdata = pd.read_csv(filepath_to_beddingdata)
    
    # turn segment bearing of AP data into segment orientation
    AP_sb_df["segment_orientation"] = remove_directionality(AP_sb_df["segment_bearing_whole_segment"].tolist())
    
    # put segment orientations of VRP and AP on [-90,90) (= [W,E)) 
    AP_sb_df["segment_orientation_WE"] = orientations_from_NS_to_WE(AP_sb_df["segment_orientation"].tolist())
    VRP_sb_df["segment_orientation_WE"] = orientations_from_NS_to_WE(VRP_sb_df["segment_orientation"].tolist())
    
    # 2. CALCULATE MEAN STRIKE AND DIP ORIENTATIONS ON [-90, 90) ( = [W, E) )
    unique_beddingdata["strike_orientation"] = unique_beddingdata.apply(lambda row: remove_directionality(row.Azimuth), axis = 1)
    unique_beddingdata["dip_orientation"] = unique_beddingdata.apply(lambda row: remove_directionality(row.Dip_Direct), axis = 1) 
    mean_strike = np.rad2deg(stats.circmean(np.deg2rad(unique_beddingdata["strike_orientation"].to_numpy()), high=np.pi/2, low=-np.pi/2))
    mean_dip = np.rad2deg(stats.circmean(np.deg2rad(unique_beddingdata["dip_orientation"].to_numpy()), high=np.pi/2, low=-np.pi/2))
    
    print(mean_strike, mean_dip)
    # 3. FOR VRP WITH BEDDING MEASUREMENTS, FIND MEANS OF difference_segment_orientation_and_mean_strike AND difference_segment_orientation_and_mean_dip_orientation
    max_stream_order = VRP_sb_df["junction_stream_order"].max()
    stream_orders = np.linspace(0, max_stream_order, max_stream_order+1).astype(int)
    
    d_so_ms_mean_perorder = []   
    d_so_md_mean_perorder = []
    for stream_order in stream_orders:
        if stream_order == 0:
            # I made stream order 0 to refer to segments of all orders
            d_so_ms = VRP_sb_df[(VRP_sb_df["difference_segment_orientation_and_mean_strike"] >= -45) & (VRP_sb_df["difference_segment_orientation_and_mean_strike"] <= 45)]["difference_segment_orientation_and_mean_strike"].to_numpy()
            d_so_md = VRP_sb_df[(VRP_sb_df["difference_segment_orientation_and_mean_dip_orientation"] >= -45) & (VRP_sb_df["difference_segment_orientation_and_mean_dip_orientation"] <= 45)]["difference_segment_orientation_and_mean_dip_orientation"].to_numpy()                       
        else:  
            d_so_ms = VRP_sb_df[(VRP_sb_df["difference_segment_orientation_and_mean_strike"] >= -45) & (VRP_sb_df["difference_segment_orientation_and_mean_strike"] <= 45) & (VRP_sb_df["junction_stream_order"] == stream_order)]["difference_segment_orientation_and_mean_strike"].to_numpy()
            d_so_md = VRP_sb_df[(VRP_sb_df["difference_segment_orientation_and_mean_dip_orientation"] >= -45) & (VRP_sb_df["difference_segment_orientation_and_mean_dip_orientation"] <= 45) & (VRP_sb_df["junction_stream_order"] == stream_order)]["difference_segment_orientation_and_mean_dip_orientation"].to_numpy()                       
        mean_d_so_ms_thisorder = np.rad2deg(stats.circmean(np.deg2rad(d_so_ms), high = np.pi/2, low=-np.pi/2))
        mean_d_so_md_thisorder = np.rad2deg(stats.circmean(np.deg2rad(d_so_md), high = np.pi/2, low=-np.pi/2))
        d_so_ms_mean_perorder.append(mean_d_so_ms_thisorder)  
        d_so_md_mean_perorder.append(mean_d_so_md_thisorder)
    
    # convert means of difference segment orientation with mean dip orientation to the reference frame of strike
    d_so_md_mean_perorder = [mean - 90 if mean >= 0 else mean + 90 for mean in d_so_md_mean_perorder]
    
    
    # 4. PLOT FIGURE
    fig = plt.figure(figsize=figsize)
    # set up gridspec
    gs_outer = gridspec.GridSpec(2,1, height_ratios = [1,21]) # vertical grid, upper one is for text labels, lower one is for all the plots
    width_ratios = [1,6,6,6,1]
    gs_upper = gridspec.GridSpecFromSubplotSpec(1,5, subplot_spec=gs_outer[0], width_ratios=width_ratios) # horizontal grid for text labels
    gs_lower = gridspec.GridSpecFromSubplotSpec(1,5, subplot_spec=gs_outer[1], width_ratios=width_ratios) # horizontal grid for columns of halfdials
    
    # set up text - headers
    ax = plt.subplot(gs_upper[1])
    y = (ax.get_ylim()[0]+ax.get_ylim()[1])/2
    x = (ax.get_xlim()[0]+ax.get_xlim()[1])/2
    ax.text(x,y, s = "Appalachian Plateau", ha='center', va='center')
    ax.set_axis_off()
 
    ax = plt.subplot(gs_upper[2])
    y = (ax.get_ylim()[0]+ax.get_ylim()[1])/2
    x = (ax.get_xlim()[0]+ax.get_xlim()[1])/2
    ax.text(x,y, s = "Valley and Ridge", ha='center', va='center')
    ax.set_axis_off()
    
    ax = plt.subplot(gs_upper[3])
    y = (ax.get_ylim()[0]+ax.get_ylim()[1])/2
    x = (ax.get_xlim()[0]+ax.get_xlim()[1])/2
    ax.text(x,y, s = "Difference between \n orientation of segment \n and strike in Valley \n and Ridge", ha='center', va='center')
    ax.set_axis_off()
                
    # plot half dials
    AP_N_list = []
    VRP_N_list = []
    
    for column in [1,2,3]: # 1 = AP, 2 = VRP, 3 = VRP_d
        # make column gridspec
        gs_column = gridspec.GridSpecFromSubplotSpec(max_stream_order+1,1, subplot_spec=gs_lower[column], hspace = hspace)
        # loop through stream orders 
        for stream_order in stream_orders:
            # make axes object
            ax = plt.subplot(gs_column[stream_order], projection = "polar")
            # get data in
            if column == 1: 
                if stream_order == 0:
                    # when plotting segments of all orders, only get AP segments with orders up to the max order in the VRP dataset
                    data = AP_sb_df[AP_sb_df["junction_stream_order"]<=max_stream_order]["segment_orientation_WE"].to_numpy()
                else: 
                    # get segment orientations West-East for this order
                    data = AP_sb_df[AP_sb_df["junction_stream_order"]==stream_order]["segment_orientation_WE"].to_numpy()
            elif column == 2:
                if stream_order == 0:
                    data = VRP_sb_df["segment_orientation_WE"].to_numpy()
                else:
                    data = VRP_sb_df[VRP_sb_df["junction_stream_order"]==stream_order]["segment_orientation_WE"].to_numpy()
            elif column == 3:
                if stream_order == 0:
                    data = VRP_sb_df["difference_segment_orientation_and_mean_strike"].to_numpy()
                else:
                    data = VRP_sb_df[VRP_sb_df["junction_stream_order"]==stream_order]["difference_segment_orientation_and_mean_strike"].to_numpy()
           
            N_segments = len(data)
            
            # append N to respective list
            if column == 1:
                AP_N_list.append(N_segments)
            elif column == 2:
                VRP_N_list.append(N_segments)
            
            # compute the theta of each bin, that is, the starting angle of each bin, based on the parameter number_of_bins
            # e.g, for number_of_bins = 4, theta_bins = [-90, -45, 0, 45] 
            binsize = 180/number_of_bins
            theta_bins = np.linspace(-90, 90-binsize, number_of_bins)
                
            # count the number of values within each bin
            count_bins = []
            for z in range(number_of_bins):
                min_theta = theta_bins[z]
                max_theta = min_theta + binsize
                count = 0
                for orientation in data:
                    if orientation >= min_theta and orientation < max_theta:
                        count+=1
                count_bins.append(count)
                
            # now determine the radii of the bin sectors; unlike a cartesian histogram, the area of a polar histogram should be corrected
            radii_bins = np.sqrt(np.array(count_bins)/N_segments*100)
            rmax = max(radii_bins) #use this later when setting max radius of plot
                
            # convert theta bins and binsize to radians for matplotlib
            binsize_rad = np.pi/number_of_bins
            theta_bins = theta_bins*np.pi/180
                
            # theta should actually be at the center of the bin so add half a binsize to all theta values
            theta_bins += binsize_rad/2
                
            # plot the bins
            if column == 1:
                colour = "limegreen"
            else: 
                colour = "darkorange"
            ax.bar(theta_bins, radii_bins, width=binsize_rad, color=colour, alpha=0.5)
                        
            # Make the polar plot with North up
            ax.set_theta_zero_location("N")
            ax.set_theta_direction("clockwise")
            ax.grid(visible=True)
            ax.set_thetamin(-90)
            ax.set_thetamax(90)

            # set r axis
            ax.set_rmin(0) # set min r to 0
            rmax_sq = rmax**2 # square rmax (see above) to get the data % that this area represents 
            rmax_sq_rounded = round(rmax_sq,-1)
            if rmax_sq_rounded <= rmax_sq:
                rmax_sq_rounded += 10
            rmax = rmax_sq_rounded**0.5 # take the square root of rmax_sq to get the new rmax 
            ax.set_rmax(rmax)
            radii_ticks = np.sqrt(np.arange(0,100,10))[np.sqrt(np.arange(0,100,10)) <= rmax]
            radii_tickmarks = len(radii_ticks)*[""]
            ax.set_rticks(radii_ticks, radii_tickmarks)

            if column == 1 or column == 2:
                ax.set_xticks([-np.pi/2, -np.pi/4, 0, np.pi/4, np.pi/2], labels = ["W", "NW", "N", "NE", "E"])
            if column == 3:
                ax.set_xticks([-np.pi/2, -np.pi/3, -np.pi/6, 0, np.pi/6, np.pi/3, np.pi/2], labels = ["-90", "-60", "-30", "0", "30", "60", "90"])

            # plot the bedding data
            if column == 2:
                # VRP segment orientations, plot mean
                ax.plot([0, mean_strike*np.pi/180], [0, 0.9*rmax], color = "cornflowerblue", linewidth = 2)
                ax.plot([0, mean_dip*np.pi/180], [0, 0.9*rmax], color = "cornflowerblue", linewidth = 2, linestyle = "--")
                # also plot S and D
                ax.text(mean_strike*np.pi/180, rmax+1.0, s = "S", color = "cornflowerblue", ha='center', va='center')
                ax.text(mean_dip*np.pi/180, rmax+1.0, s = "D", color = "cornflowerblue", ha='center', va='center')
            # plot the mean of strike and dip mode of data in column 3
            if column == 3:
                if stream_order in include_means_d_strike_stream_orders:
                    ax.plot([0, d_so_ms_mean_perorder[stream_order]*np.pi/180], [0, 0.9*rmax], color = "#ff6200", linewidth = 2, zorder=3)
                if stream_order in include_means_d_dip_stream_orders:
                    ax.plot([0, d_so_md_mean_perorder[stream_order]*np.pi/180], [0, 0.9*rmax], color = "#ff6200", linewidth = 2, zorder=3)
                    
    # set up text - stream orders and N measurements AP
    gs_text = gridspec.GridSpecFromSubplotSpec(max_stream_order+1,1, subplot_spec=gs_lower[0], hspace = hspace) # for labels of stream order
    for stream_order in stream_orders:
        ax = plt.subplot(gs_text[stream_order])
        # stream orders
        if stream_order == 0:
            ax.text(0,0.25, s = "All segments", rotation = "vertical")
        else:
            ax.text(0,0.25, s = "Order " + str(stream_order), rotation = "vertical")
        # N masurments AP
        ax.text(0.7, 0.25, s = "N = " + str(AP_N_list[stream_order]), rotation = "vertical", color = "limegreen")
        ax.set_axis_off()
        
    # set up text - N measurements VRP
    gs_text = gridspec.GridSpecFromSubplotSpec(max_stream_order+1,1, subplot_spec=gs_lower[4], hspace = hspace) # for labels of stream order
    for stream_order in stream_orders:
        ax = plt.subplot(gs_text[stream_order])
        ax.text(0, 0.25, s = "N = " + str(VRP_N_list[stream_order]), rotation = "vertical", color = "darkorange")
        ax.set_axis_off()

    plt.savefig(destination_directory+"/halfdials.jpg", format = "jpg", dpi = 360)

    
####################################
# DISTRIBUTION GENERATION FUNCTION #
####################################

def generate_distribution(p1, peak, distribution, s = 30, m = 90, n = 1000):
    """
    Function to generate a distribution.
    
    Args:
        p1 (float, [0,1]): the chance of picking a value of a location given by location_of_peak
        peak (int or float): the value of the unimodal peak (given on [0,180)) or the value of the first peak of a bimodal distribution (given on [0,90))
        distribution (str): either "unimodal" or "bimodal". If "unimodal", only 1 peak will be considered. If "bimodal", a second peak will be generated +90 degrees of the value given for the first peak.
        s (int or float): the standard deviation of the distribution from which a value is drawn to modify the originally generated value
        m (int or float): the mean of the distribution from which a value is drawn to modify the originally generated value
        n (int): the number of values in the distribution
        
    Returns: 1xn list 
    """
    
    # set up list for distribution
    distribution_list = []
    
    # determine peak(s) 
    if distribution == "unimodal":
        peaks = [peak]
    elif distribution == "bimodal" or distribution == "skewed_bimodal":
        peaks = [peak, peak+90]
    else:
        print("You didn't specifiy the distribution correctly!")
    
    for _ in range(n):
        # Randomly select one of two or one of one values from 'peaks' with a probability p1, otherwise get a random uniform number on the interval [0,180)
        value = random.choice(peaks) if random.random() < p1 else random.uniform(0, 180)
        if distribution == "skewed_bimodal":
            value = peaks[0] if random.random() <p1 else peaks[1]
        # Modify the value using a normal distribution with mean m and standard deviation s
        modification = np.random.normal(m, s)
        value += modification

        # Adjust the value to be within the range [0, 18)
        if value < 0:
            value += 180
        elif value >= 180:
            value -= 180
        
        # Put value in list
        distribution_list.append(value)
        
    return distribution_list
    
def calculate_kuipers_V(distribution1, distribution2):
    """
    Function to calculate Kuiper's D statistic from two distributions.
    """

    # Create a common set of values for the x-axis
    x_values = np.linspace(0, 180, 1000)
    
    # Calculate the 'empirical' cumulative distribution functions (ECDFs)
    ecdf1 = np.searchsorted(np.sort(distribution1), x_values, side='right') / len(distribution1)
    ecdf2 = np.searchsorted(np.sort(distribution2), x_values, side='right') / len(distribution2)

    # Calculate the differences between ECDFs
    differences = ecdf1 - ecdf2

    # Calculate the maximum positive and maximum negative differences
    max_positive_difference = np.max(np.clip(differences, 0, None))
    max_negative_difference = np.abs(np.min(np.clip(differences, None, 0)))

    # Calculate the combined value V
    V = max_positive_difference + max_negative_difference

    return V

    
#####################################################
# INTERNAL FUNCTIONS ALSO USED BY Dataprocessing.py #
#####################################################

def remove_directionality(directional_measurements):
    """
    Function to remove directionality from a set of directional measurements (on [0,360),in a list)
    
    Returns a list of orientations, on [0,180)
    """
    
    directional_measurements = np.array(directional_measurements)
    
    # We want to remove directionality to get just orientations (see "Statistics & data analysis in Geology", describing method developed by Krumbein (1939))
    #   double all angles to remove directionality
    directional_measurements = directional_measurements*2 
    
    #   all values should now be in the same direction but we need to take out multiples of 360 degrees (i.e. 700 degrees needs to be 700%360=340)
    directional_measurements = directional_measurements % 360
    
    #   values should now be halved again. this process automatically brings all the values to the interval [0,180)
    directional_measurements = directional_measurements/2
    
    #   save as orientations and make list
    orientations = directional_measurements.tolist()
    
    # return orientations
    return orientations

def get_strikes_and_dip_orientations_from_dip_directions(dip_directions):
    """
    Converts dip directions (degrees) into strikes (degrees) on the interval [0,180) 
    
    Arguments:
        dip_directions (list): dip directions in degrees
        
    Returns:
        strikes (list): strikes in degrees on [0,180)
        dip_orientations (list): dip orientations in degrees on [0,180)
    """
    
    # remove directionality from dip directions
    dip_orientations = remove_directionality(dip_directions)
    
    # calculate strikes 
    dip_orientations_arr = np.array(dip_orientations)
    strikes_arr = dip_orientations_arr + 90
    strikes = strikes_arr.tolist()
    strikes = remove_directionality(strikes) # this brings the measurements back to [0,180)
    
    # return the strikes and dip orientations
    return strikes, dip_orientations

def orientations_from_NS_to_WE(orientations_list):
    """
    Puts a list of orientations given on [0,180) (North-South) into orientations on [-90,90) (West-East)
    """
    new_orientations = [orientation - 180 if orientation >= 90 else orientation for orientation in orientations_list]
    return new_orientations
    

def find_diff_segm_or_and_mean_strike_or_dip_or(segment_orientation, mean_strike_or_dip_orientation):
    """
    Calculates the difference between mean strike and segment bearing in degrees, on the interval [-90, 90]
    
    Arguments:
        segment_orientation (float): the orientation of the segment on the interval [0,180)
        mean_strike_or_dip_orientation (float): mean_strike or mean_dip_orientation of bedding in degrees on the interval [0,180)
        
    Returns:
        d (float): the difference between the two inputs
    """
    # calculate difference
    d = segment_orientation - mean_strike_or_dip_orientation
    
    # correct d such that result lies on [-90, 90] instead of [-180, 180]
    if d>90:
        d -= 180
    elif d<-90:
        d += 180
    
    return d
        
def fix_list(lst):
    """
    This short function turns a string representation of a list into a list.
    Should be used when importing a csv with columns containing lists (i.e. the _bedding_stats.csv files)
    When importing as dataframe, you can call this function:
    gdf['column_name'] = gdf.apply(lambda row: fix_list(row.column_name), axis = 1)
    To fix the lists.
    """
    
    lst=lst.replace("'","")   #removes any ' from the data in case they are there. 
    
    if lst == "[]": # deals with empty lists
        lst = [] # set list
    else: # deals with non empty lists
        lst = list(map(int, map(float, lst.strip('][').split(',')))) #set list
    return lst