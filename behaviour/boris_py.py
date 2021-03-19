# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 09:44:17 2019
BORIS csv reader
@author: Chris Boccia
"""

"""
    In the GitHub repo, a list of BORIS files is already provided
    The following instructions only apply if you've reanalyzed videos files
    Or changed the BORIS input csvs

    To generate a list of BORIS files..
    NOTE: Ensure that the file folder only contains BORIS csvs that
    you wish to analyze; otherwise you will need to remove those
    files from the resulting text file by hand. The goal here is to 
    create a list of file names with extensions, one per line
    
   In Windows: Open a command prompt window one folder level above
   your BORIS csv output file folder and run the following command 
   (replace <csv_directory> with your folder name)
   
   dir <csv_directory> /B > boris_names.txt
   
   In *NIX/Mac: Open a terminal/Bash shell one folder level above
   your BORIS csv output file folder, and run the following command 
   (replace <csv_directory> with your folder name, and remember 
   to unalias ls if you're using additonal output flags by default )
   
   ls <csv_directory> > boris_names.txt
"""

"""Run this python script from within your working directory (or alter
    the paths in the lines below to suit your filesystem)"""

"""Script consists of two functions that are each run once"""   

"""get list of BORIS files from text file"""
def get_names():
    f_l = list()
    """Change path to suit wd """
    """reads through file, removes end of line characters, assigns each file name to a list entry"""
    
    n_file = open("./boris_csvs/boris_names.txt")
    
    for x in n_file:
        f_l.append(x.rstrip('\n'))
        """print(f_l) uncomment if you need to check files"""
    """close file and return list of cleaned file names"""
    n_file.close()
    return(f_l)
        
"""run above function to get file names/paths"""
f_list = get_names()

"""Function to process a list of BORIS files"""
def csv_test(f_list):
    
    """import packages """
    import csv
    import statistics
    import os
    
    """open output file"""
    fpr = open("./boris_summary.csv", "w")
    
    """write output csv file header"""
    fpr.write("VideoID, ID,Start,End,Duration,num_all_rb,rb_full_rate,rb_full_int,first_rb,last_rb,rb_period,rate_adj_rb,low_rb, low_rb_rate, num_br,rate_br,num_gul,rate_gul,num_ch,rate_ch,coov,toov,hoov,dur_exh_tot,exh_pct,num_exh, exh_avg,soh,toh,hbn,lbn\n")
    """open directory containing all files; adjust to suit your file system; currently only need for output"""
    os.chdir("./boris_csvs")
    tot_intc = 0
    
    for e in f_list:
        """print(e)
        uncomment print statement if you need to identify which csv is causing
        errors to occur"""
        
        """initialize variables; need to have a string for printing to file
        that will summarize the data for the file being processed; to generate
        this data, need to assemble vectors of times, event IDs, event states,
        comments; also need a placeholder for individual ID to sync metadata later"""
        hold_string = "0"
        times = []
        events = []
        state = []
        comments = []
        indiv = "placehold"
        """open file for processing; following lines extract all events, times etc"""
        with open (e) as fil:
            for q in range(16):
                fil.readline()
            csvReader = csv.reader(fil)
            for row in csvReader:
                indiv = row[4]
                times.append(row[0])
                events.append(row[5])
                comments.append(row[7])
                state.append(row[8])
            fil.close()
        
        """add individual to data string"""
        hold_string = e +","+indiv+","
        
        """Following lines calculate the amount of time the head was out of view.
        These events occur in pairs (start, end) so the script takes the 
        cumulative amounts of time mapped out from the 1st member of a pair to
        the second. Possible cases: >1 set of pairs, 1 pair, no pairs"""
        hoov_in = [i for i, x in enumerate(events) if x == "head out of view"]
        if(len(hoov_in) >2):
            s_hoov = 0
            for u in range((len(hoov_in)-1)):
                s_hoov += float(times[hoov_in[u+1]]) - float(times[hoov_in[u]])
        elif(len(hoov_in) > 0):
            s_hoov = float(times[hoov_in[1]]) - float(times[hoov_in[0]])
        else:
            s_hoov = 0
        
        """get time trial started from event vector"""
        s = events.index("start of trial")
        start = float(times[s])
        
        """print start time to data string"""
        hold_string = hold_string+times[s]+","
        
        """get time trial ended from event vector"""
        s = events.index("end of trial")
        end = float(times[s])
        
        """if events[-1]!="end of trial" and times[-1] != times[-2]:
            prob.write("end prob:" + e + "\n")"""
        
        """print end time to data string"""
        hold_string = hold_string+times[s]+","
        
        """calculate duration, print to data string"""
        dur = end-start
        hold_string = hold_string+str(dur)+","
        
        """get indices for each type of reinhale"""
        indices = [i for i, x in enumerate(events) if x == "reinhale"]
        indices_lbn = [i for i, x in enumerate(comments) if x == "lbn"]
        indices_soh = [i for i, x in enumerate(comments) if x == "soh"]
        indices_hbn = [i for i, x in enumerate(comments) if x == "hbn"]
        indices_toh = [i for i, x in enumerate(comments) if x == "toh"]
        
        """calculate numbers of each type"""
        low_rein = len(indices_lbn)
        soh = len(indices_soh)
        hbn = len(indices_hbn)
        toh = len(indices_toh)
        lbn = low_rein
        
        """remove lbn reinhales from data set"""
        if (len(indices_lbn) > 0):
            n_v = len(indices) - len(indices_lbn)
            if(n_v==0):
                indices = ''
            else:
                indices = set(indices) - set(indices_lbn)
                indices = list(indices)
                
        """calculate number of total reinhales (without lbn), print to file"""
        rein = len(indices)
        hold_string = hold_string+str(rein)+","
        
        """if head was out of view for some duration, subtract that from consideration 
        time; calculate and print rebreathing rate"""
        if((dur-s_hoov) != 0):
            rb_r = rein/(dur-s_hoov)
        else:
            rb_r = "NA"
        hold_string = hold_string+str(rb_r)+","
        
        """get all intervals between rebreathing events; get average; print average
        or absence of data to string"""
        intervals = []
        cc=0
        while cc < (len(indices)-1):
            intervals.append(abs(float(times[indices[cc+1]])-float(times[indices[cc]])))
            cc+=1
        if(len(intervals) > 0):
            tot_intc += len(intervals)
            m_i = statistics.mean(intervals)
        else:
            m_i = "NA"
        hold_string = hold_string+str(m_i)+","
        """find first and last reinhales (if >1 exist); calculate
        and print rebreathing period (duration of time between first and last reinhale"""
        if(rein > 1):
            rf_1 = times[indices[0]]
            rf_last = times[indices[(rein-1)]]
            rb_p = abs(float(rf_last)-float(rf_1))
        elif(rein == 1):
            rf_1 = times[indices[0]]
            rf_last = "NA"
            rb_p = 0  
        else:
            rf_1 = "NA"
            rf_last = "NA"
            rb_p = 0    
        hold_string = hold_string+str(rf_1)+","
        hold_string = hold_string+str(rf_last)+","
        
        """calculate and print rebreathing rate during rebreating period only"""
        if(rb_p > 0):
                cor = rein/float(rb_p)
        else:
            cor = "NA"
        if(rb_p == 0):
            rb_p = "NA"
        else:
            rb_p = rb_p
        hold_string = hold_string+str(rb_p)+","
        hold_string = hold_string+str(cor)+","
        
        """calculate and print lbn number, rate; remove hov time from rate"""
        if((dur-s_hoov) != 0):
            rel_low_rein = low_rein / (dur-s_hoov)
        else:
            rel_low_rein = "NA"
        hold_string = hold_string+str(low_rein)+","
        hold_string = hold_string+str(rel_low_rein)+","
        
        """get indices of all bubble releases; calculate and print
        number and rate, removing hov time"""
        bub_r_indices = [i for i, x in enumerate(events) if x == "bubble release"]
        tot_bub_r = len(bub_r_indices)
        if((dur-s_hoov) != 0):
            bub_r_rel = tot_bub_r / (dur-s_hoov)
        else:
            bub_r_rel = "NA"
        hold_string = hold_string+str(tot_bub_r)+","
        hold_string = hold_string+str(bub_r_rel)+","
    
        """get time throat was out of view; see hov above for explanation"""
        toov_in = [i for i, x in enumerate(events) if x == "throat out of view"]
        if(len(toov_in) >2):
            s_toov = 0
            for u in range((len(toov_in)-1)):
                s_toov += float(times[toov_in[u+1]]) - float(times[toov_in[u]])
        elif(len(toov_in) > 0):
            s_toov = float(times[toov_in[1]]) - float(times[toov_in[0]])
        else:
            s_toov = 0
        """get gular pumping instances, calculate and print number and rate, 
        removing tov time"""
        gul_indices = [i for i, x in enumerate(events) if x == "gular pump"]
        tot_gul = len(gul_indices)
        if((dur-s_toov)!=0):
            gul_rel = tot_gul/ (dur-s_toov)
        else:
            gul_rel = "NA"
        hold_string = hold_string+str(tot_gul)+","
        hold_string = hold_string+str(gul_rel)+","
    
        """get chest out of view time; see hov for explanation"""
        coov_in = [i for i, x in enumerate(events) if x == "chest out of view"]
        if(len(coov_in) >2):
            s_coov = 0
            for u in range((len(coov_in)-1)):
                s_coov += float(times[coov_in[u+1]]) - float(times[coov_in[u]])
        elif(len(coov_in)==2):
            s_coov = float(times[coov_in[1]]) - float(times[coov_in[0]])
        else:
            s_coov = 0
        
        """get all chest movement events, calculate and print number
        and rate, removing cov time; also print all out of view times"""
        ch_indices = [i for i, x in enumerate(events) if x == "chest"]
        tot_ch = len(ch_indices)
        if((dur-s_coov) != 0):
            ch_rel = tot_ch/ (dur-s_coov)
        else:
            ch_rel = "NA"
        hold_string = hold_string+str(tot_ch)+","
        hold_string = hold_string+str(ch_rel)+","
        hold_string = hold_string+str(s_coov)+","        
        hold_string = hold_string+str(s_toov)+","        
        hold_string = hold_string+str(s_hoov)+","
        
        """get exhale indices, calculate sum of time a bubble was left exhaled
        by taking the time between each exhale start/end state event. Calculate
        and print sum of exhalation time, as well as percent of trial with
        a bubble exhaled (removing hov time); also calculate number of exhales
        and average duration of each exhale"""
        ex_indices = [i for i, x in enumerate(events) if x == "exhale"]
        ex_intervals = []
        sum_ex = 0
        cc=0
        while cc < (len(ex_indices)-1):
            val = (float(times[ex_indices[cc+1]])-float(times[ex_indices[cc]]))
            sum_ex += val
            ex_intervals.append(val)
            cc+=2
        if(len(ex_intervals) > 0):
            m_ex = statistics.mean(ex_intervals)
        else:
            m_ex = "NA"
        if((dur-s_hoov) > 0):
            ex_pct = sum_ex / (dur-s_hoov)
        else:
            ex_pct = "NA"
        hold_string = hold_string+str(sum_ex)+","
        hold_string = hold_string+str(ex_pct)+","
        hold_string = hold_string+str(len(ex_indices)/2)+","                
        hold_string = hold_string+str(m_ex)+","
        """print numbers of each type of rebreathing bubble to file"""
        hold_string = hold_string+str(soh)+"," + str(toh)+"," + str(hbn)+"," + str (lbn)+"\n"
            
        """write data to file, then move on to next file"""    
        fpr.write(hold_string)
    
    """close output file"""
    fpr.close()
    
    """the following lines of code generate rebreathing intervals for each individual
    for distribution analysis"""
    """Interval data were not included in the paper"""
    f2 = open("./intervals.csv", "w")
    fprint = "ID\n"
    f2.write(fprint)
    
    for e in f_list:
        """initialize variables"""
        hold_string = "0"
        times = []
        events = []
        state = []
        comments = []
        indiv = "placehold"
        """open file for processing"""
        with open (e) as fil:
            for q in range(16):
                fil.readline()
            csvReader = csv.reader(fil)
            for row in csvReader:
                indiv = row[4]
                times.append(row[0])
                events.append(row[5])
                comments.append(row[7])
                state.append(row[8])
        fil.close()
        """get reinhale events only"""
        indices = [i for i, x in enumerate(events) if x == "reinhale"]
        """find lbn reinhales"""
        indices_lbn = [i for i, x in enumerate(comments) if x == "lbn"]
        """get number of lbns"""
        low_rein = len(indices_lbn)
        """remove lbns from consideration from interval calculations"""
        if (len(indices_lbn) > 0):
            """if all rebreathing events are lbns, don't calculate intervals"""
            n_v = len(indices) - len(indices_lbn)
            if(n_v==0):
                indices = ''
            else:
                indices = set(indices) - set(indices_lbn)
                indices = list(indices)
        intervals = []
        cc=0
        """add individual ID, all intervals to a string for output"""
        fprint = str(indiv)
        while cc < (len(indices)-1):
            intervals.append(abs(float(times[indices[cc+1]])-float(times[indices[cc]])))
            fprint = fprint + "," + str(abs(float(times[indices[cc+1]])-float(times[indices[cc]])))
            cc+=1
        fprint = fprint + "\n"
        f2.write(fprint)
    """close interval file"""    
    f2.close()
    
"""run function"""
csv_test(f_list)



