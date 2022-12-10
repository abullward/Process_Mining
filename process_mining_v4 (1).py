#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pm4py as pm4py
from datetime import datetime
import seaborn as sns
import matplotlib.pyplot as plt
from statistics import mean
from statistics import mode
from statistics import median
from pm4py.algo.filtering.log.variants import variants_filter
from pm4py.statistics.traces.log import case_statistics
import numpy as np

def load_events(Cancer_code_icd10):
    ###build an event file from cancer sim
    ###load files
    av_tumour = pd.read_csv('C:/Users/Ali/Documents/Mres R/simulacrum_release_v1.2.0.2017/simulacrum_release_v1.2.0.2017/data/sim_av_tumour.csv') 
    sact_cycle = pd.read_csv('C:/Users/Ali/Documents/Mres R/simulacrum_release_v1.2.0.2017/simulacrum_release_v1.2.0.2017/data/sim_sact_cycle.csv') 
    sact_cycle= sact_cycle.rename(columns={"MERGED_PATIENT_ID":"PATIENTID"})
    deaths = pd.read_csv('C:/Users/Ali/Documents/Mres R/simulacrum_release_v1.2.0.2017/simulacrum_release_v1.2.0.2017/data/sim_av_patient.csv') 
    sact_regimen = pd.read_csv('C:/Users/Ali/Documents/Mres R/simulacrum_release_v1.2.0.2017/simulacrum_release_v1.2.0.2017/data/sim_sact_regimen.csv')
    sact_regimen= sact_regimen.rename(columns={"MERGED_PATIENT_ID":"PATIENTID"})

    #################################filter by icd10
    ###filter table to ICD10 code
    icd10 = av_tumour['SITE_ICD10_O2_3CHAR'] == Cancer_code_icd10 #brain
    #now filter table
    filter_tab = av_tumour[icd10] ### need uncomment to start filtering again
    #display(filter_tab)

    # creat a patient index for given icd10
    patients_df = filter_tab[['PATIENTID']] ### 
    #display(patients_df)

    ####start building events
    #for diagnosis events
    diagnosis_events = patients_df.merge(filter_tab,on='PATIENTID',how='left')            
    #display(diagnosis_events)
    #only select relevent columns
    diagnosis_events = diagnosis_events[['PATIENTID','DIAGNOSISDATEBEST']]
    #add event type column
    diagnosis_events['Activity']='Diagnosis'
    ##rename columns
    diagnosis_events = diagnosis_events.rename(columns={"PATIENTID":"case_id","DIAGNOSISDATEBEST":"timestamp"})

    #now do for first surgery

    #left join to table
    firstsurgery_events = patients_df.merge(filter_tab,on='PATIENTID',how='left')            
    #display(diagnosis_events)
    #only select relevent columns
    firstsurgery_events = firstsurgery_events[['PATIENTID','DATE_FIRST_SURGERY']]
    #add event type column
    firstsurgery_events['Activity']='First_Surgery'
    ##rename columns
    firstsurgery_events = firstsurgery_events.rename(columns={"PATIENTID":"case_id","DATE_FIRST_SURGERY":"timestamp"})

    #display(firstsurgery_events)

    ###now add sact cycles

    #left join to table
    sact_events = patients_df.merge(sact_cycle,on='PATIENTID',how='left')            
    #display(sact_events)
    #only select relevent columns
    sact_events = sact_events[['PATIENTID','START_DATE_OF_CYCLE']]
    #add event type column
    sact_events['Activity']='Sact_Cycle'
    ##rename columns
    sact_events = sact_events.rename(columns={"PATIENTID":"case_id","START_DATE_OF_CYCLE":"timestamp"})

    #display(sact_events)

    ###now add deaths
    #display(deaths)

    #left join to table
    death_events = patients_df.merge(deaths,on='PATIENTID',how='left')            
    #remove alive patients
    death_only = death_events['NEWVITALSTATUS']=='D'
    death_events = death_events[death_only]
    #display(death_events)


    #only select relevent columns
    death_events = death_events[['PATIENTID','VITALSTATUSDATE']]
    #add event type column
    death_events['Activity']='Death'
    ##rename columns
    death_events = death_events.rename(columns={"PATIENTID":"case_id","VITALSTATUSDATE":"timestamp"})
    #filter out live

    #display(death_events)

    ###now add alive patients
    #display(deaths)

    #left join to table
    alive_events = patients_df.merge(deaths,on='PATIENTID',how='left')            
    #remove alive patients
    alive_only = alive_events['NEWVITALSTATUS']=='A'
    alive_events = alive_events[death_only]
    #display(death_events)


    #only select relevent columns
    alive_events = alive_events[['PATIENTID','VITALSTATUSDATE']]
    #add event type column
    alive_events['Activity']='Alive'
    ##rename columns
    alive_events = alive_events.rename(columns={"PATIENTID":"case_id","VITALSTATUSDATE":"timestamp"})
    #filter out live

    #display(alive_events)

    ###now add sact regimens decision to treat

    #left join to table
    decision_to_treat_events = patients_df.merge(sact_regimen,on='PATIENTID',how='left')            

    #display(diagnosis_events)
    #only select relevent columns
    decision_to_treat_events = decision_to_treat_events[['PATIENTID','DATE_DECISION_TO_TREAT']]
    #add event type column
    decision_to_treat_events['Activity']='Decision_to_treat'
    decision_to_treat_events= decision_to_treat_events.rename(columns={"PATIENTID":"case_id","DATE_DECISION_TO_TREAT":"timestamp"})
    #display(decision_to_treat_events)

    ###now add sact 	START_DATE_OF_REGIMEN

    #left join to table
    START_DATE_OF_REGIMEN_events = patients_df.merge(sact_regimen,on='PATIENTID',how='left')            

    #display(diagnosis_events)
    #only select relevent columns
    START_DATE_OF_REGIMEN_events = START_DATE_OF_REGIMEN_events[['PATIENTID','START_DATE_OF_REGIMEN']]
    #add event type column
    START_DATE_OF_REGIMEN_events['Activity']='start_date_of_regimen'
    START_DATE_OF_REGIMEN_events= START_DATE_OF_REGIMEN_events.rename(columns={"PATIENTID":"case_id","START_DATE_OF_REGIMEN":"timestamp"})
    #display(START_DATE_OF_REGIMEN_events)

    #stack tables
    #cancer_events = pd.concat([diagnosis_events,firstsurgery_events,sact_events,decision_to_treat_events,START_DATE_OF_REGIMEN_events])
    cancer_events = pd.concat([diagnosis_events,firstsurgery_events,sact_events,death_events,decision_to_treat_events,START_DATE_OF_REGIMEN_events])
    #remove nas
    cancer_events = cancer_events[cancer_events['timestamp'].notna()]
    #display(cancer_events)
    cancer_events.to_csv('C:/Users/Ali/Documents/process mining/events.csv')

#### A function for plotting process maps
def process_map():
    ##process map
    event_log = pd.read_csv('C:/Users/Ali/Documents/process mining/events.csv', sep=',')
    num_events = len(event_log)
    num_cases = len(event_log.case_id.unique())
    print("Number of events: {}\nNumber of cases: {}".format(num_events, num_cases))

    df = pm4py.format_dataframe(event_log, case_id='case_id', activity_key='Activity', timestamp_key='timestamp')

    dfg, start_activities, end_activities = pm4py.discover_dfg(df)
    pm4py.view_dfg(dfg, start_activities, end_activities)
    
    process_tree = pm4py.discover_tree_inductive(event_log)
    #bpmn_model = pm4py.convert_to_bpmn(process_tree)
    #pm4py.view_bpmn(bpmn_model)
    
    ##### get most common traces by length
    variants = variants_filter.get_variants(event_log)
  
    print(f"There are :{len(variants)} variants in the log")

    variants_count = case_statistics.get_variant_statistics(event_log)
    variants_count = sorted(variants_count, key=lambda x: x['count'], reverse=True)
    ## Printing the top 10 variants by case number
    variants = variants_count#[:500] 
    ##split out counts for a hist and other stats
    hist_plot = []
    event_list = []
    output_for_csv = []
    for z in variants:
       # print(z)
        split = str(z).split(":")[2][:-1].strip()
        event_list = str(z).split(":")[1].split("'")[1].strip()
        #print(event_list)
        #{'variant': 'Diagnosis,First_Surgery,Death', 'count': 6545}
        hist_plot.append(split)
        output_for_csv.append(event_list+":" + split)
        #print(output_for_csv)
    #plt.hist(hist_plot, bins=500, alpha=0.5)
   # sns.distplot(a=hist_plot,bins=5).set_title("Distribution of Traces")
    #plt.show
    #pd.DataFrame(output_for_csv).to_csv('C:/Users/Ali/Documents/process mining/traces.csv')

  #  plt.show()
    ###
    print("#### event to event stats")
    print(pm4py.discover_dfg(event_log))
    filtered = pm4py.filter_start_activities(event_log, {'Decision_to_treat'})
    print(filtered)
       
### a function for density plots
def expl_plots(dates_from,dates_to,title_text):
    ###explore density distributions between dates
    density_fact = dates_from.merge(dates_to,on='case_id',how='left') 
    #strip out nans
    density_fact = density_fact[density_fact['timestamp_y'].notna()]
    density_fact = density_fact[density_fact['timestamp_x'].notna()]
    list_of_days = []
    ##add calculation for days
    #print(density_fact)
    for i, row in density_fact.iterrows():
   
        list_of_days.append(datetime.strptime(row['timestamp_y'],'%Y-%m-%d')- datetime.strptime(row['timestamp_x'],'%Y-%m-%d'))
   

    for_plot = []
    count = 0
    string_val = []
    for a in list_of_days:
    #print(str.lstrip(str(list_of_days[count]))[0:3])
        string_val = str.lstrip(str(list_of_days[count]))[0:2]
        if string_val != "0:":
            string_val =for_plot.append(string_val)
      #  print(string_val)  
        count= count+1
   
    ###plt density chart
    plt.figure()
    sns.distplot(a=for_plot).set_title(title_text)
    plt.show
    print("lenght of",title_text,len(for_plot))
    res = [int(item) for item in for_plot]
    print("Average for", title_text, sum(res)/len(res))
 
def transition_table(title):
    ########### A function for generating the transition tables from the event log

    ###load event log
    event_log = pd.read_csv('C:/Users/Ali/Documents/process mining/events.csv', sep=',')

    df = pm4py.format_dataframe(event_log, case_id='case_id', activity_key='Activity', timestamp_key='timestamp')

    dfg, start_activities, end_activities = pm4py.discover_dfg(df)

    events_to_events = pm4py.discover_dfg(event_log)

    ### create empty table

    heat_map = pd.DataFrame(columns=['Start','Diagnosis','Decision_to_treat','First_Surgery','start_date_of_regimen','Sact_Cycle','Death'],index=range(8))


    ###### clean events for plotting
    list_of_stuff = str(events_to_events[0])

    list_of_stuff = list_of_stuff.replace('(', '')
    list_of_stuff = list_of_stuff.replace(')', '')
    list_of_stuff = list_of_stuff.replace('{', '')
    list_of_stuff = list_of_stuff.replace('}', '')
    list_of_stuff = list_of_stuff.replace(':', ',')
    list_of_stuff = list_of_stuff.replace('"', '')
    list_of_stuff = list_of_stuff.replace("'", '')
    list_of_stuff = list_of_stuff.replace(" ", '')
    #print(list_of_stuff)
    ## replace with refs
    event_types = ['Start','Diagnosis','Decision_to_treat','First_Surgery','start_date_of_regimen','Sact_Cycle','Death']
    n=0
    #print(list_of_stuff)
    for a in event_types:
        #print(a)
        list_of_stuff = list_of_stuff.replace(a, str(n))
        n=n+1

    list_of_stuff = list_of_stuff.split(",")

    #print(list_of_stuff)

    ##now write to table
    x = 0
    y=1
    val =2
    itterations = len(list_of_stuff)/3
    for z in range(int(itterations)):
        heat_map.iat[int(list_of_stuff[x]),int(list_of_stuff[y]) ]= int(list_of_stuff[val])
        x = x+3
        y= y+3
        val = val+3

    heat_map = heat_map.fillna(0)
    #display(heat_map)
    heat_map = heat_map.drop(labels =7, axis=0)
    heat_map.index= ['Start','Diagnosis','Decision_to_treat','First_Surgery','start_date_of_regimen','Sact_Cycle','Death']
    #print(heat_map)
    
    ##populate start stats
    from pm4py.algo.filtering.log.start_activities import start_activities_filter
    log_start = start_activities_filter.get_start_activities(df)
    print(log_start)
    #{'Diagnosis': 200077, 'Sact_Cycle': 2294, 'Decision_to_treat': 16409, 'First_Surgery': 60, 'start_date_of_regimen': 379}

    heat_map.at['Start','Diagnosis'] = int(str(log_start).split(',')[0].split(':')[1])
    heat_map.at['Start','Decision_to_treat'] = int(str(log_start).split(',')[2].split(':')[1])
    heat_map.at['Start','First_Surgery'] = int(str(log_start).split(',')[3].split(':')[1])
    
    heat_map.at['Start','start_date_of_regimen'] = int(str(log_start).split(',')[4].split(':')[1].strip()[0:-1])
    heat_map.at['Start','Sact_Cycle'] = int(str(log_start).split(',')[1].split(':')[1])
    heat_map.at['Start','Death'] = 0#int(str(log_start).split(',')[1].split(':')[1])
  
   
    ### construct heat map
    s = sns.heatmap(heat_map, annot=True,fmt="d", cmap="YlGnBu")
    s.set_xlabel('To', fontsize=10)
    s.set_ylabel('From', fontsize=10)
    s.set_title(title)
    display(s)
    print("Total Transitions:", heat_map.sum())

def traces(label):
    ###some code that actually generates traces
     ##### get most common traces by length

    ##load event log
    event_log = pd.read_csv('C:/Users/Ali/Documents/process mining/events.csv', sep=',')
    df = pm4py.format_dataframe(event_log, case_id='case_id', activity_key='Activity', timestamp_key='timestamp')
    dfg, start_activities, end_activities = pm4py.discover_dfg(df)


    ###generate traces
    variants = variants_filter.get_variants(df)

    print(f"There are have:{len(variants)} unique traces")

    variants_count = case_statistics.get_variant_statistics(df)
    variants_count = sorted(variants_count, key=lambda x: x['count'], reverse=True)
    ## Printing the top 10 variants by case number
    print(f"Top 10 Traces are:")
    top_ten_traces = variants_count[:10] 
    print(top_ten_traces)

    hist_plot = []
    event_list = []
    output_for_csv = []
    for z in variants_count:
        print(z)
        #{'variant': 'Diagnosis,First_Surgery,Death', 'count': 6545}
        split = str(z).split(":")[2].strip("}")
        print(split)
        event_list = str(z).split(":")[1].split("'")[1].strip()
        print(event_list)
        #{'variant': 'Diagnosis,First_Surgery,Death', 'count': 6545}
        hist_plot.append(split)
        output_for_csv.append(event_list+":" + split)
        #print(output_for_csv)
        #plt.hist(hist_plot, bins=500, alpha=0.5)
       # sns.distplot(a=hist_plot,bins=5).set_title("Distribution of Traces")
        #plt.show
    pd.DataFrame(output_for_csv).to_csv('C:/Users/Ali/Documents/process mining/traces_'+label+'.csv')
    print(f"All Traces have been printed to csv")

####Check Heights script
# This is a script that loads the "sim_sact_regimen" table and identifies number of pathways where height has shrun over a pathway over a certain threshold (set to 5cm) 
# Last modified 13/12/2022


##load sact table
sact_regimen = pd.read_csv('C:/Users/Ali/Documents/Mres R/simulacrum_release_v1.2.0.2017/simulacrum_release_v1.2.0.2017/data/sim_sact_regimen.csv')
sact_regimen= sact_regimen.rename(columns={"MERGED_PATIENT_ID":"PATIENTID"})
ids = sact_regimen['PATIENTID']
#print(ids)
###extract patients with over xx events

from collections import Counter
counts = Counter(ids)
dupids = [id for id in ids if counts[id] > 2] ##  This defines the min lenghth of pathway to be considered.  
dupids = list(set(dupids)) 

#print(dupids)
#print(len(dupids))
shrink = 0
shrink_count = 0
# convert start date of regiment to date time

sact_regimen["START_DATE_OF_REGIMEN"] = pd.to_datetime(sact_regimen["START_DATE_OF_REGIMEN"], errors = 'coerce') # coerce returns dates before 1600 as NaT
sact_regimen.sort_values(by='START_DATE_OF_REGIMEN', inplace=True)

#display(filter_tab)
for a in dupids:
    search = sact_regimen['PATIENTID'] == a
    filter_tab = sact_regimen[search]
    #display(filter_tab)
   # filter_tab["START_DATE_OF_REGIMEN"] = pd.to_datetime(filter_tab["START_DATE_OF_REGIMEN"])
    #filter_tab.sort_values(by='START_DATE_OF_REGIMEN', inplace=True)
    
    # filter to one patient
    list_of_patients = filter_tab['PATIENTID']

    #print(len(list_of_patients))
    
    heights = filter_tab['HEIGHT_AT_START_OF_REGIMEN']
    heights  = heights[~np.isnan(heights)].tolist() #turn nans to 0
    heights = [i for i in heights if i != 0] ##remove 0s
    
    #list(heights)
    #print(heights)
    count = 0
    shrink = 0
    for z in heights:   
        for b in heights[count:len(heights)]:
           # print("a",a)
           # print("b",b)
            if z-b>0.05:   #####adjust here for shrinkage threshold in meters
                shrink = shrink+1
                
                break          ##break loop on finding one shrink
                #print("#################################shrink")
        count = count+1
        
       # print("#################")
    if shrink > 0:
        shrink_count = shrink_count+1
    #print("shrink count",shrink_count)
    
   # print(heights, "shrink",shrink_count)
    
print("out of ", len(dupids), "pathways there are ", shrink_count, "pathways with a shrink. This is", (shrink_count/len(dupids))*100,"%" )

## function for identifying deaths out out of place in a pathway
def pathway_stats(label):
    pathways = pd.read_csv('C:/Users/Ali/Documents/process mining/traces_'+label+'.csv', sep=',')
    a =pathways['0']
    pathway_length = []
    deaths_outofplace = []
    for z in a:
        pathway_length.append(z.count(',')+1)
        deaths_outofplace.append(z.split(':')[0][:-5].count('Death'))
       # print(z.split(':')[0][:-5])
    
    #print(label+ ' Number of Events',sum(pathway_length))  
    #print(label+ ' Number of Unique Pathways',len(pathway_length))
    #print(label+ ' Max pathway Length',max(pathway_length))
    #print(label+ ' Pathway Min',min(pathway_length))
    #print(label+ ' Pathway Length (mean)',mean(pathway_length))
    #print(label+ ' Pathway Length (mode)',mode(pathway_length))
    #print(label+ ' Pathway Length (median)',median(pathway_length))
    #print(label+ ' Standard Deviation',np.std(pathway_length))
    
    print(label+ ' Death out of place',sum(deaths_outofplace),'out of',len(deaths_outofplace),'pathways which is ', 100*(sum(deaths_outofplace)/len(deaths_outofplace)),'%')
    
    
    sns.distplot(a=pathway_length,bins=20).set_title("Distribution of Pathway Lengths for "+label)
    #display(pathways)


########################## Below is some sample analysis using the defined functions above ###
## the below generates pathway stats for deaths out of place for C50. Can be changed for other top 5 tumour types C44,D09,C34,C61
load_events('D06')
transition_table('D06')



######################################top 5 cancer analyis


# In[19]:


#load_events('C50')
#traces('C50')

load_events('D06')
traces('D06')

#load_events('C34')
#traces('C34')

#load_events('C61')
#traces('C61')


# In[8]:


###### events and process
load_events('C44')
process_map()
transition_table('MALIGNANT NEOPLASM OF SKIN (C44)')
#traces()

### A block on diagnosis
#expl_plots(diagnosis_events,decision_to_treat_events,"Days Between Diagnosis and Treat Events") 
#expl_plots(diagnosis_events,firstsurgery_events,"Days Between Diagnosis and First Surgery")
#expl_plots(diagnosis_events,death_events,"Days Between Diagnosis and Death")
#expl_plots(diagnosis_events,START_DATE_OF_REGIMEN_events,"Days Between Diagnosis and Start of Regimen")


## other useful plots
#expl_plots(firstsurgery_events,death_events,"Days Between First Surgery and Death")
#expl_plots(decision_to_treat_events,death_events,"Days Between Decision to Treat and Death")

#[diagnosis_events,firstsurgery_events,sact_events,death_events,decision_to_treat_events,START_DATE_OF_REGIMEN_events


# In[9]:


load_events('C50')
process_map()
transition_table('MALIGNANT NEOPLASM OF BREAST (C50)')


# In[30]:


load_events('C61')
process_map()
transition_table('MALIGNANT NEOPLASM OF PROSTATE(C61) ')


# In[31]:


load_events('C34')
process_map()
transition_table('MALIGNANT NEOPLASM OF LUNG(C34)')


# In[36]:


load_events('D06')
process_map()
transition_table('CARCINOMA-IN-SITU OF CERVIX UTERI(D06)')


# In[46]:


###get some summary stats on pathways


# In[47]:


pathway_stats('C44')


# In[48]:


pathway_stats('C44')
pathway_stats('C50')
pathway_stats('C61')
pathway_stats('C34')
pathway_stats('D06')


# In[ ]:


counts,values = pd.Series(df).value_counts().values, pd.Series(df).value_counts().index
df_results = pd.DataFrame(list(zip(values,counts)),columns=["value","count"])
display(df_results)
df_results.to_csv('C:/Users/Ali/Documents/process mining/num_cancer_tumours.csv')


# In[3]:





# In[30]:


###load event log
event_log = pd.read_csv('C:/Users/Ali/Documents/process mining/events.csv', sep=',')

df = pm4py.format_dataframe(event_log, case_id='case_id', activity_key='Activity', timestamp_key='timestamp')

#dfg, start_activities, end_activities = pm4py.discover_dfg(df)

#events_to_events = pm4py.discover_dfg(event_log)

from pm4py.algo.filtering.log.start_activities import start_activities_filter
log_start = start_activities_filter.get_start_activities(df)


# In[13]:


print(log_start)
print(str(log_start).split(',')[0].split(':')[1])


# In[ ]:


process_map()


# In[7]:


# a bit of discovery on tumour table
 ###plt density chart
plt.figure()
sns.distplot(filter_tab[['AGE']],bins=19).set_title("Age at Diagnosis for Breast Cancer (C50)")
plt.show

##put into age buckets
ages = filter_tab[['AGE']]
bins = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95]
labels = ["0","5","10","15","20","25","30","35","40","45","50","55","60","65","70","75","80","85","90"]
ages['agerange']= pd.cut(ages.AGE,bins, labels=labels, include_lowest=True)

grouped = ages.groupby(by="agerange").sum()
grouped.to_csv('C:/Users/Ali/Documents/process mining/grouped_c50.csv')
print(grouped)
#filter_tab.columns

