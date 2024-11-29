#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: C. Martinez-Sebastian

GOAL : Creating the summary table 
"""

import csv #import to use the csv module
from csv import DictWriter
import numpy as np

import pandas as pd 
path = '/media/cms/Elements/PhD/DATA'
def add_star(star_id):
    ###################### Charging the list to correlate #####################
    RVpp = pd.read_excel(path+'/FASTWIND/star_analysis/FINAL_TABLE_TO_PLOTS_mar2023_ALL_parameters.xlsx',usecols=[1,8,10,12,15])
    Gonzalos_parameters = pd.read_csv(path+'/SPECTRA/BDD_IACOB_030222-BDD_IACOB_100121.csv')
    best_models_txt = pd.read_csv(path+'/SPECTRA/Tabla_fwlab_OStars_V1.txt',delimiter=r"\s+",skiprows=[1])
    ###################### Charging the list to correlate #####################
    
    ################## Creating the dictionaries to correlate #################
    RVpp_dict,Gonzalos_parameters_dict,best_models = {},{},{}
    for i in range(len(RVpp['HD2'])):
        RVpp_dict[RVpp['HD2'][i]] = [RVpp['vsini_final'][i],RVpp['Teff'][i],\
                                    RVpp['Rvpp_final'][i],RVpp['logLsp(Lsp=Teff/g)'][i]]
    
    for i in range(len(Gonzalos_parameters['Name'])):
        Gonzalos_parameters_dict[Gonzalos_parameters['Name'][i]] = \
            Gonzalos_parameters['Teff'][i],Gonzalos_parameters['vsini_F'][i],Gonzalos_parameters['vmac_F'][i]

    for i in range(len(best_models_txt['Name_HD'])):
        star,model = best_models_txt['Name_HD'][i],best_models_txt['fwlab'][i]
        best_models[star] = model
    ################## Creating the dictionaries to correlate #################
  
    Msun, Lsun = 2e30, 3.83e26 # S.I.
    G, sigma = 6.67e-11, 5.67e-8 # S.I.
    brk = False
    try:
        ################### Proobing if a star is consider ####################
        with open(path+'/short_stellar_data.csv', 'rt') as f:
             reader = csv.reader(f, delimiter=',') 
             for row in reader:
               if row[0] == star_id:
                   brk = True
        ################### Proobing if a star is consider ####################
    except:  # If the file does not existe
        ##################### Creating file + adding star #####################
        with open(path+'/short_stellar_data.csv', 'a') as f_object:
            star = {}
            star['STAR'] = star_id
            star['Teff'] = Gonzalos_parameters_dict[star_id][0]
            star['g'] = round(4*np.log10(float(Gonzalos_parameters_dict[star_id][0])*1e3) \
                - RVpp_dict[star_id][3] +np.log10(4*np.pi*sigma*G*Msun/Lsun) + 2,2)
            try: star['vsini'] = RVpp_dict[star_id][0]#Gonzalos_parameters_dict[star_id][1]
            except: star['vsini'] = Gonzalos_parameters_dict[star_id][1]
            
            star['vmac'] = Gonzalos_parameters_dict[star_id][2]
            star['RVpp'] = RVpp_dict[star_id][2]
            try: 
                if int(best_models[star_id][13:-3])<127: Q = '127'
                else: Q = best_models[star_id][13:-3]
                BM = best_models[star_id][:10]+'10Q'+Q+'0b10'
                star['FWmodel'] = BM
            except: star['FWmodel'] = None
            
            dictwriter_object = DictWriter(f_object,star.keys())
            dictwriter_object.writeheader()
            dictwriter_object.writerow(star)
            brk = True
            # Close the file object
            f_object.close()
        ##################### Creating file + adding star #####################
    if brk == False: # If the star is not in the file
        ########################### Adding the star ###########################
        star = {}
        star['STAR'] = star_id
        star['Teff'] = Gonzalos_parameters_dict[star_id][0]
        star['g'] = round(4*np.log10(float(Gonzalos_parameters_dict[star_id][0])*1e3) \
            - RVpp_dict[star_id][3] +np.log10(4*np.pi*sigma*G*Msun/Lsun) + 2,2)
        try: star['vsini'] = RVpp_dict[star_id][0]#Gonzalos_parameters_dict[star_id][1]
        except: star['vsini'] = Gonzalos_parameters_dict[star_id][1]
        star['vmac'] = Gonzalos_parameters_dict[star_id][2]
        star['RVpp'] = RVpp_dict[star_id][2]
        try: 
            if int(best_models[star_id][13:-3])<127: Q = '127'
            else: Q = best_models[star_id][13:-3]
            BM = best_models[star_id][:10]+'10Q'+Q+'0b10'
            star['FWmodel'] = BM
        except: star['FWmodel'] = None
        
        with open(path+'/short_stellar_data.csv', 'a') as f_object:
            dictwriter_object = DictWriter(f_object,star.keys())
            dictwriter_object.writerow(star)
         
            # Close the file object
            f_object.close()
        ########################### Adding the star ###########################
    else: print(star_id+' already existing')
    return None
################################ Adding stars #################################
# Simon-Diaz+10 added manually 

#%%
### Carneiro+19 ###
add_star('HD36512')
add_star('HD34078')
add_star('HD46202')
add_star('HD214680')
add_star('HD97848')
add_star('HD46966')
add_star('HD12993')
add_star('HD303311')
add_star('HD96715')
add_star('HD195592')
add_star('HD152249')
add_star('HD71304')
add_star('HD207198')
add_star('HD225160')
add_star('HD171589')
add_star('HD151515')
add_star('HD169582')
add_star('HD93222')
### Carneiro+19 ###

### Simon-Diaz+10 ###
add_star('HD37020')
add_star('HD36960')
add_star('HD37042')
add_star('HD36591')
add_star('HD36959')
### Simon-Diaz+10 ###

### Britavskiy ###
add_star('HD152590')
add_star('HD36486AaAb')
add_star('HD96622')
add_star('HD154643')
add_star('HD152424')
add_star('HD164438')
add_star('HD24431')
add_star('HD46202')
### Britavskiy ###

### vsini<45 & Teff<32 ###
add_star('HD55879')
add_star('HD195592')
add_star('HD68450')
add_star('HD118198')
add_star('HD101545A')
### vsini<45 & Teff<32 ###

### vsini<35 & 32<=Teff<=37 ###
add_star('HD57682')
add_star('HD206183')
add_star('HD54879')
add_star('BD+60499')
add_star('HD202214')
add_star('HD207538')
add_star('HD46149')
add_star('CPD-546791')
add_star('HD5005D')
### vsini<35 & 32<=Teff<=37 ###

### vsini<35 & Teff>37 ###
add_star('ALS12619')
add_star('HD44811')
add_star('HDE338916')
add_star('BD+61411')
add_star('HD242935A')
### vsini<35 & Teff>37 ###

### vsini<35 & Teff>37 ###
add_star('BD+60586')
add_star('BD+622078')
add_star('HD166546')
add_star('CPD-582611')
add_star('HD164492A')
add_star('HD35619')
add_star('HD47839AaAb')
add_star('HD218195A')
add_star('HD216898')
add_star('HD192001')
### vsini<35 & Teff>37 ###


### vsini<50 ###
add_star('HD16832'),
add_star('HD93027'),
add_star('HD193595'),
add_star('HDE227245'),
add_star('HDE344777'),
add_star('CPD-595634'),
add_star('HD305619'),
add_star('HD101190'),
add_star('HD42088')
### vsini<50 ###

'HD36486'
### vsini<60 ###
add_star('HD209975'),
add_star('HD99897'),
add_star('ALS8294'),
add_star('HD99546'),
add_star('HD46223'),
add_star('HD91824'),
add_star('BD+602635'),
add_star('HD305532'),
add_star('HD76341'),
add_star('HD237211'),
add_star('HD149038'),
add_star('HD36861'),
add_star('HD5005A'),
add_star('HD199579'),
add_star('HD76968'),
add_star('HD163800'),
add_star('HD162978'),
add_star('HD188209'),
add_star('HD101223'),
add_star('HD90273'),
add_star('HD105056'),
add_star('HD104565'),
add_star('HD191978'),
add_star('HD305523'),
add_star('HD93843'),
add_star('HD93128'),
add_star('HD152405'),
add_star('HD114737')
### vsini<60 ###


### Alex's TFM ###
add_star('HD85953'),
add_star('HD3360'),
add_star('HD35039'),
add_star('HD42690'),
add_star('HD157056'),
add_star('HD886'),
add_star('HD253049'),
add_star('HD60553'),
add_star('HD223128'),
add_star('HD829'),
add_star('HD77770'),
add_star('HD162374'),
add_star('HD174298'),
add_star('HD64365'),
add_star('HD158304'),
add_star('HD37356'),
add_star('HD220172'),
add_star('HD35299'),
add_star('HD197512'),
add_star('HD35337'),
add_star('HD43112'),
add_star('CPD573524B'),
add_star('HD201795'),
add_star('HD34816'),
add_star('HD36822'),
add_star('HD119069'),
add_star('HD46328'),
add_star('HD166523'),
add_star('HD59882'),
add_star('HD255134'),
add_star('HD34989'),
add_star('HD37209'),
add_star('HD54025'),
add_star('HD129557'),
add_star('HD251847'),
add_star('HD201638'),
add_star('HD211880'),
add_star('HD164704'),
add_star('HD165016'),
add_star('HD192039'),
add_star('HD166539'),
add_star('HD161789')
### Alex's TFM ###


### vsini<70 ###
add_star('HD46150'),
add_star('HD91572'),
add_star('HD93146'),
add_star('HD135591'),
add_star('BD+62424'),
add_star('HDE227018'),
add_star('HDE227465'),
add_star('HD168504'),
add_star('HD186980'),
add_star('HD173010'),
add_star('HD156154'),
add_star('HD123008'),
add_star('HD164794'),
add_star('HD152233'),
add_star('HD155756'),
add_star('HD218915'),
add_star('HD154368'),
add_star('HD156738'),
add_star('HD344784'),
add_star('HD152003'),
add_star('HD190864'),
add_star('TRUMPLER14-9'),
add_star('HD93129'),
add_star('HD305525'),
add_star('CPD-472963'),
add_star('HD34656'),
add_star('HD225146'),
add_star('HD256725'),
add_star('HD151018'),
add_star('HD168076'),
add_star('HD322417'),
add_star('HD188001'),
add_star('HD319699'),
add_star('HD164019')
 ### vsini<70 ###

### vsini<80 ###
add_star('HD210809'),
add_star('HD93250'),
add_star('HD97253'),
add_star('HD101298'),
add_star('HD167659'),
add_star('HD167264'),
add_star('HD151804'),
add_star('HD96946'),
add_star('CYGOB2-4'),
add_star('HD69464'),
add_star('HD193514'),
add_star('HD152723'),
add_star('BD-114586'),
add_star('CYGOB2-7'),
add_star('CYGOB2-11'),
add_star('HD64568'),
add_star('HD163758'),
add_star('HD94963'),
add_star('HD15629'),
add_star('HD46573'),
add_star('HD201345')
### vsini<80 ###

### vsini<90 ###
add_star('HD101413'),
add_star('HD15570'),
add_star('HD93161B'),
add_star('HD192639'),
add_star('HD152247'),
add_star('HD152623'),
add_star('HD168941'),
add_star('ALS4880'),
add_star('HD209339'),
add_star('HD113904B'),
add_star('HD303308AB'),
add_star('ALS7833'),
add_star('HD148546'),
add_star('HD93632'),
add_star('HD75222'),
add_star('HD303492'),
add_star('HD326329'),
add_star('HD202124'), 
add_star('HD189957'),
add_star('HD242926'),
add_star('HD173783')        
### vsini<90 ###

### vsini<100 ###
add_star('HD190429'),
add_star('BD+391328'),
add_star('HD152147'),
add_star('HD110360'),
add_star('HD47432'),
add_star('HD13022'),
add_star('ALS15210'),
add_star('BD-134927'),
add_star('HD226868')
### vsini<100 ###

### vsini<110 ###
add_star('HD120521'), 
add_star('HD17603'), 
add_star('HD191781'), 
add_star('HD298429'), 
add_star('HD168112')
### vsini<110 ###

### vsini<120 ###
add_star('HD93204'),
add_star('HD38666'),
add_star('CYGOB2-9'),
add_star('HD37743'),
add_star('HD30614'),
add_star('HD14947'),
add_star('HD157857'),
add_star('HD159176'),
add_star('HD116852'),
add_star('BD+60498'),
add_star('HD125241')
### vsini<120 ###

### vsini<125 ###
add_star('HD93160'),
add_star('HD14633'),
add_star('HD12323'),
add_star('HD10125'),
add_star('HD37742'),
add_star('HD154811'),
add_star('CPD-592551'),
add_star('HD112244')
### vsini<125###

### vsini<140 ###
add_star('CPD-592600'),
add_star('HD167633'),
add_star('HD48279A'),
add_star('HD18409'),
add_star('HD101191'),
add_star('CPD-262716')
### vsini<140 ###

### vsini<140 ###
add_star('HD105627')
add_star('HD75211')
add_star('HD73882')
add_star('HD117797')
### vsini<140 ###
#%% Trying to change vsini to Mijail's values 
data_base = pd.read_excel(path+'/FASTWIND/star_analysis/FINAL_TABLE_TO_PLOTS_mar2023_ALL_parameters.xlsx',usecols=[1,8,10,12,15])
stellar_data = pd.read_csv(path+"/short_stellar_data.csv")

for star in stellar_data['STAR']:
    if star in data_base['HD2'].values:
        #print(stellar_data['vsini'][stellar_data['STAR']==star].values[0],\
        #      data_base['vsini_final'][data_base['HD2']==star].values[0])
        stellar_data['vsini'][stellar_data['STAR']==star] = \
              data_base['vsini_final'][data_base['HD2']==star].values
    else:
        print(star)
        
stellar_data['vsini']
stellar_data.to_csv(path+"/short_stellar_data.csv",index=False)