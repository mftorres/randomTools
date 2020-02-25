# python
# code by M.F.Torres.J, 2020
# Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)


####################################################################################################################################
# go to http://wcsp.science.kew.org/reportbuilder.do and generate your list
# save the list as HTML: cursor over the list, then right click > inspect > search for main ID='main' > righ click > copy element
# open a new document, paste the element and save as HTML

# file saved with UTF-8
# removed symbol *

# replace 'Native_species_to_btocountries.html' with your file
####################################################################################################################################

import pandas as pd
import re

kew_html_path='./Native_species_to_btocountries.html'
species_areas_dict={}
species_synon_dict={}
for line in open(kew_html_path,'r',encoding="utf8"):
    # detects the species count line and doesn't take it as species line (it is genus line)
    speciescountmatch=re.search('rptDetNameSpcCount',line)
    genuscountmatch=re.search('Arecaceae',line)
    # skips lines with irrelevant information
    if speciescountmatch or genuscountmatch:
        pass
    else:
        detnamematch=re.search('rptDetName',line)
        if detnamematch:
            # fetching information
            speciesmatch=re.search('<i><b>([A-Z][a-z]+)</b><b> ([a-z-]+)</b></i>',line) # note the space
            areamatch=re.search('"rptDetNameDescription">(.+?)</p>',line)
            synonmatch=re.findall('("rptDetNameSynonyms"> <i>[A-Z][a-z]+</i><i> [a-z-]+</i>.+</p>)',line)
            
            # asigning
            species='%s_%s'%(speciesmatch.group(1),speciesmatch.group(2)) if speciesmatch else None
            species_areas_dict[species]=re.findall('([A-Z][A-Z][A-Z])+',areamatch.group(0)) if areamatch else None
            
            synonyms=re.findall('([A-Z][a-z]+</i><i> [a-z-]+)</i> ',line)# if synonmatch else None
            species_synon_dict[species]=[x.replace('</i><i> ','_') for x in synonyms] if synonmatch else 'none'
            
data=pd.DataFrame()

data['species']=list(species_areas_dict.keys())
data['areas']=data['species'].map(species_areas_dict)
data['synonyms']=data['species'].map(species_synon_dict)
data.to_csv('./Plants_countries_synonyms.csv')
