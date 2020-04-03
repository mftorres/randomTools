  
# python
# code by M.F.Torres.J, 2020
# Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)

####################################################################################################################################
# The script stores the HTML body and pieces from the International Union for Conservation of Nature
# it searchers for specific text within the HTML. It was written as a quick example on how to do a task someone else needed.

# needs opening with every website visited. Pretty aggressive but...
# download lates version of ChromeDriver that suits your chrome version
# at https://chromedriver.storage.googleapis.com/index.html?path=77.0.3865.40/ for Chrome v81

## URL looking like
## https://www.iucnredlist.org/species/14519/97215090
## https://www.iucnredlist.org/species/18597/46364962
#####


import pandas as pd
import urllib3
from bs4 import BeautifulSoup
import requests
import re
from selenium import webdriver
import time
import sys
print('Python version: %s'%(sys.version))

driver=webdriver.Chrome('./chromedriver.exe')
driver.get('https://www.iucnredlist.org/species/7951/45171204')
driver.find_elements_by_class_name('layout-card--split__minor')
html=driver.page_source

# each entry should be '[0-9]+/[0-9]+'
page_list=['',''] # fill up with the numbers to add in the URL

dic_sppyears={}
for page in page_list:
    driver=webdriver.Chrome('C:/Users/xtorrm/Desktop/chromedriver_win32/chromedriver.exe')
    driver.get('https://www.iucnredlist.org/species/%s'%(page))
    driver.find_elements_by_class_name('layout-card--split__minor')
    html=driver.page_source
    # print(html)


    # <title>Neofelis nebulosa</title>

    # creates an empty dictionary with species name as keys and generation times as values
    for line in html.split('\n'):
        # starts assuming there's no data
        years_data='absent'
        species_match=re.search('title\>([A-Za-z]+ [a-z]+).+\<\/title',line)
        generation_match=re.search('Generation length \(years\)</h3><p class="card__data card__data--std card__data--accent">([0-9].+) (years)</p>',line)
        # species name should appear before than generation time
        # if species name is found, asigns to variable
        if species_match:
            spp_name=species_match.group(1)
        # continues down the stream searching for generation time, if True, creates the accession in the dictionary
        elif generation_match:
            dic_sppyears[spp_name]=generation_match.group(1)
            years_data='present'
        # continues through lines after finding the data in the line
        # passes until end of html
        elif years_data == 'present':
            pass
        # goes to the end and asigns the key species name a value of no information
        else:
            dic_sppyears[spp_name]='no_data'
        time.sleep(5) # 1 second I think?

# print the dic:
print(dic_sppyears)
## access dictionary by
# dic_sppyears['species_name']=information
