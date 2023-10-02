#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra
import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
from pandas import DataFrame
import plotly
import plotly.express as px
import plotly.graph_objs as go
import sklearn
import seaborn as sns


# In[2]:


comparison = pd.DataFrame()
get_ipython().run_line_magic('store', '-r _59m_absolute_fluxes')
get_ipython().run_line_magic('store', '-r heya8_absolute_fluxes')
get_ipython().run_line_magic('store', '-r ov56_absolute_fluxes')
get_ipython().run_line_magic('store', '-r caov3_absolute_fluxes')
get_ipython().run_line_magic('store', '-r cov318_absolute_fluxes')
get_ipython().run_line_magic('store', '-r oaw28_absolute_fluxes')
comparison['reaction'] = _59m_absolute_fluxes['reaction']
comparison['59m_absolute_flux'] = _59m_absolute_fluxes['flux']
comparison['heya8_absolute_flux'] = heya8_absolute_fluxes['flux']
comparison['ov56_absolute_flux'] = ov56_absolute_fluxes['flux']
comparison['caov3_absolute_flux'] = caov3_absolute_fluxes['flux']
comparison['cov318_absolute_flux'] = cov318_absolute_fluxes['flux']
comparison['oaw28_absolute_flux'] = oaw28_absolute_fluxes['flux']
comparison


# In[46]:


for r in range(len(comparison['reaction'])):
    if comparison.iloc[r,0] == 'HMR_4501':
        print(comparison.iloc[r,:])


# In[28]:


(1.805663+3.922173)/3


# In[29]:


0.071612/3


# In[103]:


for r in model.reactions: 
    print(r.id)


# In[623]:


get_ipython().run_line_magic('store', 'comparison')


# In[3]:


get_ipython().run_line_magic('store', '-r comparison_5')


# In[4]:


comparison_5


# In[7]:


zero_rows = []
for n in range(len(comparison['reaction'])):
    if comparison.iloc[n,1] == 0:
        if comparison.iloc[n,2] == 0:
            if comparison.iloc[n,3] == 0:
                if comparison.iloc[n,4] == 0:
                    if comparison.iloc[n,5] == 0:
                        if comparison.iloc[n,6] == 0:
                            zero_rows.append(n)
print('proportion of rows where all fluxes are zero:',(len(zero_rows)/len(comparison['reaction'])))


# In[9]:


(1-(len(zero_rows)/len(comparison['reaction'])))*100


# In[10]:


comparison_2 = comparison.drop(comparison.index[zero_rows])


# In[11]:


len(comparison_2)


# In[15]:


comparison_2


# In[28]:


105/244


# In[13]:


difference = []
ratio = []
for n in range(len(comparison_2['reaction'])):
    lg = comparison_2.iloc[n,1:4].to_list()
    hg = comparison_2.iloc[n,4:].to_list()
    avg_low = ((lg[0] + lg[1] + lg[2])/3)
    avg_high = ((hg[0] + hg[1] + hg[2])/3)
    if avg_high > avg_low:
        if avg_low < avg_high*0.9:
            difference.append('sig_diff')
        else:
            difference.append('no_sig_diff')
    if avg_low > avg_high:
        if avg_high < avg_low*0.9:
            difference.append('sig_diff')
        else:
            difference.append('no_sig_diff')
    if avg_low == avg_high:
        difference.append('no_sig_diff')
comparison_2['difference'] = difference

sig_rows = []
for n in range(len(comparison_2['reaction'])): #remove rows with less than 10% diff (calculated above)
    if comparison_2.iloc[n,7] == 'sig_diff':
        sig_rows.append(n)
insig_rows = []
for n in range(len(comparison_2['reaction'])):
    if comparison_2.iloc[n,7] == 'no_sig_diff':
        insig_rows.append(n)
        
comparison_3 = comparison_2.drop(comparison_2.index[insig_rows])
comparison_3


# In[17]:


len(comparison_3)/len(comparison_2)


# In[77]:


insig_rows_2 = []       
for n in range(len(comparison_3['reaction'])): #remove rows with more than 67% zeros.
    number_zeros = []
    if comparison_3.iloc[n,1] == 0:
        number_zeros.append(n)
    if comparison_3.iloc[n,2] == 0:
        number_zeros.append(n)
    if comparison_3.iloc[n,3] == 0:
        number_zeros.append(n)
    if comparison_3.iloc[n,4] == 0:
        number_zeros.append(n)
    if comparison_3.iloc[n,5] == 0:
        number_zeros.append(n)
    if comparison_3.iloc[n,6] == 0:
        number_zeros.append(n)
    if len(number_zeros) >=5:
        insig_rows_2.append(n)
        
comparison_4 = comparison_3.drop(comparison_3.index[insig_rows_2])
comparison_4


# In[79]:


low_rows = []
for n in range(len(comparison_4['reaction'])):
    if comparison_4.iloc[n,1] < 0.5:
        if comparison_4.iloc[n,2] < 0.5:
            if comparison_4.iloc[n,3] < 0.5:
                if comparison_4.iloc[n,4] < 0.5:
                    if comparison_4.iloc[n,5] < 0.5:
                        if comparison_4.iloc[n,6] < 0.5:
                            low_rows.append(n)

comparison_5 = comparison_4.drop(comparison_4.index[low_rows])
comparison_5


# In[137]:


#find the flux through a reaction for each model
reaction = 'HMR_6367'
for n in range(len(comparison['reaction'])):
    if comparison.iloc[n,0] == reaction:
        print(comparison.iloc[n,:])


# In[102]:


from mewpy.simulation import solvers
from mewpy.simulation import set_default_solver
set_default_solver('glpk')

from cobra.io.sbml import read_sbml_model
model = read_sbml_model('Human-GEM-annotated.xml')


# In[626]:


model


# In[107]:


model.reactions.get_by_id(comparison_5.iloc[0,0]).gene_reaction_rule


# In[102]:


significant_reactions = comparison_5.iloc[:,0].to_list()
significant_reactions


# In[103]:


overall_rule_list = []
for r in significant_reactions:
    rule_list = []
    rule = model.reactions.get_by_id(r).gene_reaction_rule
    rule_list.append(rule.split())
    rule_genes = []
    for n in rule_list:
        for num in n:
            if 'ENS' in num:
                rule_genes.append(num)
    overall_rule_list.append(rule_genes)


# In[629]:


(1-92/336)*100


# In[105]:


comparison_5['rules'] = overall_rule_list
comparison_5


# In[109]:


comparison_5.to_csv('comparison_5.csv')


# In[111]:


comparison_5 = pd.read_csv(r'comparison_5.csv')
comparison_5


# In[112]:


comparison_5.keys()


# In[116]:


comparison_5 = comparison_5.sort_values(by=['rules'])
comparison_5


# In[156]:


rules_1 = []
for n in comparison_5['rules']:
    rules_1.append(n)
    
rules_2 = [x for x in rules_1 if str(x) != 'nan']

from collections import Counter
counter = Counter(rules_2)


# In[856]:


for num in range(len(comparison_5['reaction'])):
    if 'HMR_4802' in comparison_5.iloc[num,1]:
        print(comparison_5.iloc[num,1])
        print(comparison_5.iloc[num,2:8])


# In[854]:


n_reactions = []
for n in range(len(comparison_5['rules'])):
    if 'ENSG00000112759' == comparison_5.iloc[n,9]:
        n_reactions.append(comparison_5.iloc[n,1])
for n in n_reactions:
    for num in range(len(comparison_5['reaction'])):
        if comparison_5.iloc[num,1] == n:
            print(n)
            print(comparison_5.iloc[num,2:8])


# In[850]:


comparison_5


# In[157]:


counter.most_common()


# In[154]:


rules_3 = []
for n in comparison_5.iloc[:,10]:
    rules_3.append(n)
    
rules_4 = [x for x in rules_3 if str(x) != 'nan']

from collections import Counter
counter_2 = Counter(rules_4)


# In[155]:


counter_2.most_common()


# In[159]:


rules_5 = []
for n in comparison_5.iloc[:,11]:
    rules_5.append(n)
    
rules_6 = [x for x in rules_5 if str(x) != 'nan']

counter_3 = Counter(rules_6)
counter_3.most_common()


# # CPT1A: a regulator of lipid metabolism

# In[112]:


#upload dictionaries for each cell line so I can scan gene expression
get_ipython().run_line_magic('store', '-r caov3_gene_exp_dict')
get_ipython().run_line_magic('store', '-r cov318_gene_exp_dict')
get_ipython().run_line_magic('store', '-r oaw28_gene_exp_dict')
get_ipython().run_line_magic('store', '-r heya8_gene_exp_dict')
get_ipython().run_line_magic('store', '-r ov56_gene_exp_dict')
get_ipython().run_line_magic('store', '-r _59m_gene_exp_dict')


# In[129]:


print(_59m_gene_exp_dict['ENSG00000111775'])
print(heya8_gene_exp_dict['ENSG00000111775'])
print(ov56_gene_exp_dict['ENSG00000111775'])
print(caov3_gene_exp_dict['ENSG00000111775'])
print(cov318_gene_exp_dict['ENSG00000111775'])
print(oaw28_gene_exp_dict['ENSG00000111775'])


# In[123]:


# create a dataset
y1 = (_59m_gene_exp_dict['ENSG00000074800']+_59m_gene_exp_dict['ENSG00000111674']+_59m_gene_exp_dict['ENSG00000108515'])
y2 = (heya8_gene_exp_dict['ENSG00000074800']+heya8_gene_exp_dict['ENSG00000111674']+heya8_gene_exp_dict['ENSG00000108515'])
y3 = (ov56_gene_exp_dict['ENSG00000074800']+ov56_gene_exp_dict['ENSG00000111674']+ov56_gene_exp_dict['ENSG00000108515'])
y4 = (caov3_gene_exp_dict['ENSG00000074800']+caov3_gene_exp_dict['ENSG00000111674']+caov3_gene_exp_dict['ENSG00000108515'])
y5 = (cov318_gene_exp_dict['ENSG00000074800']+cov318_gene_exp_dict['ENSG00000111674']+cov318_gene_exp_dict['ENSG00000108515'])
y6 = (oaw28_gene_exp_dict['ENSG00000074800']+oaw28_gene_exp_dict['ENSG00000111674']+oaw28_gene_exp_dict['ENSG00000108515'])

x = [y1,y2,y3,y4,y5,y6]
bars = ('59M','HEYA8','OV56','CAOV3','COV318','OAW28')
x_pos = np.arange(len(bars))

# Create bars with different colors
plt.bar(x_pos, x, color=[plotly.colors.qualitative.Dark24[3], plotly.colors.qualitative.G10[8]
                         , plotly.colors.qualitative.Light24[0], plotly.colors.qualitative.G10[7], 
                         plotly.colors.qualitative.Plotly[7],plotly.colors.qualitative.D3[2]], edgecolor = 'black')

plt.xticks(x_pos, bars)
plt.setp(plt.gca().get_xticklabels(), horizontalalignment='center')
plt.ylabel('CCLE; normalised gene expression')
plt.title('gene expression of enolase')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Show graph
plt.savefig('ENO_summed.png', dpi=300)
plt.show()


# In[125]:


y1 = (_59m_gene_exp_dict['ENSG00000074800'])
y2 = (heya8_gene_exp_dict['ENSG00000074800'])
y3 = (ov56_gene_exp_dict['ENSG00000074800'])
y4 = (caov3_gene_exp_dict['ENSG00000074800'])
y5 = (cov318_gene_exp_dict['ENSG00000074800'])
y6 = (oaw28_gene_exp_dict['ENSG00000074800'])

x = [y1,y2,y3,y4,y5,y6]
bars = ('59M','HEYA8','OV56','CAOV3','COV318','OAW28')
x_pos = np.arange(len(bars))

# Create bars with different colors
plt.bar(x_pos, x, color=[plotly.colors.qualitative.Dark24[3], plotly.colors.qualitative.G10[8]
                         , plotly.colors.qualitative.Light24[0], plotly.colors.qualitative.G10[7], 
                         plotly.colors.qualitative.Plotly[7],plotly.colors.qualitative.D3[2]], edgecolor = 'black')

plt.xticks(x_pos, bars)
plt.setp(plt.gca().get_xticklabels(), horizontalalignment='center')
plt.ylabel('CCLE; normalised gene expression')
plt.title('gene expression of ENO1')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Show graph
plt.savefig('ENO1.png', dpi=300)
plt.show()


# In[126]:


y1 = (_59m_gene_exp_dict['ENSG00000111674'])
y2 = (heya8_gene_exp_dict['ENSG00000111674'])
y3 = (ov56_gene_exp_dict['ENSG00000111674'])
y4 = (caov3_gene_exp_dict['ENSG00000111674'])
y5 = (cov318_gene_exp_dict['ENSG00000111674'])
y6 = (oaw28_gene_exp_dict['ENSG00000111674'])

x = [y1,y2,y3,y4,y5,y6]
bars = ('59M','HEYA8','OV56','CAOV3','COV318','OAW28')
x_pos = np.arange(len(bars))

# Create bars with different colors
plt.bar(x_pos, x, color=[plotly.colors.qualitative.Dark24[3], plotly.colors.qualitative.G10[8]
                         , plotly.colors.qualitative.Light24[0], plotly.colors.qualitative.G10[7], 
                         plotly.colors.qualitative.Plotly[7],plotly.colors.qualitative.D3[2]], edgecolor = 'black')

plt.xticks(x_pos, bars)
plt.setp(plt.gca().get_xticklabels(), horizontalalignment='center')
plt.ylabel('CCLE; normalised gene expression')
plt.title('gene expression of ENO2')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Show graph
plt.savefig('ENO2.png', dpi=300)
plt.show()


# In[124]:


y1 = (_59m_gene_exp_dict['ENSG00000108515'])
y2 = (heya8_gene_exp_dict['ENSG00000108515'])
y3 = (ov56_gene_exp_dict['ENSG00000108515'])
y4 = (caov3_gene_exp_dict['ENSG00000108515'])
y5 = (cov318_gene_exp_dict['ENSG00000108515'])
y6 = (oaw28_gene_exp_dict['ENSG00000108515'])

x = [y1,y2,y3,y4,y5,y6]
bars = ('59M','HEYA8','OV56','CAOV3','COV318','OAW28')
x_pos = np.arange(len(bars))

# Create bars with different colors
plt.bar(x_pos, x, color=[plotly.colors.qualitative.Dark24[3], plotly.colors.qualitative.G10[8]
                         , plotly.colors.qualitative.Light24[0], plotly.colors.qualitative.G10[7], 
                         plotly.colors.qualitative.Plotly[7],plotly.colors.qualitative.D3[2]], edgecolor = 'black')

plt.xticks(x_pos, bars)
plt.setp(plt.gca().get_xticklabels(), horizontalalignment='center')
plt.ylabel('CCLE; normalised gene expression')
plt.title('gene expression of ENO3')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Show graph
plt.savefig('ENO3.png', dpi=300)
plt.show()


# In[170]:


lg_mean = (y1+y2+y3/3)
hg_mean = (y4+y5+y6/3)
print('lg mean:',lg_mean)
print('hg mean:',hg_mean)

lg_std = np.std([y1,y2,y3])
hg_std = np.std([y4,y5,y6])
print('lg std:',lg_std)
print('hg std:',hg_std)


# In[180]:


print('upper limit of lg:',lg_mean+lg_std,'lower limit of hg:',hg_mean-hg_std)

if (lg_mean+lg_std)<(hg_mean-hg_std):
    print('significant difference between expressions')


# In[184]:


groups = ['low-grade', 'high-grade']
x_pos = np.arange(len(groups))
means = [lg_mean, hg_mean]
error = [lg_std, hg_std]

fig, ax = plt.subplots()
ax.bar(x_pos, means, yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10, color=['red','green'])
ax.set_ylabel('CCLE; normalised gene expression')
ax.set_xticks(x_pos)
ax.set_xticklabels(groups)
ax.set_title('gene expression of CPT1A')
ax.yaxis.grid(True)

plt.tight_layout()
plt.savefig('cpt1a_expression.png')
plt.show()


# In[221]:


cpt1a_reactions = []
cpt1a_iloc_indexes = []
for n in range(len(comparison_5['reaction'])):
    for num in comparison_5.iloc[n,:]:
        if num == 'ENSG00000110090':
            print(comparison_5.iloc[n,1],':',n,model.reactions.get_by_id(comparison_5.iloc[n,1]).name)
            cpt1a_reactions.append(comparison_5.iloc[n,1])
            cpt1a_iloc_indexes.append(n)


# In[249]:


cpt1b_reactions = []
cpt1b_iloc_indexes = []
for n in range(len(comparison_5['reaction'])):
    for num in comparison_5.iloc[n,:]:
        if num == 'ENSG00000205560':
            print(comparison_5.iloc[n,1],':',n,model.reactions.get_by_id(comparison_5.iloc[n,1]).name)
            cpt1b_reactions.append(comparison_5.iloc[n,1])
            cpt1b_iloc_indexes.append(n)


# In[250]:


cpt1c_reactions = []
cpt1c_iloc_indexes = []
for n in range(len(comparison_5['reaction'])):
    for num in comparison_5.iloc[n,:]:
        if num == 'ENSG00000169169':
            print(comparison_5.iloc[n,1],':',n,model.reactions.get_by_id(comparison_5.iloc[n,1]).name)
            cpt1c_reactions.append(comparison_5.iloc[n,1])
            cpt1c_iloc_indexes.append(n)


# In[199]:


for n in range(len(comparison_5['reaction'])):
    if comparison_5.iloc[n,1] in cpt1a_reactions:
        print(comparison_5.iloc[n,1:8])
        print('\n')


# In[200]:


#take means and plot 2 bars with std
cpt1a_reactions


# In[209]:


comparison_5.iloc[123,2:5]


# In[226]:


r0437_lg_mean = sum(comparison_5.iloc[123,2:5].to_list())/3
r0437_lg_std = np.std(comparison_5.iloc[123,2:5].to_list())
print(r0437_lg_mean)
print(r0437_lg_std)
r0437_hg_mean = sum(comparison_5.iloc[123,5:8].to_list())/3
r0437_hg_std = np.std(comparison_5.iloc[123,5:8].to_list())
print(r0437_hg_mean)
print(r0437_hg_std)
if r0437_lg_mean<r0437_hg_mean:
    if (r0437_lg_mean+r0437_lg_std)<(r0437_hg_mean-r0437_hg_std):
        print('significant difference between expressions')
    else:
        print('no significant difference between expressions')
if r0437_lg_mean>r0437_hg_mean:
    if (r0437_lg_mean-r0437_lg_std)>(r0437_hg_mean+r0437_hg_std):
        print('significant difference between expressions')
    else:
        print('no significant difference between expressions')
    
HMR_0159_lg_mean = sum(comparison_5.iloc[124,2:5].to_list())/3
HMR_0159_lg_std = np.std(comparison_5.iloc[124,2:5].to_list())
print(HMR_0159_lg_mean)
print(HMR_0159_lg_std)
HMR_0159_hg_mean = sum(comparison_5.iloc[124,5:8].to_list())/3
HMR_0159_hg_std = np.std(comparison_5.iloc[124,5:8].to_list())
print(HMR_0159_hg_mean)
print(HMR_0159_hg_std)
if HMR_0159_lg_mean<HMR_0159_hg_mean:
    if (HMR_0159_lg_mean+HMR_0159_lg_std)<(HMR_0159_hg_mean-HMR_0159_hg_std):
        print('significant difference between expressions')
    else:
        print('no significant difference between expressions')
if HMR_0159_lg_mean>HMR_0159_hg_mean:
    if (HMR_0159_lg_mean-HMR_0159_lg_std)>(HMR_0159_hg_mean+HMR_0159_hg_std):
        print('significant difference between expressions')
    else:
        print('no significant difference between expressions')
    
HMR_2626_lg_mean = sum(comparison_5.iloc[125,2:5].to_list())/3
HMR_2626_lg_std = np.std(comparison_5.iloc[125,2:5].to_list())
print(HMR_2626_lg_mean)
print(HMR_2626_lg_std)
HMR_2626_hg_mean = sum(comparison_5.iloc[125,5:8].to_list())/3
HMR_2626_hg_std = np.std(comparison_5.iloc[125,5:8].to_list())
print(HMR_2626_hg_mean)
print(HMR_2626_hg_std)
if HMR_2626_lg_mean<HMR_2626_hg_mean:
    if (HMR_2626_lg_mean+HMR_2626_lg_std)<(HMR_2626_hg_mean-HMR_2626_hg_std):
        print('significant difference between expressions')
    else:
        print('no significant difference between expressions')
if HMR_2626_lg_mean>HMR_2626_hg_mean:
    if (HMR_2626_lg_mean-HMR_2626_lg_std)>(HMR_2626_hg_mean+HMR_2626_hg_std):
        print('significant difference between expressions')
    else:
        print('no significant difference between expressions')
    
HMR_2591_lg_mean = sum(comparison_5.iloc[126,2:5].to_list())/3
HMR_2591_lg_std = np.std(comparison_5.iloc[126,2:5].to_list())
print(HMR_2591_lg_mean)
print(HMR_2591_lg_std)
HMR_2591_hg_mean = sum(comparison_5.iloc[126,5:8].to_list())/3
HMR_2591_hg_std = np.std(comparison_5.iloc[126,5:8].to_list())
print(HMR_2591_hg_mean)
print(HMR_2591_hg_std)
if HMR_2591_lg_mean<HMR_2591_hg_mean:
    if (HMR_2591_lg_mean+HMR_2591_lg_std)<(HMR_2591_hg_mean-HMR_2591_hg_std):
        print('significant difference between expressions')
    else:
        print('no significant difference between expressions')
if HMR_2591_lg_mean>HMR_2591_hg_mean:
    if (HMR_2591_lg_mean-HMR_2591_lg_std)>(HMR_2591_hg_mean+HMR_2591_hg_std):
        print('significant difference between expressions')
    else:
        print('no significant difference between expressions')
    
HMR_2636_lg_mean = sum(comparison_5.iloc[127,2:5].to_list())/3
HMR_2636_lg_std = np.std(comparison_5.iloc[127,2:5].to_list())
print(HMR_2636_lg_mean)
print(HMR_2636_lg_std)
HMR_2636_hg_mean = sum(comparison_5.iloc[127,5:8].to_list())/3
HMR_2636_hg_std = np.std(comparison_5.iloc[127,5:8].to_list())
print(HMR_2636_hg_mean)
print(HMR_2636_hg_std)
if HMR_2636_lg_mean<HMR_2636_hg_mean:
    if (HMR_2636_lg_mean+HMR_2636_lg_std)<(HMR_2636_hg_mean-HMR_2636_hg_std):
        print('significant difference between expressions')
    else:
        print('no significant difference between expressions')
if HMR_2636_lg_mean>HMR_2636_hg_mean:
    if (HMR_2636_lg_mean-HMR_2636_lg_std)>(HMR_2636_hg_mean+HMR_2636_hg_std):
        print('significant difference between expressions')
    else:
        print('no significant difference between expressions')


# In[222]:


print(cpt1a_reactions)
print(cpt1a_iloc_indexes)


# In[238]:


#HMR_0159 and HMR_2591 have significant difference between LG and HG
groups = ['LG HMR_0159', 'HG HMR_0159', 'LG HMR_2591', 'HG HMR_2591']
x_pos = np.arange(len(groups))
means = [HMR_0159_lg_mean, HMR_0159_hg_mean, HMR_2591_lg_mean, HMR_2591_hg_mean]
error = [HMR_0159_lg_std, HMR_0159_hg_std, HMR_2591_lg_std, HMR_2591_hg_std]

fig, ax = plt.subplots()
ax.bar(x_pos, means, yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10, color=[plotly.colors.qualitative.Dark24[3],plotly.colors.qualitative.Light24[8], plotly.colors.qualitative.Dark24[3],plotly.colors.qualitative.Light24[8]])
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
ax.set_ylabel('reaction rate (mmol/gDW/hour)')
ax.set_xticks(x_pos)
ax.set_xticklabels(groups)
ax.set_title('flux through CPT1A-regulated reactions')
ax.yaxis.grid(True)

plt.tight_layout()
plt.savefig('cpt1a_reaction_flux.png')
plt.show()


# In[239]:


y1 = (_59m_gene_exp_dict['ENSG00000205560']+_59m_gene_exp_dict['ENSG00000110090']+_59m_gene_exp_dict['ENSG00000169169'])
y2 = (heya8_gene_exp_dict['ENSG00000205560']+heya8_gene_exp_dict['ENSG00000110090']+heya8_gene_exp_dict['ENSG00000169169'])
y3 = (ov56_gene_exp_dict['ENSG00000205560']+ov56_gene_exp_dict['ENSG00000110090']+ov56_gene_exp_dict['ENSG00000169169'])
y4 = (caov3_gene_exp_dict['ENSG00000205560']+caov3_gene_exp_dict['ENSG00000110090']+caov3_gene_exp_dict['ENSG00000169169'])
y5 = (cov318_gene_exp_dict['ENSG00000205560']+cov318_gene_exp_dict['ENSG00000110090']+cov318_gene_exp_dict['ENSG00000169169'])
y6 = (oaw28_gene_exp_dict['ENSG00000205560']+oaw28_gene_exp_dict['ENSG00000110090']+oaw28_gene_exp_dict['ENSG00000169169'])


# In[246]:


lg_mean = (y1+y2+y3/3)
hg_mean = (y4+y5+y6/3)
print('lg mean:',lg_mean)
print('hg mean:',hg_mean)

lg_std = np.std([y1,y2,y3])
hg_std = np.std([y4,y5,y6])
print('lg std:',lg_std)
print('hg std:',hg_std)


# In[248]:


groups = ['low-grade', 'high-grade']
x_pos = np.arange(len(groups))
means = [lg_mean, hg_mean]
error = [lg_std, hg_std]

fig, ax = plt.subplots()
ax.bar(x_pos, means, yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10, color=['red','green'])
ax.set_ylabel('CCLE; normalised gene expression')
ax.set_xticks(x_pos)
ax.set_xticklabels(groups)
ax.set_title('summed gene expression of CPT1A/B/C')
ax.yaxis.grid(True)

plt.tight_layout()
plt.savefig('cpt1_expression.png')
plt.show()


# In[256]:


comparison_5.loc[[111,112]]


# In[263]:


subset = comparison_5.loc[[111,112]]

y1 = subset['59m_absolute_flux']
y2 = subset['heya8_absolute_flux']
y3 = subset['ov56_absolute_flux']
y4 = subset['caov3_absolute_flux']
y5 = subset['cov318_absolute_flux']
y6 = subset['oaw28_absolute_flux']

N=len(subset)
ind = np.arange(N) 
width = 0.1      
plt.bar(ind, y1, width, label='59m_fluxes', color=plotly.colors.qualitative.Dark24[3])
plt.bar(ind + width, y2, width,
    label='heya8_fluxes', color=plotly.colors.qualitative.G10[8])
plt.bar(ind + width + width, y3, width,
    label='ov56_fluxes', color=plotly.colors.qualitative.Light24[0])
plt.bar(ind + 3*width, y4, width,
    label='caov3_fluxes', color=plotly.colors.qualitative.G10[7])
plt.bar(ind + 4*width, y5, width,
    label='cov318_fluxes', color=plotly.colors.qualitative.Plotly[7])
plt.bar(ind + 5*width, y6, width,
    label='oaw28_fluxes', color=plotly.colors.qualitative.D3[2])

reaction_ids = subset['reaction']
plt.xticks(ind, reaction_ids)
plt.ylabel('Reaction flux (mmol_per_gDW_per_hr)')
plt.title('rate through CPT1A/B/C-regulated reactions')

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('cpt1-regulated_reaction_flux.png',bbox_inches='tight')
plt.show()


# In[282]:


# create a dataset
y1 = (_59m_gene_exp_dict['ENSG00000205560']+_59m_gene_exp_dict['ENSG00000110090']+_59m_gene_exp_dict['ENSG00000169169'])
y2 = (heya8_gene_exp_dict['ENSG00000205560']+heya8_gene_exp_dict['ENSG00000110090']+heya8_gene_exp_dict['ENSG00000169169'])
y3 = (ov56_gene_exp_dict['ENSG00000205560']+ov56_gene_exp_dict['ENSG00000110090']+ov56_gene_exp_dict['ENSG00000169169'])
y4 = (caov3_gene_exp_dict['ENSG00000205560']+caov3_gene_exp_dict['ENSG00000110090']+caov3_gene_exp_dict['ENSG00000169169'])
y5 = (cov318_gene_exp_dict['ENSG00000205560']+cov318_gene_exp_dict['ENSG00000110090']+cov318_gene_exp_dict['ENSG00000169169'])
y6 = (oaw28_gene_exp_dict['ENSG00000205560']+oaw28_gene_exp_dict['ENSG00000110090']+oaw28_gene_exp_dict['ENSG00000169169'])

x = [y1,y2,y3,y4,y5,y6]
bars = ('59M','HEYA8','OV56','CAOV3','COV318','OAW28')
x_pos = np.arange(len(bars))

# Create bars with different colors
plt.bar(x_pos, x, color=[plotly.colors.qualitative.Dark24[3], plotly.colors.qualitative.G10[8]
                         , plotly.colors.qualitative.Light24[0], plotly.colors.qualitative.G10[7], 
                         plotly.colors.qualitative.Plotly[7],plotly.colors.qualitative.D3[2]])

plt.xticks(x_pos, bars)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('summed gene expression')
plt.title('summed gene expression of CPT1A/B/C')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Show graph
plt.savefig('summed_cpt1_exp.png',bbox_inches='tight')
plt.show()


# In[267]:


model.reactions.get_by_id('HMR_2591')


# In[268]:


model.reactions.get_by_id('ACRNtm') #this reaction transport butyryl-CoA from the cytosol to the mitochondria
#is this reaction exporting butyryl-CoA to a quicker extent in the HG samples?


# In[308]:


reaction = 'HMR_8676'
for n in range(len(comparison['reaction'])):
    if comparison.iloc[n,0] == reaction:
        print(comparison.iloc[n,:]) #flux in model is zero because we don't have information on the reaction so it has been switched off


# In[307]:


subset = comparison_5.loc[[111]]

y1 = subset.iloc[0,2]
y2 = subset.iloc[0,3]
y3 = subset.iloc[0,4]
y4 = subset.iloc[0,5]
y5 = subset.iloc[0,6]
y6 = subset.iloc[0,7]

x = [y1,y2,y3,y4,y5,y6]
bars = ('59M','HEYA8','OV56','CAOV3','COV318','OAW28')
x_pos = np.arange(len(bars))

# Create bars with different colors
plt.bar(x_pos, x, color=[plotly.colors.qualitative.Dark24[3], plotly.colors.qualitative.G10[8]
                         , plotly.colors.qualitative.Light24[0], plotly.colors.qualitative.G10[7], 
                         plotly.colors.qualitative.Plotly[7],plotly.colors.qualitative.D3[2]])

plt.xticks(x_pos, bars)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('Reaction flux (mmol_per_gDW_per_hr)')
plt.title('L-carnitine [c] + butyryl-CoA [c] ⇔ CoA [c] + O-butanoylcarnitine [c]')

# Show graph
plt.savefig('HMR_0159.png',bbox_inches='tight')
plt.show()


# In[309]:


comparison.loc[[6557]]


# In[311]:


subset = comparison.loc[[6557]]

y1 = subset.iloc[0,1]
y2 = subset.iloc[0,2]
y3 = subset.iloc[0,3]
y4 = subset.iloc[0,4]
y5 = subset.iloc[0,5]
y6 = subset.iloc[0,6]

x = [y1,y2,y3,y4,y5,y6]
bars = ('59M','HEYA8','OV56','CAOV3','COV318','OAW28')
x_pos = np.arange(len(bars))

# Create bars with different colors
plt.bar(x_pos, x, color=[plotly.colors.qualitative.Dark24[3], plotly.colors.qualitative.G10[8]
                         , plotly.colors.qualitative.Light24[0], plotly.colors.qualitative.G10[7], 
                         plotly.colors.qualitative.Plotly[7],plotly.colors.qualitative.D3[2]])

plt.xticks(x_pos, bars)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('Reaction flux (mmol_per_gDW_per_hr)')
plt.title('L-carnitine [c] ⇔ L-carnitine [s]')

# Show graph
plt.savefig('HMR_8676_flux.png',bbox_inches='tight')
plt.show()


# In[291]:


y1 = (_59m_gene_exp_dict['ENSG00000178537'])
y2 = (heya8_gene_exp_dict['ENSG00000178537'])
y3 = (ov56_gene_exp_dict['ENSG00000178537'])
y4 = (caov3_gene_exp_dict['ENSG00000178537'])
y5 = (cov318_gene_exp_dict['ENSG00000178537'])
y6 = (oaw28_gene_exp_dict['ENSG00000178537'])

x = [y1,y2,y3,y4,y5,y6]
bars = ('59M','HEYA8','OV56','CAOV3','COV318','OAW28')
x_pos = np.arange(len(bars))

# Create bars with different colors
plt.bar(x_pos, x, color=[plotly.colors.qualitative.Dark24[3], plotly.colors.qualitative.G10[8]
                         , plotly.colors.qualitative.Light24[0], plotly.colors.qualitative.G10[7], 
                         plotly.colors.qualitative.Plotly[7],plotly.colors.qualitative.D3[2]])

plt.xticks(x_pos, bars)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('CCLE; normalised gene expression')
plt.title('gene expression of SLC25A20')

# Show graph
plt.savefig('SLC25A20_exp.png',bbox_inches='tight')
plt.show()


# In[300]:


reaction = 'HMR_8676'
for n in range(len(comparison['reaction'])):
    if comparison.iloc[n,0] == reaction:
        print(comparison.iloc[n,:]) #increased export of L-carnitine from the cell would explain why HMR_0159 is lower in LG


# In[302]:


model.reactions.get_by_id('HMR_8676')


# In[312]:


reaction = 'HMR_0156'
for n in range(len(comparison['reaction'])):
    if comparison.iloc[n,0] == reaction:
        print(comparison.iloc[n,:])


# In[354]:


model.reactions.get_by_id('HMR_3107')


# In[523]:


reaction = 'HMR_2660'
for n in range(len(comparison['reaction'])):
    if comparison.iloc[n,0] == reaction:
        print(comparison.iloc[n,:])


# In[614]:


for n in range(len(comparison['reaction'])):
    if 'RE0577M' in comparison.iloc[n,0]:
        print(comparison.iloc[n,:])


# In[550]:


fa_synthesis_reactions = ['DESAT16_2','DESAT18_9','DESAT20_1','DESAT20_2','DESAT22_2p','FAS160COA']
fa_oxidation_reactions = ['FAOXC2242046x','FAOXC2252053m','FAOXC2442246x','FAOXC2452253x']


# In[514]:


for r in fa_synthesis_reactions:
    print(model.reactions.get_by_id(r).gene_reaction_rule)


# In[545]:


for n in range(len(comparison['reaction'])):
    if 'fatty' in model.reactions.get_by_id(comparison.iloc[n,0]).name:
        print(comparison.iloc[n,:])


# In[547]:


fa_pool_reactions = ['HMR_0668', 'HMR_0683', 'HMR_0753', 'HMR_8720', 'HMR_9209']
for r in fa_pool_reactions:
    print(r, model.reactions.get_by_id(r).name)


# In[549]:


for r in fa_pool_reactions:
    for n in range(len(comparison['reaction'])):
        if r in comparison.iloc[n,0]:
            print(comparison.iloc[n,:])


# In[551]:


for r in fa_synthesis_reactions:
    for n in range(len(comparison['reaction'])):
        if r in comparison.iloc[n,0]:
            print(comparison.iloc[n,:])


# In[601]:


for n in range(len(comparison['reaction'])):
    if 'HMR_2636' in comparison.iloc[n,0]:
        print(comparison.iloc[n,:])


# In[594]:


model.reactions.get_by_id('HMR_0155')


# In[595]:


print(_59m_gene_exp_dict['ENSG00000135218'])
print(heya8_gene_exp_dict['ENSG00000135218'])
print(ov56_gene_exp_dict['ENSG00000135218'])
print(caov3_gene_exp_dict['ENSG00000135218'])
print(cov318_gene_exp_dict['ENSG00000135218'])
print(oaw28_gene_exp_dict['ENSG00000135218'])


# In[105]:


for n in range(len(comparison['reaction'])):
    if 'HMR_4363' in comparison.iloc[n,0]:
        print(n)


# In[109]:


get_ipython().run_line_magic('store', 'subset')


# In[107]:


subset = comparison.loc[[16]]

y1 = subset['59m_absolute_flux']
y2 = subset['heya8_absolute_flux']
y3 = subset['ov56_absolute_flux']
y4 = subset['caov3_absolute_flux']
y5 = subset['cov318_absolute_flux']
y6 = subset['oaw28_absolute_flux']

N=len(subset)
ind = np.arange(N) 
width = 0.35      
plt.bar(ind, y1, width, label='59m_fluxes', color=plotly.colors.qualitative.Dark24[3])
plt.bar(ind + width, y2, width,
    label='heya8_fluxes', color=plotly.colors.qualitative.G10[8])
plt.bar(ind + width + width, y3, width,
    label='ov56_fluxes', color=plotly.colors.qualitative.Light24[0])
plt.bar(ind + 3*width, y4, width,
    label='caov3_fluxes', color=plotly.colors.qualitative.G10[7])
plt.bar(ind + 4*width, y5, width,
    label='cov318_fluxes', color=plotly.colors.qualitative.Plotly[7])
plt.bar(ind + 5*width, y6, width,
    label='oaw28_fluxes', color=plotly.colors.qualitative.D3[2])

reaction_ids = subset['reaction']
plt.xticks(ind, reaction_ids)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('Reaction flux (mmol_per_gDW_per_hr)')
plt.title('2-phospho-D-glycerate [c] ⇔ H2O [c] + PEP [c]')

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('MAR04363.png',bbox_inches='tight')
plt.show()


# In[591]:


#fa de novo synthesis
subset = comparison.loc[[8258,8266,8267,8365]]

y1 = subset['59m_absolute_flux']
y2 = subset['heya8_absolute_flux']
y3 = subset['ov56_absolute_flux']
y4 = subset['caov3_absolute_flux']
y5 = subset['cov318_absolute_flux']
y6 = subset['oaw28_absolute_flux']

N=len(subset)
ind = np.arange(N) 
width = 0.1      
plt.bar(ind, y1, width, label='59m_fluxes', color=plotly.colors.qualitative.Dark24[3])
plt.bar(ind + width, y2, width,
    label='heya8_fluxes', color=plotly.colors.qualitative.G10[8])
plt.bar(ind + width + width, y3, width,
    label='ov56_fluxes', color=plotly.colors.qualitative.Light24[0])
plt.bar(ind + 3*width, y4, width,
    label='caov3_fluxes', color=plotly.colors.qualitative.G10[7])
plt.bar(ind + 4*width, y5, width,
    label='cov318_fluxes', color=plotly.colors.qualitative.Plotly[7])
plt.bar(ind + 5*width, y6, width,
    label='oaw28_fluxes', color=plotly.colors.qualitative.D3[2])

reaction_ids = subset['reaction']
plt.xticks(ind, reaction_ids)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('Reaction flux (mmol_per_gDW_per_hr)')
plt.title('fatty acid de novo synthesis')

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('fa_synthesis_metabolism.png',bbox_inches='tight')
plt.show()


# In[596]:


y1 = (_59m_gene_exp_dict['ENSG00000135218'])
y2 = (heya8_gene_exp_dict['ENSG00000135218'])
y3 = (ov56_gene_exp_dict['ENSG00000135218'])
y4 = (caov3_gene_exp_dict['ENSG00000135218'])
y5 = (cov318_gene_exp_dict['ENSG00000135218'])
y6 = (oaw28_gene_exp_dict['ENSG00000135218'])

x = [y1,y2,y3,y4,y5,y6]
bars = ('59M','HEYA8','OV56','CAOV3','COV318','OAW28')
x_pos = np.arange(len(bars))

# Create bars with different colors
plt.bar(x_pos, x, color=[plotly.colors.qualitative.Dark24[3], plotly.colors.qualitative.G10[8]
                         , plotly.colors.qualitative.Light24[0], plotly.colors.qualitative.G10[7], 
                         plotly.colors.qualitative.Plotly[7],plotly.colors.qualitative.D3[2]])

plt.xticks(x_pos, bars)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('CCLE; normalised gene expression')
plt.title('gene expression of CD36')

# Show graph
plt.savefig('CD36_exp.png',bbox_inches='tight')
plt.show()


# In[600]:


subset = comparison.iloc[5171,:]

y1 = subset['59m_absolute_flux']*-1
y2 = subset['heya8_absolute_flux']*-1
y3 = subset['ov56_absolute_flux']*-1
y4 = subset['caov3_absolute_flux']*-1
y5 = subset['cov318_absolute_flux']*-1
y6 = subset['oaw28_absolute_flux']*-1

x = [y1,y2,y3,y4,y5,y6]
bars = ('59M','HEYA8','OV56','CAOV3','COV318','OAW28')
x_pos = np.arange(len(bars))

# Create bars with different colors
plt.bar(x_pos, x, color=[plotly.colors.qualitative.Dark24[3], plotly.colors.qualitative.G10[8]
                         , plotly.colors.qualitative.Light24[0], plotly.colors.qualitative.G10[7], 
                         plotly.colors.qualitative.Plotly[7],plotly.colors.qualitative.D3[2]])

plt.xticks(x_pos, bars)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('CCLE; normalised gene expression')
plt.title('butyrate [c] ⇔ butyrate [s]')

# Show graph
plt.savefig('HMR_0155_flux.png',bbox_inches='tight')
plt.show()


# In[610]:


model.reactions.get_by_id('FAS160COA').gene_reaction_rule


# In[612]:


y1=(_59m_gene_exp_dict['ENSG00000099194'])
y2=(heya8_gene_exp_dict['ENSG00000099194'])
y3=(ov56_gene_exp_dict['ENSG00000099194'])
y4=(caov3_gene_exp_dict['ENSG00000099194'])
y5=(cov318_gene_exp_dict['ENSG00000099194'])
y6=(oaw28_gene_exp_dict['ENSG00000099194'])
lg_scd_avg = (y1+y2+y3)/3
hg_scd_avg = (y4+y5+y6)/3
lg_scd_std = np.std([y1,y2,y3])
hg_scd_std = np.std([y3,y4,y5])
y1=(_59m_gene_exp_dict['ENSG00000149485'])
y2=(heya8_gene_exp_dict['ENSG00000149485'])
y3=(ov56_gene_exp_dict['ENSG00000149485'])
y4=(caov3_gene_exp_dict['ENSG00000149485'])
y5=(cov318_gene_exp_dict['ENSG00000149485'])
y6=(oaw28_gene_exp_dict['ENSG00000149485'])
lg_fads1_avg = (y1+y2+y3)/3
hg_fads1_avg = (y4+y5+y6)/3
lg_fads1_std = np.std([y1,y2,y3])
hg_fads1_std = np.std([y3,y4,y5])
y1=(_59m_gene_exp_dict['ENSG00000134824'])
y2=(heya8_gene_exp_dict['ENSG00000134824'])
y3=(ov56_gene_exp_dict['ENSG00000134824'])
y4=(caov3_gene_exp_dict['ENSG00000134824'])
y5=(cov318_gene_exp_dict['ENSG00000134824'])
y6=(oaw28_gene_exp_dict['ENSG00000134824'])
lg_fads2_avg = (y1+y2+y3)/3
hg_fads2_avg = (y4+y5+y6)/3
lg_fads2_std = np.std([y1,y2,y3])
hg_fads2_std = np.std([y3,y4,y5])


# In[613]:


print(lg_scd_avg, hg_scd_avg, lg_scd_std, hg_scd_std)
print(lg_fads1_avg, hg_fads1_avg, lg_fads1_std, hg_fads1_std)
print(lg_fads2_avg, hg_fads2_avg, lg_fads2_std, hg_fads2_std)


# In[617]:


#plot significant difference in lg/hg expression of SCD, FADS1 and FADS2 (involved in de novo FA synthesis)
groups = ['LG SCD', 'HG SCD', 'LG FADS1', 'HG FADS1', 'LG FADS2', 'HG FADS2']
x_pos = np.arange(len(groups))
means = [lg_scd_avg, hg_scd_avg, lg_fads1_avg, hg_fads1_avg, lg_fads2_avg, hg_fads2_avg]
error = [lg_scd_std, hg_scd_std, lg_fads1_std, hg_fads1_std, lg_fads2_std, hg_fads2_std]

fig, ax = plt.subplots()
ax.bar(x_pos, means, yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10, color=[plotly.colors.qualitative.Dark24[3],plotly.colors.qualitative.Light24[8], plotly.colors.qualitative.Dark24[3],plotly.colors.qualitative.Light24[8], plotly.colors.qualitative.Dark24[3],plotly.colors.qualitative.Light24[8]])
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
ax.set_ylabel('CCLE; normalised gene expression')
ax.set_xticks(x_pos)
ax.set_xticklabels(groups)
ax.set_title('SCD, FADS1 and FADS2 gene expression')
ax.yaxis.grid(True)

plt.tight_layout()
plt.savefig('fa_de_novo_gene_exp.png')
plt.show()


# In[619]:


print(lg_1_avg, hg_1_avg, lg_1_std, hg_1_std)
print(lg_2_avg, hg_2_avg, lg_2_std, hg_2_std)
print(lg_3_avg, hg_3_avg, lg_3_std, hg_3_std)


# In[621]:


subset = comparison.loc[[8258]]
subset.iloc[0,1]


# In[622]:


#fa de novo synthesis
subset = comparison.loc[[8258]]
y1 = subset.iloc[0,1]
y2 = subset.iloc[0,2]
y3 = subset.iloc[0,3]
y4 = subset.iloc[0,4]
y5 = subset.iloc[0,5]
y6 = subset.iloc[0,6]
lg_1_avg = (y1+y2+y3)/3
hg_1_avg = (y4+y5+y6)/3
lg_1_std = np.std([y1,y2,y3])
hg_1_std = np.std([y3,y4,y5])

subset = comparison.loc[[8266]]
y1 = subset.iloc[0,1]
y2 = subset.iloc[0,2]
y3 = subset.iloc[0,3]
y4 = subset.iloc[0,4]
y5 = subset.iloc[0,5]
y6 = subset.iloc[0,6]
lg_2_avg = (y1+y2+y3)/3
hg_2_avg = (y4+y5+y6)/3
lg_2_std = np.std([y1,y2,y3])
hg_2_std = np.std([y3,y4,y5])

subset = comparison.loc[[8267]]
y1 = subset.iloc[0,1]
y2 = subset.iloc[0,2]
y3 = subset.iloc[0,3]
y4 = subset.iloc[0,4]
y5 = subset.iloc[0,5]
y6 = subset.iloc[0,6]
lg_3_avg = (y1+y2+y3)/3
hg_3_avg = (y4+y5+y6)/3
lg_3_std = np.std([y1,y2,y3])
hg_3_std = np.std([y3,y4,y5])

groups = ['LG DESAT16_2', 'HG DESAT16_2', 'LG DESAT18_9', 'HG DESAT18_9', 'LG DESAT20_1', 'HG DESAT20_1']
x_pos = np.arange(len(groups))
means = [lg_1_avg, hg_1_avg, lg_2_avg, hg_2_avg, lg_3_avg, hg_3_avg]
error = [lg_1_std, hg_1_std, lg_2_std, hg_2_std, lg_3_std, hg_3_std]

fig, ax = plt.subplots()
ax.bar(x_pos, means, yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10, color=[plotly.colors.qualitative.Dark24[3],plotly.colors.qualitative.Light24[8], plotly.colors.qualitative.Dark24[3],plotly.colors.qualitative.Light24[8], plotly.colors.qualitative.Dark24[3],plotly.colors.qualitative.Light24[8]])
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
ax.set_ylabel('reaction rate mmol/gDW/hour')
ax.set_xticks(x_pos)
ax.set_xticklabels(groups)
ax.set_title('rate through fatty acid de novo synthesis reactions')
ax.yaxis.grid(True)

plt.tight_layout()
plt.savefig('fa_de_novo_rates.png')
plt.show()


# # central carbon metabolism

# In[857]:


fa_metabolism_reactions = ['C161CRN2t','HMR_2952','RE2912M','HMR_0226','HMR_0263','r0437','ALDD21','HMR_2194','HMR_3311']


# In[860]:


for r in fa_metabolism_reactions:
    for n in range(len(comparison_5['reaction'])):
        if r in comparison_5.iloc[n,1]:
            print(r)
            print(comparison_5.iloc[n,2:8])


# In[861]:


carnitine_metabolism = ['HMR_2592','ACRNtm','C4tmc','HMR_0161','HMR_2809','HMR_0159','HMR_2591','HMR_2634','HMR_2638','HMR_0158']


# In[863]:


for r in carnitine_metabolism:
    for n in range(len(comparison_5['reaction'])):
        if r in comparison_5.iloc[n,1]:
            print(r)
            print(comparison_5.iloc[n,2:8])


# In[864]:


lipid_metabolism = ['HMR_0668','HMR_0683','CDSm','HMR_0607','G3PD2m','HMR_0677','HMR_0588','HMR_0779','HMR_0748','HMR_0773','HMR_0746','HMR_0758','HMR_8219']


# In[865]:


for r in lipid_metabolism:
    for n in range(len(comparison_5['reaction'])):
        if r in comparison_5.iloc[n,1]:
            print(r)
            print(comparison_5.iloc[n,2:8])


# In[866]:


central_carbon = ['HMR_4358','HMR_6412','HMR_4394','HMR_4388','HMR_4281','HMR_4363','HMR_4355','HMR_6410','HMR_4408','HMR_4456','HMR_3958','HMR_4141','HMR_4354','HMR_4052']


# In[867]:


for r in central_carbon:
    for n in range(len(comparison_5['reaction'])):
        if r in comparison_5.iloc[n,1]:
            print(r)
            print(comparison_5.iloc[n,2:8])


# In[868]:


aa_metabolism = ['HMR_6725','HMR_3804','HMR_3802','HMR_3890','HMR_8611','HMR_3838','HMR_3835','HMR_3806','HMR_8436','HMR_8435','HMR_3784','HMR_4742','HMR_3839','HMR_3843','HMR_3841','HMR_3845','HMR_4667','HMR_4288','HMR_4597','HMR_4739','HMR_4685','HMR_4687','NADQNOXR']


# In[869]:


for r in aa_metabolism:
    for n in range(len(comparison_5['reaction'])):
        if r in comparison_5.iloc[n,1]:
            print(r)
            print(comparison_5.iloc[n,2:8])


# In[870]:


nucleotide_metabolism = ['HMR_4802','HMR_4799','HMR_4806','HMR_8461','HMR_4083','HMR_8477','HMR_7806','HMR_7810','HMR_6334','HMR_6332','HMR_5031','HMR_8476','HMR_5036','HMR_4843','HMR_6328','HMR_4611','r1432','r1431','HMR_4615','HMR_4171','HMR_4573','HMR_4802','HMR_4004','HMR_4814','HMR_4038','HMR_4042','HMR_4810','HMR_4406','HMR_4135','HMR_4808','HMR_4651','HMR_6612','HMR_4637','DTMPKm','HMR_3931','r1433']


# In[871]:


for r in nucleotide_metabolism:
    for n in range(len(comparison_5['reaction'])):
        if r in comparison_5.iloc[n,1]:
            print(r)
            print(comparison_5.iloc[n,2:8])


# In[872]:


fa_metabolism_reactions


# # individual plots

# In[874]:


#fa metabolism
for r in fa_metabolism_reactions:
    for n in range(len(comparison_5['reaction'])):
        if r in comparison_5.iloc[n,1]:
            print(comparison_5.iloc[n,:8])


# In[997]:


for n in range(len(comparison['reaction'])):
    if 'HMR_2193' in comparison.iloc[n,0]:
        print(comparison.iloc[n,:8])


# In[876]:


subset = comparison_5.loc[[279,107,309,104,105,294,277,109,130]]

y1 = subset['59m_absolute_flux']
y2 = subset['heya8_absolute_flux']
y3 = subset['ov56_absolute_flux']
y4 = subset['caov3_absolute_flux']
y5 = subset['cov318_absolute_flux']
y6 = subset['oaw28_absolute_flux']

N=len(subset)
ind = np.arange(N) 
width = 0.1      
plt.bar(ind, y1, width, label='59m_fluxes', color=plotly.colors.qualitative.Dark24[3])
plt.bar(ind + width, y2, width,
    label='heya8_fluxes', color=plotly.colors.qualitative.G10[8])
plt.bar(ind + width + width, y3, width,
    label='ov56_fluxes', color=plotly.colors.qualitative.Light24[0])
plt.bar(ind + 3*width, y4, width,
    label='caov3_fluxes', color=plotly.colors.qualitative.G10[7])
plt.bar(ind + 4*width, y5, width,
    label='cov318_fluxes', color=plotly.colors.qualitative.Plotly[7])
plt.bar(ind + 5*width, y6, width,
    label='oaw28_fluxes', color=plotly.colors.qualitative.D3[2])

reaction_ids = subset['reaction']
plt.xticks(ind, reaction_ids)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('Reaction flux (mmol_per_gDW_per_hr)')
plt.title('fatty acid metabolism')

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()


# In[888]:


model.reactions.get_by_id('HMR_3311')


# In[1401]:


get_ipython().run_line_magic('store', 'comparison_5')


# In[889]:


for r in carnitine_metabolism:
    for n in range(len(comparison_5['reaction'])):
        if r in comparison_5.iloc[n,1]:
            print(comparison_5.iloc[n,:8])


# In[890]:


subset = comparison_5.loc[[118,276,311,115,124,111,112,121,122,205]]

y1 = subset['59m_absolute_flux']
y2 = subset['heya8_absolute_flux']
y3 = subset['ov56_absolute_flux']
y4 = subset['caov3_absolute_flux']
y5 = subset['cov318_absolute_flux']
y6 = subset['oaw28_absolute_flux']

N=len(subset)
ind = np.arange(N) 
width = 0.1      
plt.bar(ind, y1, width, label='59m_fluxes', color=plotly.colors.qualitative.Dark24[3])
plt.bar(ind + width, y2, width,
    label='heya8_fluxes', color=plotly.colors.qualitative.G10[8])
plt.bar(ind + width + width, y3, width,
    label='ov56_fluxes', color=plotly.colors.qualitative.Light24[0])
plt.bar(ind + 3*width, y4, width,
    label='caov3_fluxes', color=plotly.colors.qualitative.G10[7])
plt.bar(ind + 4*width, y5, width,
    label='cov318_fluxes', color=plotly.colors.qualitative.Plotly[7])
plt.bar(ind + 5*width, y6, width,
    label='oaw28_fluxes', color=plotly.colors.qualitative.D3[2])

reaction_ids = subset['reaction']
plt.xticks(ind, reaction_ids)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('Reaction flux (mmol_per_gDW_per_hr)')
plt.title('carnitine metabolism')

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()


# In[1096]:


for n in range(len(comparison['reaction'])):
    if 'HMR_7160' in comparison.iloc[n,0]:
        print(comparison.iloc[n,:8])


# In[5]:


get_ipython().run_line_magic('store', '-r comparison')


# In[53]:


for n in range(len(comparison['reaction'])):
    if 'HMR_7160' in comparison.iloc[n,0]:
        print(comparison.iloc[n,:8])


# In[1077]:


(13.179972+11.228799+12.096528)/(12.40578+12.287718+11.604815)


# In[894]:


model.reactions.get_by_id('HMR_0161')


# In[895]:


for r in lipid_metabolism:
    for n in range(len(comparison_5['reaction'])):
        if r in comparison_5.iloc[n,1]:
            print(comparison_5.iloc[n,:8])


# In[907]:


subset = comparison_5.loc[[126,127,281,149,286,145,144]]

y1 = subset['59m_absolute_flux']
y2 = subset['heya8_absolute_flux']
y3 = subset['ov56_absolute_flux']
y4 = subset['caov3_absolute_flux']
y5 = subset['cov318_absolute_flux']
y6 = subset['oaw28_absolute_flux']

N=len(subset)
ind = np.arange(N) 
width = 0.1      
plt.bar(ind, y1, width, label='59m_fluxes', color=plotly.colors.qualitative.Dark24[3])
plt.bar(ind + width, y2, width,
    label='heya8_fluxes', color=plotly.colors.qualitative.G10[8])
plt.bar(ind + width + width, y3, width,
    label='ov56_fluxes', color=plotly.colors.qualitative.Light24[0])
plt.bar(ind + 3*width, y4, width,
    label='caov3_fluxes', color=plotly.colors.qualitative.G10[7])
plt.bar(ind + 4*width, y5, width,
    label='cov318_fluxes', color=plotly.colors.qualitative.Plotly[7])
plt.bar(ind + 5*width, y6, width,
    label='oaw28_fluxes', color=plotly.colors.qualitative.D3[2])

reaction_ids = subset['reaction']
plt.xticks(ind, reaction_ids)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('Reaction flux (mmol_per_gDW_per_hr)')
plt.title('lipid metabolism')

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()


# In[906]:


subset = comparison_5.loc[[139,136,138,135,137,141]]

y1 = subset['59m_absolute_flux']
y2 = subset['heya8_absolute_flux']
y3 = subset['ov56_absolute_flux']
y4 = subset['caov3_absolute_flux']
y5 = subset['cov318_absolute_flux']
y6 = subset['oaw28_absolute_flux']

N=len(subset)
ind = np.arange(N) 
width = 0.1      
plt.bar(ind, y1, width, label='59m_fluxes', color=plotly.colors.qualitative.Dark24[3])
plt.bar(ind + width, y2, width,
    label='heya8_fluxes', color=plotly.colors.qualitative.G10[8])
plt.bar(ind + width + width, y3, width,
    label='ov56_fluxes', color=plotly.colors.qualitative.Light24[0])
plt.bar(ind + 3*width, y4, width,
    label='caov3_fluxes', color=plotly.colors.qualitative.G10[7])
plt.bar(ind + 4*width, y5, width,
    label='cov318_fluxes', color=plotly.colors.qualitative.Plotly[7])
plt.bar(ind + 5*width, y6, width,
    label='oaw28_fluxes', color=plotly.colors.qualitative.D3[2])

reaction_ids = subset['reaction']
plt.xticks(ind, reaction_ids)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('Reaction flux (mmol_per_gDW_per_hr)')
plt.title('sphingolipid metabolism')

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()


# In[905]:


model.reactions.get_by_id('HMR_8219')


# In[908]:


for r in central_carbon:
    for n in range(len(comparison_5['reaction'])):
        if r in comparison_5.iloc[n,1]:
            print(comparison_5.iloc[n,:8])


# In[912]:


subset = comparison_5.loc[[1,0,2,93,92,15,14]]

y1 = subset['59m_absolute_flux']
y2 = subset['heya8_absolute_flux']
y3 = subset['ov56_absolute_flux']
y4 = subset['caov3_absolute_flux']
y5 = subset['cov318_absolute_flux']
y6 = subset['oaw28_absolute_flux']

N=len(subset)
ind = np.arange(N) 
width = 0.1      
plt.bar(ind, y1, width, label='59m_fluxes', color=plotly.colors.qualitative.Dark24[3])
plt.bar(ind + width, y2, width,
    label='heya8_fluxes', color=plotly.colors.qualitative.G10[8])
plt.bar(ind + width + width, y3, width,
    label='ov56_fluxes', color=plotly.colors.qualitative.Light24[0])
plt.bar(ind + 3*width, y4, width,
    label='caov3_fluxes', color=plotly.colors.qualitative.G10[7])
plt.bar(ind + 4*width, y5, width,
    label='cov318_fluxes', color=plotly.colors.qualitative.Plotly[7])
plt.bar(ind + 5*width, y6, width,
    label='oaw28_fluxes', color=plotly.colors.qualitative.D3[2])

reaction_ids = subset['reaction']
plt.xticks(ind, reaction_ids)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('Reaction flux (mmol_per_gDW_per_hr)')
plt.title('central carbon metabolism 2')

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()


# In[910]:


subset = comparison_5.loc[[3,7,5,4,6,96,91]]

y1 = subset['59m_absolute_flux']
y2 = subset['heya8_absolute_flux']
y3 = subset['ov56_absolute_flux']
y4 = subset['caov3_absolute_flux']
y5 = subset['cov318_absolute_flux']
y6 = subset['oaw28_absolute_flux']

N=len(subset)
ind = np.arange(N) 
width = 0.1      
plt.bar(ind, y1, width, label='59m_fluxes', color=plotly.colors.qualitative.Dark24[3])
plt.bar(ind + width, y2, width,
    label='heya8_fluxes', color=plotly.colors.qualitative.G10[8])
plt.bar(ind + width + width, y3, width,
    label='ov56_fluxes', color=plotly.colors.qualitative.Light24[0])
plt.bar(ind + 3*width, y4, width,
    label='caov3_fluxes', color=plotly.colors.qualitative.G10[7])
plt.bar(ind + 4*width, y5, width,
    label='cov318_fluxes', color=plotly.colors.qualitative.Plotly[7])
plt.bar(ind + 5*width, y6, width,
    label='oaw28_fluxes', color=plotly.colors.qualitative.D3[2])

reaction_ids = subset['reaction']
plt.xticks(ind, reaction_ids)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('Reaction flux (mmol_per_gDW_per_hr)')
plt.title('central carbon metabolism 1')

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()


# In[ ]:


#HMR_4406 is in the nucleotide reactions


# In[938]:


nucleotide_metabolism


# In[989]:


for r in central_carbon:
    for n in range(len(comparison_5['reaction'])):
        if r in comparison_5.iloc[n,1]:
            print(comparison_5.iloc[n,:8])


# In[994]:


for n in range(len(comparison['reaction'])):
    if 'HMR_4381' in comparison.iloc[n,0]:
        print(comparison.iloc[n,:8])


# In[983]:


(12.76+12.29+11.6)/(13.18+11.23+12.1)


# In[934]:


#upregulated in LG
LG_1 = [32,31,52,35,17,18,34,21,33] #1,2.8,2.25 pattern
LG_2 = [224,298,29,19,30] #4,5.9,4.3 pattern
LG_3 = [234,26,297,45,299] #>10
LG_4 = [55,57,242,243,236,235,260,213,16,47,313]#all others
HG_1 = [51,223,20,25,48]


# In[1346]:


for r in aa_metabolism:
    if 'ENSG00000091140' in model.reactions.get_by_id(r).gene_reaction_rule:
        print(r) 


# In[1399]:


for n in range(len(comparison['reaction'])):
    if 'HMR_4127' in comparison.iloc[n,0]:
        print(comparison.iloc[n,:8])


# In[1249]:


print((0.000246+0.000508+0.000036)/3)
print((0.000245+0.000228+0.000249)/3)


# In[935]:


subset = comparison_5.loc[LG_1]

y1 = subset['59m_absolute_flux']
y2 = subset['heya8_absolute_flux']
y3 = subset['ov56_absolute_flux']
y4 = subset['caov3_absolute_flux']
y5 = subset['cov318_absolute_flux']
y6 = subset['oaw28_absolute_flux']

N=len(subset)
ind = np.arange(N) 
width = 0.1      
plt.bar(ind, y1, width, label='59m_fluxes', color=plotly.colors.qualitative.Dark24[3])
plt.bar(ind + width, y2, width,
    label='heya8_fluxes', color=plotly.colors.qualitative.G10[8])
plt.bar(ind + width + width, y3, width,
    label='ov56_fluxes', color=plotly.colors.qualitative.Light24[0])
plt.bar(ind + 3*width, y4, width,
    label='caov3_fluxes', color=plotly.colors.qualitative.G10[7])
plt.bar(ind + 4*width, y5, width,
    label='cov318_fluxes', color=plotly.colors.qualitative.Plotly[7])
plt.bar(ind + 5*width, y6, width,
    label='oaw28_fluxes', color=plotly.colors.qualitative.D3[2])

reaction_ids = subset['reaction']
plt.xticks(ind, reaction_ids)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('Reaction flux (mmol_per_gDW_per_hr)')
plt.title('nucleotide metabolism 1')

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()


# In[930]:


subset = comparison_5.loc[LG_2]

y1 = subset['59m_absolute_flux']
y2 = subset['heya8_absolute_flux']
y3 = subset['ov56_absolute_flux']
y4 = subset['caov3_absolute_flux']
y5 = subset['cov318_absolute_flux']
y6 = subset['oaw28_absolute_flux']

N=len(subset)
ind = np.arange(N) 
width = 0.1      
plt.bar(ind, y1, width, label='59m_fluxes', color=plotly.colors.qualitative.Dark24[3])
plt.bar(ind + width, y2, width,
    label='heya8_fluxes', color=plotly.colors.qualitative.G10[8])
plt.bar(ind + width + width, y3, width,
    label='ov56_fluxes', color=plotly.colors.qualitative.Light24[0])
plt.bar(ind + 3*width, y4, width,
    label='caov3_fluxes', color=plotly.colors.qualitative.G10[7])
plt.bar(ind + 4*width, y5, width,
    label='cov318_fluxes', color=plotly.colors.qualitative.Plotly[7])
plt.bar(ind + 5*width, y6, width,
    label='oaw28_fluxes', color=plotly.colors.qualitative.D3[2])

reaction_ids = subset['reaction']
plt.xticks(ind, reaction_ids)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('Reaction flux (mmol_per_gDW_per_hr)')
plt.title('nucleotide metabolism 2')

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()


# In[931]:


subset = comparison_5.loc[LG_3]

y1 = subset['59m_absolute_flux']
y2 = subset['heya8_absolute_flux']
y3 = subset['ov56_absolute_flux']
y4 = subset['caov3_absolute_flux']
y5 = subset['cov318_absolute_flux']
y6 = subset['oaw28_absolute_flux']

N=len(subset)
ind = np.arange(N) 
width = 0.1      
plt.bar(ind, y1, width, label='59m_fluxes', color=plotly.colors.qualitative.Dark24[3])
plt.bar(ind + width, y2, width,
    label='heya8_fluxes', color=plotly.colors.qualitative.G10[8])
plt.bar(ind + width + width, y3, width,
    label='ov56_fluxes', color=plotly.colors.qualitative.Light24[0])
plt.bar(ind + 3*width, y4, width,
    label='caov3_fluxes', color=plotly.colors.qualitative.G10[7])
plt.bar(ind + 4*width, y5, width,
    label='cov318_fluxes', color=plotly.colors.qualitative.Plotly[7])
plt.bar(ind + 5*width, y6, width,
    label='oaw28_fluxes', color=plotly.colors.qualitative.D3[2])

reaction_ids = subset['reaction']
plt.xticks(ind, reaction_ids)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('Reaction flux (mmol_per_gDW_per_hr)')
plt.title('nucleotide metabolism 3')

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()


# In[932]:


subset = comparison_5.loc[LG_4]

y1 = subset['59m_absolute_flux']
y2 = subset['heya8_absolute_flux']
y3 = subset['ov56_absolute_flux']
y4 = subset['caov3_absolute_flux']
y5 = subset['cov318_absolute_flux']
y6 = subset['oaw28_absolute_flux']

N=len(subset)
ind = np.arange(N) 
width = 0.1      
plt.bar(ind, y1, width, label='59m_fluxes', color=plotly.colors.qualitative.Dark24[3])
plt.bar(ind + width, y2, width,
    label='heya8_fluxes', color=plotly.colors.qualitative.G10[8])
plt.bar(ind + width + width, y3, width,
    label='ov56_fluxes', color=plotly.colors.qualitative.Light24[0])
plt.bar(ind + 3*width, y4, width,
    label='caov3_fluxes', color=plotly.colors.qualitative.G10[7])
plt.bar(ind + 4*width, y5, width,
    label='cov318_fluxes', color=plotly.colors.qualitative.Plotly[7])
plt.bar(ind + 5*width, y6, width,
    label='oaw28_fluxes', color=plotly.colors.qualitative.D3[2])

reaction_ids = subset['reaction']
plt.xticks(ind, reaction_ids)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('Reaction flux (mmol_per_gDW_per_hr)')
plt.title('nucleotide metabolism 4')

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()


# In[933]:


subset = comparison_5.loc[HG_1]

y1 = subset['59m_absolute_flux']
y2 = subset['heya8_absolute_flux']
y3 = subset['ov56_absolute_flux']
y4 = subset['caov3_absolute_flux']
y5 = subset['cov318_absolute_flux']
y6 = subset['oaw28_absolute_flux']

N=len(subset)
ind = np.arange(N) 
width = 0.1      
plt.bar(ind, y1, width, label='59m_fluxes', color=plotly.colors.qualitative.Dark24[3])
plt.bar(ind + width, y2, width,
    label='heya8_fluxes', color=plotly.colors.qualitative.G10[8])
plt.bar(ind + width + width, y3, width,
    label='ov56_fluxes', color=plotly.colors.qualitative.Light24[0])
plt.bar(ind + 3*width, y4, width,
    label='caov3_fluxes', color=plotly.colors.qualitative.G10[7])
plt.bar(ind + 4*width, y5, width,
    label='cov318_fluxes', color=plotly.colors.qualitative.Plotly[7])
plt.bar(ind + 5*width, y6, width,
    label='oaw28_fluxes', color=plotly.colors.qualitative.D3[2])

reaction_ids = subset['reaction']
plt.xticks(ind, reaction_ids)
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('Reaction flux (mmol_per_gDW_per_hr)')
plt.title('nucleotide metabolism 5')

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()


# In[1027]:


for r in aa_metabolism:
    for n in range(len(comparison_5['reaction'])):
        if r in comparison_5.iloc[n,1]:
            print(comparison_5.iloc[n,:8])


# In[1400]:


len(aa_metabolism)


# In[366]:


model.reactions.get_by_id('HMR_6296')


# In[572]:


for n in range(len(comparison['reaction'])):
    if 'HMR_2654' in comparison.iloc[n,0]:
        print(comparison.iloc[n,:8])


# In[566]:


hg_avg_ol = []
for n in range(len(comparison['reaction'])):
    if 'HMR_2663' in comparison.iloc[n,0]:
        print(comparison.iloc[n,:8])
        hg_avg_ol.append(comparison.iloc[n,4])
        hg_avg_ol.append(comparison.iloc[n,5])
        hg_avg_ol.append(comparison.iloc[n,6])


# In[565]:


hg_avg_ol_2 = []
for n in range(len(comparison['reaction'])):
    if 'HMR_2827' in comparison.iloc[n,0]:
        print(comparison.iloc[n,:8])
        hg_avg_ol_2.append(comparison.iloc[n,4])
        hg_avg_ol_2.append(comparison.iloc[n,5])
        hg_avg_ol_2.append(comparison.iloc[n,6])


# In[568]:


(sum(hg_avg_ol_2))/3


# In[563]:


(sum(hg_avg_ol)*-1)/3


# In[518]:


glutamine_agrees = ['HMR_4034','HMR_4300','HMR_4406','HMR_4808','HMR_5136']
glutamine_disagrees = ['HMR_3890','HMR_3903']
for n in range(len(comparison['reaction'])):
    for r in glutamine_agrees:
        if r in comparison.iloc[n,0]:
            print(comparison.iloc[n,:8])


# In[519]:


for n in range(len(comparison['reaction'])):
    for r in glutamine_disagrees:
        if r in comparison.iloc[n,0]:
            print(comparison.iloc[n,:8])


# In[520]:


glutamine_produced = ['HMR_3890']
glutamine_consumed = ['HMR_4034','HMR_4300','HMR_4406','HMR_4808','HMR_5136','HMR_3903']


# In[527]:


_59m_produced = []
heya8_produced = []
ov56_produced = []
caov3_produced = []
cov318_produced = []
oaw28_produced = []
for n in range(len(comparison['reaction'])):
    for r in glutamine_produced:
        if r in comparison.iloc[n,0]:
            _59m_produced.append(comparison.iloc[n,1])
            heya8_produced.append(comparison.iloc[n,2])
            ov56_produced.append(comparison.iloc[n,3])
            caov3_produced.append(comparison.iloc[n,4])
            cov318_produced.append(comparison.iloc[n,5])
            oaw28_produced.append(comparison.iloc[n,6])


# In[546]:


print(sum(_59m_produced))
print(sum(heya8_produced))
print(sum(ov56_produced))
print(sum(caov3_produced))
print(sum(cov318_produced))
print(sum(oaw28_produced))
print('\n')
print(sum(_59m_consumed))
print(sum(heya8_consumed))
print(sum(ov56_consumed))
print(sum(caov3_consumed))
print(sum(cov318_consumed))
print(sum(oaw28_consumed))


# In[548]:


(sum(_59m_consumed)+sum(heya8_consumed)+sum(ov56_consumed))/3


# In[549]:


(sum(caov3_consumed)+sum(cov318_consumed)+sum(oaw28_consumed))/3


# In[529]:


_59m_consumed = []
heya8_consumed = []
ov56_consumed = []
caov3_consumed = []
cov318_consumed = []
oaw28_consumed = []
for n in range(len(comparison['reaction'])):
    for r in glutamine_consumed:
        if r in comparison.iloc[n,0]:
            _59m_consumed.append(comparison.iloc[n,1])
            heya8_consumed.append(comparison.iloc[n,2])
            ov56_consumed.append(comparison.iloc[n,3])
            caov3_consumed.append(comparison.iloc[n,4])
            cov318_consumed.append(comparison.iloc[n,5])
            oaw28_consumed.append(comparison.iloc[n,6])


# In[534]:


_59m_total = _59m_produced[0] - sum(_59m_consumed)
heya8_total = heya8_produced[0] - sum(heya8_consumed)
ov56_total = ov56_produced[0] - sum(ov56_consumed)
caov3_total = caov3_produced[0] - sum(caov3_consumed)
cov318_total = cov318_produced[0] - sum(cov318_consumed)
oaw28_total = oaw28_produced[0] - sum(oaw28_consumed)


# In[538]:


print(_59m_total, heya8_total, ov56_total, caov3_total, cov318_total, oaw28_total)


# # correlate metabolomics with flux predictions

# In[272]:


x = 'glutamine'
rxn_list = []
for n in range(len(comparison['reaction'])):
    if x in model.reactions.get_by_id(comparison.iloc[n,0]).name:
        if 'c' in model.reactions.get_by_id(comparison.iloc[n,0]).compartments:
            if (comparison.iloc[n,1:]).to_list() != [0,0,0,0,0,0]:
                rxn_list.append(comparison.iloc[n,0])
                print(model.reactions.get_by_id(comparison.iloc[n,0]).name, comparison.iloc[n,:])
                print('\n')


# In[273]:


rxn_list


# In[573]:


ex = []
for e in model.exchanges():
    ex.append(e)
print(len(e))


# In[ ]:




