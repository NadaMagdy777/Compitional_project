#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pyopenms import *
enter = []
seq =""
f=FASTAFile()

f.load("file.fasta",enter)
for i in enter:
    seq=seq+i.sequence
dig = ProteaseDigestion()
dig.getEnzymeName() 
bsa = AASequence.fromString(seq)
result = []
dig.digest(bsa, result)
for i in result:
   print(i.toString())


# In[2]:


from pyopenms import *
exp = MSExperiment()
MzMLFile().load("Fusion_180220_40.mzML", exp)
spectrum = exp.getSpectra()


# In[3]:


from pyopenms import *
enter = []
seq =""
f=FASTAFile()
f.load("file.fasta",enter)
for i in enter :
    seq=seq+i.sequence
dig = ProteaseDigestion()
dig.getEnzymeName()
bsa = AASequence.fromString(seq)
result = []
dig.digest(bsa, result)
peptides = [AASequence.fromString(i.toString()) for i in result]

for peptide in peptides:
    tsg = TheoreticalSpectrumGenerator()
    theo_spectrum = MSSpectrum()
    spec1 = MSSpectrum()

    p = Param()
    p.setValue("add_b_ions", "true")
    p.setValue("add_y_ions", "true")
    p.setValue("add_losses", "true")
    p.setValue("add_metainfo", "true")
    tsg.setParameters(p)
    tsg.getSpectrum(spec1, peptide, 1, 1) 
    print("Spectrum 1 of", peptide, "has", spec1.size(), "peaks.")
    for ion, peak in zip(spec1.getStringDataArrays()[0], spec1):
        print(ion.decode(), "is generated at m/z", peak.getMZ())


# In[4]:


for i in range(0,1):
    alignment =[]
    matchpeaks =[]
    observed_spectrum = spectrum[i]
    tsg = TheoreticalSpectrumGenerator()
    spa = SpectrumAlignment()
    p = spa.getParameters()
    p.setValue("tolerance", 0.5)
    p.setValue("is_relative_tolerance", "false")
    spa.setParameters(p)
    for j in range(len(result)):
        theo_spectrum = MSSpectrum()
        p = spa.getParameters()
        p.setValue("add_y_ions", "true")
        p.setValue("add_b_ions", "true")
        p.setValue("add_metainfo", "true")
        tsg.setParameters(p)
        tsg.getSpectrum(theo_spectrum, result[j], 1, 2)
        alignment =[]
        spa.getSpectrumAlignment(alignment, theo_spectrum, observed_spectrum)
        matchpeaks.append(len(alignment))
    print("Number of matched peaks: " + str(max(matchpeaks)))
    print("ion\ttheo. m/z\tobserved m/z")
    theo_spectrum = MSSpectrum()
    tsg.getSpectrum(theo_spectrum, result[matchpeaks.index(max(matchpeaks))], 1 , 2)
    spa.getSpectrumAlignment(alignment, theo_spectrum, observed_spectrum)
    for theo_idx, obs_idx in alignment:
       ion_name = theo_spectrum.getStringDataArrays()[0][theo_idx].decode()
       ion_charge = theo_spectrum.getIntegerDataArrays()[0][theo_idx]
       print(ion_name + "\t" + str(ion_charge) + "\t"
          + str(theo_spectrum[theo_idx].getMZ())
          + "\t" + str(observed_spectrum[obs_idx].getMZ()))
 
    theo_mz, theo_int, obs_mz, obs_int = [], [], [], []
    for theo_idx, obs_idx in alignment:
       theo_mz.append(theo_spectrum[theo_idx].getMZ())
       theo_int.append(theo_spectrum[theo_idx].getIntensity())
       obs_mz.append(observed_spectrum[obs_idx].getMZ())
       obs_int.append(observed_spectrum[obs_idx].getIntensity())
  
        


# In[11]:


print(str(result[matchpeaks.index(max(matchpeaks))]))


# In[6]:


print(len(spectrum)) 


# In[7]:


import numpy as np
from matplotlib import pyplot as plt
for r in result:
    peptide = AASequence.fromString(r.toString())
    tsg.getSpectrum(spec1,r,1,1)
plt.bar(spec1.get_peaks()[0], spec1.get_peaks()[1], snap=False) 
plt.xlabel("m/z")
plt.ylabel("intensity")
plt.show()


# In[8]:


import numpy as np
from matplotlib import pyplot as plt
def mirror_plot(obs_mz, obs_int, theo_mz, theo_int, title):
    obs_int = [element / max(obs_int) for element in obs_int] 
    theo_int = [element * -1 for element in theo_int] 
    plt.figure(figsize=(12,8))
    plt.bar(obs_mz, obs_int, width = 3.0)
    plt.bar(theo_mz, theo_int, width = 3.0)
    plt.title(title)
    plt.ylabel('intensity')
    plt.xlabel('m/z')

obs_mz, obs_int = observed_spectrum.get_peaks()
print(min(obs_mz)) 
print(max(obs_mz)) 


theo_mz, theo_int = [], []
for mz, intensity in zip(*theo_spectrum.get_peaks()):
    if mz >= min(obs_mz) and mz <= max(obs_mz):
        theo_mz.append(mz)
        theo_int.append(intensity)

title = 'Observed vs theoretical spectrum'
mirror_plot(obs_mz, obs_int, theo_mz, theo_int, title)


# In[10]:


for i in range(0,1):
    alignment =[]
    matchpeaks =[]
    observed_spectrum = spectrum[i]
    tsg = TheoreticalSpectrumGenerator()
    spa = SpectrumAlignment()
    p = spa.getParameters()
    p.setValue("tolerance", 0.5)
    p.setValue("is_relative_tolerance", "false")
    spa.setParameters(p)
    for j in range(len(result)):
        theo_spectrum = MSSpectrum()
        p = spa.getParameters()
        p.setValue("add_y_ions", "true")
        p.setValue("add_b_ions", "true")
        p.setValue("add_metainfo", "true")
        tsg.setParameters(p)
        tsg.getSpectrum(theo_spectrum, result[j], 1, 2)
        alignment =[]
        spa.getSpectrumAlignment(alignment, theo_spectrum, observed_spectrum)
        matchpeaks.append(len(alignment))
    print("Number of matched peaks: " + str(max(matchpeaks)))
    print("ion\ttheo. m/z\tobserved m/z")
    theo_spectrum = MSSpectrum()
    tsg.getSpectrum(theo_spectrum, result[matchpeaks.index(max(matchpeaks))], 1 , 2)
    spa.getSpectrumAlignment(alignment, theo_spectrum, observed_spectrum)
    for theo_idx, obs_idx in alignment:
       ion_name = theo_spectrum.getStringDataArrays()[0][theo_idx].decode()
       ion_charge = theo_spectrum.getIntegerDataArrays()[0][theo_idx]
       print(ion_name + "\t" + str(ion_charge) + "\t"
          + str(theo_spectrum[theo_idx].getMZ())
          + "\t" + str(observed_spectrum[obs_idx].getMZ()))
 
    theo_mz, theo_int, obs_mz, obs_int = [], [], [], []
    for theo_idx, obs_idx in alignment:
       theo_mz.append(theo_spectrum[theo_idx].getMZ())
       theo_int.append(theo_spectrum[theo_idx].getIntensity())
       obs_mz.append(observed_spectrum[obs_idx].getMZ())
       obs_int.append(observed_spectrum[obs_idx].getIntensity())
    title = 'Observed vs theoretical spectrum(aligned)'
    mirror_plot(obs_mz, obs_int, theo_mz, theo_int, title)


# In[ ]:




