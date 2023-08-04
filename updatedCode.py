###################################                            
##Author: Mehdi Manjoura AKA 98mm -- Supervisor: Erik Van Heesbeen -- Lab: Philippe Campeau
###################################

###################################
###################################
import statistics                ##
import numpy as np               ##
import matplotlib.pyplot as plt  ##
import matplotlib.patches as mpatches
import pandas as pd              ##
import seaborn as sns            ##
import glob                      ##
import re                        ##
import os                        
from pathlib import Path         ##
from scipy import stats          ##
import itertools                 ##    
import math                      ##
###################################
###################################


bin = 12
fi = 'raw.csv'
cut = 100

quant = 'newsmth.csv'
annot = 'data.csv'
###################################
###################################

#global variables

listt =[]
LFC = []
pVal = []
bMean = []
LFCBM = []
Mlis =[]
mLis = []
sil =[]
listFreqq = []
lis = []
#bin = int(input("How many bins are you trying to work with:"))
#print("The annotated file, is the output of chipseeker, later scripts from chipseeker and Seacr will be added to this code!")
#file = input("Please input your annotated peaks file, in the format CHDx_xIGGx.csv:")
print("Warning::::::: Our Deseq2 file is named raw.csv in this directory!")
################################
################################


def cutoff(fi,cut):

    print("Hello to this DEG tool, please make sure you have a file under the following format : filename.csv!")
    print("*"*100)
    #Taking the input file (the ouptut of Deseq2/EdgeR)   
    #fi = input("Can you input your filename.csv:")
    df = pd.read_csv(fi, sep = ',')
    pd.set_option("display.width", None)
    pd.set_option("display.max_rows", None)
    #taking the cutoff value
    #cutoff = int(input("What is your basemean cutoff value:"))
    #dropping all values less than our cutoff value
    
    df.drop(df[df['baseMean'] < cut].index, inplace =True)
    #sorting our basemean
    df = df.sort_values(['baseMean'], ascending = False)
    df.reset_index(drop = True, inplace = True)
    #Cleaning data and saving only genes with a basemean > 0 
    df = df[df["GeneId"].notna()]
    #saving to a file that will delete later,just keeping it for V purposes
    df.to_csv(r"./Module1Files/output.csv", index = False)
    


def log2FC():
    
    df1 = pd.read_csv("./Module1Files/output.csv", sep = ",")
    df1.sort_values(['log2FoldChange'], ascending= False, inplace = True)
    df1.to_csv(r'./Module1Files/outputLFC.csv', index = False)
    


def pvalueSort():
    df2 = pd.read_csv("./Module1Files/outputLFC.csv", sep = ",")
    #taking half of the length of the file
    a = int(len(df2)/2)
    #splitting our file into 2 chunks
    df2 = df2.truncate(before = 0, after = a)
    #sorting pvalues from low to high 
    df2 = df2.sort_values("pvalue", ascending = True)
    
    #resetting index, so that we dont get extra columns
    df2.reset_index(inplace = True, drop = True)
    #saving the first chunk to a file
    df2.to_csv(r"./Module1Files/Trunc1.csv", index = False)
    df3 = pd.read_csv("./Module1Files/outputLFC.csv", sep = ",")
    #taking the second chunk of the file
    df3 = df3.truncate(before = a, after = len(df3))
    #sort pvalue from high to low
    df3 = df3.sort_values("pvalue", ascending = False)
    
    #saving to a file
    df3.to_csv(r"./Module1Files/Trunc2.csv", index = False)
    #combining the 2 sorted chunks
    comb = df2.append(df3)
    comb.to_csv(r"./Module1Files/finalFile.csv", index = False)
    


def binbar(bin):
    df4 = pd.read_csv("./Module1Files/finalFile.csv", sep = ",") 
    #fill the values with 0 pvalues or LFC with 0s
    df4.drop(df4.columns.difference(['GeneId', 'baseMean', 'log2FoldChange','pvalue']), axis = 1, inplace = True)
    df4 = df4.dropna()
    df4 = df4.reset_index(drop = True)
    df4.to_csv(r"./Module1Files/finalModified.csv",index = False)
    interval = int(len(df4)/ bin)
    #if there is a remainder from the division, we will append it to the last bi
    #saving the frequency of peaks per genes into 20 csv files, they also represent the genes in our bins respectively 
    for i, binn in enumerate(pd.read_csv("./Module1Files/finalFile.csv", chunksize = interval)):
        binn = binn["GeneId"].value_counts()
        binn = binn.reset_index()
        binn.to_csv(r"./Module1Files/GeneBinFreq{}.csv".format(i), index = False)
   
    print("*"*100)
    start = 0
    print("*"*100)
    stop = len(df4) 
    
    #getting our bins into lists, then append to a list of list, doing it with numpy array, will be efficient in time
    for i in range(start, len(df4), interval):
        a = df4["log2FoldChange"][i:(i+interval)].tolist()
        b = df4["pvalue"][i:(i+interval)].tolist()
        c = df4["baseMean"][i:(i+interval)].tolist()
        
        
        
        LFC.append(a)
        pVal.append(b)
        bMean.append(c)
    lastL = len(LFC[-1])
    for i in range(lastL):
        LFC[-2].append(LFC[-1][i])
        pVal[-2].append(pVal[-1][i])
        bMean[-2].append(bMean[-1][i])

    LFC.pop(-1)
    pVal.pop(-1)
    bMean.pop(-1)
    LFCBM =[]
    pv =[]
    for i in range(len(pVal)):
        pv.append(statistics.mean(pVal[i]))
    for i in range(len(LFC)):
        LFCBM.append(LFC[i])
        LFCBM.append(bMean[i])
    print("*"*100)
    print("Fasten your belts, we are plotting!")
    print("*"*100)
    plt.figure(figsize=(6,6))
    sns.scatterplot(data = pv)
    sns.despine(left=True, bottom=True, right = True, top = True)
    a = mpatches.Patch(color='Red', label='DownRegulated')
    b = mpatches.Patch(color='Green', label='Upregulated')
    plt.legend(handles=[a,b], loc = 'best', title = 'Gene Expression')
    #plt.figure(figsize=(5,5))
    plt.xlabel("Bins")
    plt.ylabel("Pvalue mean per bin")
    plt.title("Pvalue per gene per bin")
    plt.savefig(r"./Module1Graphs/PVV.png")

    plt.figure(figsize=(6,6))
    sns.boxplot(data = LFC,showfliers = False,showcaps = True, width = 0.6, notch = True, medianprops = {'color': 'red'}, showmeans = True, meanprops = {"markerfacecolor": "red"}, palette = 'RdYlGn_r', whiskerprops = dict(linestyle = '--', linewidth = 2));
    sns.despine(left=True, bottom=True, right = True, top = True)
    a = mpatches.Patch(color='Red', label='DownRegulated')
    b = mpatches.Patch(color='Green', label='Upregulated')
    plt.legend(handles=[a,b], loc = 'best', title = 'Gene Expression')
    #plt.figure(figsize=(5,5))
    plt.xlabel("Bins")
    plt.ylabel("Log2 Fold Change per gene")
    plt.title("Log2Fold Change per gene per bin")
    plt.savefig(r"./Module1Graphs/L2FC.png")
    plt.figure(figsize=(6,6)) 
    sns.boxplot(data= pVal, showfliers = False, showcaps = True, width =0.6 , notch = True, medianprops = {'color': 'red'}, showmeans =True, meanprops = {"markerfacecolor":"red"}, palette = 'RdYlGn_r', whiskerprops = dict(linestyle = '--', linewidth = 2 ));
    sns.despine(left=True, bottom=True, right = True, top = True)
    a = mpatches.Patch(color='Red', label='DownRegulated')
    b = mpatches.Patch(color='Green', label='Upregulated')
    plt.legend(handles=[a,b], loc = 'best', title = 'Gene Expression')

    plt.xlabel("Bins")
    plt.ylabel("pValue per gene")
    plt.title("pValue per gene per bin")
    plt.savefig(r"./Module1Graphs/PV.png")
    plt.figure(figsize=(6,6)) 
    sns.boxplot(data= bMean, showfliers = False,showcaps = True, width = 0.6, notch = True, medianprops = {'color': 'red'}, showmeans = True, meanprops = {"markerfacecolor":"red"}, palette = 'RdYlGn_r', whiskerprops = dict(linestyle ='--', linewidth =2));
    sns.despine(left=True, bottom=True, right = True, top = True)
    a = mpatches.Patch(color='Red', label='DownRegulated')
    b = mpatches.Patch(color='Green', label='Upregulated')
    plt.legend(handles=[a,b], loc = 'best', title = 'Gene Expression')

    plt.xlabel("Bins")
    plt.ylabel("Basemean per gene")
    plt.title("Basemean per gene per bin")
    plt.savefig(r"./Module1Graphs/BM.png")
    #this part is for grouping graphs in one sheet
    #fig, axes = plt.subplots(2, 2, figsize=(18, 10))

    #fig.suptitle('Graphs')
 
    
     
    #sns.boxplot(ax=axes[0, 0], data=LFC, showfliers = False,width = 0.4, notch = True,medianprops={"color": "coral"}, showmeans = True, meanprops = {"markerfacecolor": "red"},palette = "Greens", whiskerprops = dict(linestyle = "--", linewidth = 2))
    #sns.boxplot(ax=axes[0, 1], data=pVal, showfliers = True,width = 0.4, notch = True, palette = "Greens",medianprops={"color": "coral"} ,showmeans = True, meanprops = {"markerfacecolor":"red"}, whiskerprops = dict(linestyle = "--", linewidth = 2 ))
    #sns.boxplot(ax=axes[1, 0], data=bMean, showfliers = False, width = 0.4, notch = True, palette= 'Greens', medianprops={"color": "coral"},showmeans = True, meanprops = {"markerfacecolor": "red"}, whiskerprops = dict(linestyle = "--", linewidth = 2))
    
    #fig.show()
    
    #plt.savefig("sm.png")
    #plt.show()
    

#before this step, bedgraph files were thrown to SEACR for peak calling, then the peaks were annotated in Galaxy using ChipSeeker.
#will see how to integrate the scripts for peak calling and annotation in my code
#Please save your annotated peak file, under the format "CHDx_xIGGx.csv" where x is the number of your sample!

def peaksGenes(annot,bin):
    
    print("PPP: Peaks Plotting Part!")
    print("*"*100)
    print("*"*100)
    #listFreq = []
    #file = ['CHD3_1IGG1.csv','CHD3_2IGG2.csv','CHD3_3IGG3.csv','CHD4_1IGG1.csv','CHD4_2IGG2.csv','CHD4_3IGG3.csv']
    #file = input('Can you please input your annotated file:')
    #for i in range(6):
    listFreq = []
        #file = ['CHD3_1IGG1.csv','CHD3_2IGG2.csv','CHD3_3IGG3.csv','CHD4_1IGG1.csv','CHD4_2IGG2.csv','CHD4_3IGG3.csv']
    df1 = pd.read_csv(annot, sep=",")
    pd.set_option("display.width", None)
        #df1.to_csv('annotatedFile.csv')
        #df2 = pd.read_csv('annotatedFile.csv')   
    

        #it calculates the number of peaks per gene
    df1 = df1["SYMBOL"].value_counts()
    df1 = df1.reset_index()
    df1.rename(columns = {"index":"GeneId", "SYMBOL":"Freq1"}, inplace = True)
    df1.to_csv(r"./Module1Files/GenePeaks.csv", index = False)
        
    for i in range(bin):
        df3 = pd.read_csv("./Module1Files/GeneBinFreq"+str(i)+".csv")
        df4 = pd.read_csv("./Module1Files/GenePeaks.csv")
        df3.rename(columns={"index":"GeneId", "GeneId":"Freq2"},inplace = True)
        df3.drop(["Freq2"],axis = 1,inplace=True)
        df5 = pd.merge(df3,df4, on = "GeneId", how = "left")    
        
        df5 = df5.fillna(0)
            #df = df.dropna()
        df5 = df5.drop(df5[df5['Freq1'] <= 0].index)
            #sorting our basemean

        df5 = df5.reset_index(drop = True)
        df5.to_csv(r"./Module1Files/Lastbin"+str(i)+".csv")
        a = df5["Freq1"].tolist()
        listFreq.append(a)
          
    
    for i in range(bin):
        df6 = pd.read_csv("./Module1Files/Lastbin"+str(i)+".csv", sep = ',')
   
    
    
        df6 = df6.fillna(0)
            #df = df.dropna()
        df6 = df6.drop(df6[df6['Freq1'] <= 0].index)
            #sorting our basemean

        df6 = df6.reset_index(drop = True)
        a = df6['Freq1'].tolist()

        b = statistics.mean(a)
        lis.append(b)

      
    plt.figure(figsize=(6,6))
    sns.boxplot(data = listFreq, showfliers = False, showcaps = True,width = 0.6, notch = True, palette = 'RdYlGn_r', medianprops={"color": "green"},showmeans = True, meanprops = {"markerfacecolor": "green"}, whiskerprops = dict(linestyle = "--", linewidth = 2));
    sns.despine(left=True, bottom=True, right = True, top = True)
    a = mpatches.Patch(color='Red', label='DownRegulated')
    b = mpatches.Patch(color='Green', label='Upregulated')
    plt.legend(handles=[a,b], loc = 'best', title = 'Gene Expression')


    plt.xlabel("Bins")
    plt.ylabel("Num of peaks/gene")
    plt.title("Number of peaks/gene per bin")
    plt.savefig(r"./Module1Graphs/Npeaks.png")
    
    
    

def transcriptL(quant,bin):
    df = pd.read_csv('newsmth.csv', delimiter = '\t')
    df.rename(columns= {'Name':'GeneId'} , inplace = True)
      
    df1 = pd.read_csv("./Module1Files/finalModified.csv", sep=",")
    pd.set_option("display.width",None)
    pd.set_option("display.max_rows", None)


    df["TPM*EffL"] = df["TPM"]*df["EffectiveLength"]
    df2 = df.groupby("GeneId").sum()

    df2.reset_index(inplace =True)
    df2["WELG"] = df2["TPM*EffL"]/df2["TPM"]
    
    df1.drop(['baseMean', 'log2FoldChange','pvalue'], axis=1, inplace = True)

    df2.drop([ 'Length', 'EffectiveLength', 'TPM', 'NumReads','TPM*EffL'], axis=1, inplace= True)
      
    df3 = pd.merge(df1,df2, on="GeneId", how="inner")


    df4 = df3.dropna()
    df4.to_csv("FinalTr.csv", index = False)


    interval = int(len(df4)/bin)
    listt = []


    for i in range(0,len(df4), interval):
        a = df4['WELG'][i:(i+interval)].tolist()
        listt.append(a)

    #for i in range(len(listt[-1])):
        #listt[-2].append(listt[-1][i])

    
    #listt.pop(-1)
    for i in range(len(listt)):
        print(len(listt[i]))
    
    #lisi =[]
    #for i in range(len(listt)):
        #Mlis.append(statistics.mean(listt[i]))
    #for i in range(bin):
        #lisi.append(i)
   
    
    
    #plt.plot(lisi, Mlis, c = "#FF0000")
    plt.figure(figsize=(8,8)) 
    sns.boxplot(data = listt, showfliers = False,showcaps = True, width = 0.6, notch = True, palette = 'RdYlGn_r', medianprops={"color": "coral"},showmeans = True, meanprops = {"markerfacecolor": "red"}, whiskerprops = dict(linestyle = "--", linewidth = 2));
    sns.despine(left=True, bottom=True, right = True, top = True)
    a = mpatches.Patch(color='Red', label='DownRegulated')
    b = mpatches.Patch(color='Green', label='Upregulated')
    plt.legend(handles=[a,b], loc = 'best', title = 'Gene Expression')
    plt.xlabel("Bins")
    plt.ylabel("Weighted effective transcipt length/gene")
    plt.title("Weighted effective transcript length/gene per bin")
    plt.savefig(r"./Module1Graphs/WELG.png")
    
    

def GeneL(annot,bin):

    # This script has to be run on every annotated sample, example : CHD3_1IGG1.csv, CHD4_1IGG1.csv


    #annot = input("Please input your annotated peaks file, in the format CHDx_xIGGx.csv:")
    df = pd.read_csv(annot, sep = ',')
    pd.set_option("display.max_rows", None)
    pd.set_option("display.width", None)
    
    # We take the longest annotation for the gene

    df1 = df.groupby('SYMBOL')['geneLength'].max()
    df1 = df1.reset_index()
    df1.to_csv(r"./Module1Files/GeneL.csv", index = False)

    df2 = pd.read_csv("./Module1Files/finalModified.csv", sep = ",")
    df3 = pd.read_csv("./Module1Files/GeneL.csv", sep = ",")

    df2.drop(['baseMean', 'log2FoldChange','pvalue'], axis = 1, inplace = True)
    df3.rename(columns = {'SYMBOL': 'GeneId'}, inplace = True)

    df4 = pd.merge(df2,df3, on = "GeneId", how = "inner")

    #lis = []
    #sil = []
    interval = int(len(df4)/bin)


    for i in range(0, len(df4), interval):
        a = df4["geneLength"][i:(i+interval)].tolist()
        sil.append(a)

    #for i in range(len(sil[-1])):
        #sil[-2].append(sil[-1][i])

    #sil.pop(-1)
    for i in range(len(sil)):
        print(len(sil[i]))    
    ren = []
    #for i in range(len(sil)):
        #mLis.append(statistics.mean(sil[i]))
    #for i in range(bin):
        #ren.append(i)
    
    #plt.plot(ren,mLis, c ="#FF0000")
    plt.figure(figsize=(8,8)) 
    sns.boxplot(data = sil,showmeans = True, showfliers = False,showcaps = True,palette = "RdYlGn_r",width = 0.6,notch = True,medianprops={"color": "coral"},meanprops = {"markerfacecolor": "red"} ,whiskerprops = dict(linestyle = "--", linewidth = 2));
    sns.despine(left=True, bottom=True, right = True, top = True)
    a = mpatches.Patch(color='Red', label='DownRegulated')
    b = mpatches.Patch(color='Green', label='Upregulated')
    plt.legend(handles=[a,b], loc = 'best', title = 'Gene Expression')
    plt.ylabel("Gene Length/gene")
    plt.xlabel("Bins")
    plt.title("Gene Length/gene per bin")
    plt.savefig(r"./Module1Graphs/GL.png")
    

def TSS(annot, bin):

    
    d = pd.read_csv(annot, sep = ',')

    d.drop(d.columns.difference(['SYMBOL', 'distanceToTSS', 'annotation']), axis =1, inplace = True)

    d['ann'] = d['annotation'].str.split("(", expand = True)[0]
    d.drop(['annotation'], axis = 1,inplace = True )
    d.rename(columns = {'SYMBOL': 'GeneId'}, inplace = True)


    #bin = int(input('How many bins are you using:'))


    maxi =  []
    aa =  []
    bb = []
    mini =  []
    mean =  []
    Ratio = []
    for i in range(bin):
        o = pd.read_csv('./Module1Files/Lastbin'+str(i)+'.csv')
        o.drop(['Unnamed: 0', 'Freq1'], axis =1, inplace = True)
        dd = pd.merge(o,d, on = 'GeneId', how = 'inner')
        dd1 = dd.drop(dd[dd['ann'] != 'Distal Intergenic'].index)
        dd1.reset_index(drop = True)
        dd2 = dd.drop(dd[dd['ann'] != 'Intron '].index)
        dd2.reset_index(drop = True)
        dd3 = dd.drop(dd[dd['ann'] != 'Promoter '].index)
        dd3.reset_index(drop = True)
        dd1.rename(columns= {'ann': 'Distal Intergenic R'}, inplace = True)
        dd2.rename(columns= {'ann': 'Intron R'}, inplace = True)
        dd3.rename(columns = {'ann': 'Promoter R'}, inplace = True)
        dd1.drop(['distanceToTSS'], axis= 1, inplace = True)
        dd2.drop(['distanceToTSS'], axis= 1, inplace = True)
        dd3.drop(['distanceToTSS'], axis =1, inplace = True)
        dd4 = dd1.groupby('GeneId')['Distal Intergenic R'].count()
        dd5 = dd2.groupby('GeneId')['Intron R'].count()
        dd6 = dd3.groupby('GeneId')['Promoter R'].count()
        dd7 = pd.merge(dd4,dd5, on = 'GeneId', how = 'outer')
        dd8 = pd.merge(dd7, dd6, on = 'GeneId', how = 'outer')
        #dd6 = dd5.fillna(0)
        dd9 = dd8.fillna(0)
        dd10 = dd9[(dd9['Intron R'] == 0) & (dd9['Distal Intergenic R'] == 0)]
        dd11 = dd9[(dd9['Intron R'] == 0) & (dd9['Distal Intergenic R'] != 0)]
        dd12 = dd9[(dd9['Intron R'] != 0) & (dd9['Distal Intergenic R'] == 0)]
        dd13 = dd9[(dd9['Intron R'] != 0) & (dd9['Distal Intergenic R'] != 0)]
        dd13 ['Ratio'] = dd13['Intron R']/dd13['Distal Intergenic R']
        a = dd13['Ratio'].tolist()
        Ratio.append(a)
    plt.figure(figsize=(4,6)) 
    sns.boxplot(data = Ratio,showcaps = True, width = 0.6, showfliers = False, showmeans = True, notch = True, whiskerprops = dict(linestyle = '--'), color = 'blue')
    plt.xlabel('Bins')
    plt.title('Nbr of intronic peaks per gene/Number of distal Intergenic peaks per gene')
    plt.savefig(r'./Module3Graphs/#intronicPks#DisIntergenicPks.svg')
    

def lastP(annot, bin):
    #bin = int(input('What is your number of bins:'))

    d = pd.read_csv(annot, sep = ',')
    pd.set_option('display.max_rows',None) 
    pd.set_option('display.width',None)

    d['ann'] = d['annotation'].str.split("(", expand = True)[0]

    #d.drop(['Chrom', 'Start', 'End', 'chr1.3000700.3000900','geneChr','transcriptId','geneStrand', 'geneStart', 'geneEnd', 'geneName'], axis = 1, inplace = True)
    #d.drop(d.columns.difference(['SYMBOL', 'ann','distanceToTSS', 'X3.88']), axis =1, inplace = True)
    d.drop(d.columns[[0,1,2,3,4,5,6,8,10,11,12,13,14,15,18]],axis = 1, inplace = True)
    d.rename(columns = {'SYMBOL': 'geneId'}, inplace = True)
    
    d1 = d.drop(d[d['ann'] != 'Intron '].index)
    d1.reset_index(drop = True)
    d2 = d1.groupby('geneId')['ann'].count()
    ddd = d1.groupby('geneId')[d1.columns[0]].sum()
    d21 = d1.groupby('geneId')[d1.columns[0]].max()
    d2.reset_index(drop = True)


    dd = d.drop(d[d['distanceToTSS'] > 0].index)
    dd.reset_index(drop = True)
    d3 = dd.drop(dd[dd['ann'] != 'Distal Intergenic'].index)


    d3.reset_index(drop = True)
    d3.to_csv('data1.csv', index = False)
    d4 = d3.groupby('geneId')['ann'].count()
    #was mean before
    print(d3.columns)
    dddd = d3.groupby('geneId')[d3.columns[0]].sum()
    #max peak depth upstream
    ddddd = d3.groupby('geneId')[d3.columns[0]].max()

    d7 = d.drop(d[d['distanceToTSS'] < 0].index)
    d7.reset_index(drop = True)
    d8 = d7.drop(d7[d7['ann'] != 'Distal Intergenic'].index)
    d8.reset_index(drop = True)
    d9 = d8.groupby('geneId')['ann'].count()
    #was mean before
    d10 = d8.groupby('geneId')[d8.columns[0]].sum()
    #max peak depth downstream
    dddddd = d8.groupby('geneId')[d8.columns[0]].max()

    d11 = d.drop(d[d['ann'] != 'Promoter '].index)
    d11.reset_index(drop = True)

    d12 = d11.groupby('geneId')['ann'].count()

    d13 = d11.groupby('geneId')[d11.columns[0]].max()
    d14 = d11.groupby('geneId')[d11.columns[0]].sum()

    d21.to_csv(r'./Module4Files/IntronMX.csv')
    d2.to_csv(r'./Module4Files/countIntrons.csv')
    ddd.to_csv(r'./Module4Files/SumIntrons.csv')
    d4.to_csv(r'./Module4Files/countDistalUp.csv')
    dddd.to_csv(r'./Module4Files/SumDistalUp.csv')
    d9.to_csv(r'./Module4Files/countDistalDown.csv')
    d10.to_csv(r'./Module4Files/SumDistalDown.csv')
    d12.to_csv(r'./Module4Files/countPromoters.csv')
    d13.to_csv(r'./Module4Files/meanPromoters.csv')
    d14.to_csv(r'./Module4Files/sumPromDep.csv')
    ddddd.to_csv(r'./Module4Files/upstreamMX.csv')
    dddddd.to_csv(r'./Module4Files/downstreamMX.csv')

    d22 = pd.read_csv('./Module4Files/IntronMX.csv')
    s = pd.read_csv('./Module4Files/sumPromDep.csv',sep = ',')
    v = pd.read_csv('./Module4Files/upstreamMX.csv', sep = ',')
    h = pd.read_csv('./Module4Files/downstreamMX.csv', sep = ',')
    x = pd.read_csv('./Module4Files/countIntrons.csv', sep = ',')
    y = pd.read_csv('./Module4Files/SumIntrons.csv', sep = ',')
    a = pd.read_csv('./Module4Files/countDistalUp.csv', sep = ',')
    b = pd.read_csv('./Module4Files/SumDistalUp.csv', sep = ',')
    c = pd.read_csv('./Module4Files/countDistalDown.csv', sep = ',')
    e = pd.read_csv('./Module4Files/SumDistalDown.csv', sep = ',')
    j = pd.read_csv('./Module4Files/countPromoters.csv', sep = ',')
    k = pd.read_csv('./Module4Files/meanPromoters.csv', sep = ',')
    d22.rename(columns={d22.columns[1]:'MaxIntrPeak'}, inplace = True)
    s.rename(columns={s.columns[1]:'SumPromDep'}, inplace = True)
    v.rename(columns={v.columns[1]:'MaxUpPeak'},inplace = True)
    h.rename(columns={h.columns[1]:'MaxDownPeak'},inplace = True)
    j.rename(columns={'ann':'#PromotersPeaks'},inplace = True)
    k.rename(columns={k.columns[1]:'MaxDepthProm'},inplace = True)
    x.rename(columns={'ann':'#IntronsPeaks'},inplace = True)
    y.rename(columns={y.columns[1]:'Sum_Dpth_Pks_Introns'},inplace = True)
    a.rename(columns={'ann':'#DistalUpPeaks'},inplace = True)
    b.rename(columns={b.columns[1]:'Sum_Dpth_Pks_D_Up'},inplace = True)
    c.rename(columns={'ann':'#DistalDownPeaks'},inplace = True)
    e.rename(columns={e.columns[1]:'Sum_Dpth_Pks_D_Down'},inplace = True)

    w = pd.merge(v,h, on = 'geneId', how='outer')
    u = pd.merge(j,k, on='geneId', how ='outer')
    t = pd.merge(x,y, on= 'geneId', how = 'outer')
    p = pd.merge(a,b, on= 'geneId', how = 'outer')
    o = pd.merge(c,e, on= 'geneId', how = 'outer')


    r = pd.merge(t,p, on='geneId', how='outer')
    n = pd.merge(r,o, on='geneId', how='outer')
    m = pd.merge(n,u, on='geneId', how='outer')
    l = pd.merge(w,m, on = 'geneId', how='outer')
    d20 = pd.merge(l,s, on ='geneId', how = 'outer')

    d100 = d3.groupby(['geneId','distanceToTSS'])[d3.columns[0]].max()
    d101 = d8.groupby(['geneId','distanceToTSS'])[d8.columns[0]].max()
    d100.to_csv(r'./Module4Files/maxDown.csv')
    d101.to_csv(r'./Module4Files/maxUP.csv')

    d1000 = pd.read_csv('./Module4Files/maxUP.csv')
    d1001 = pd.read_csv('./Module4Files/maxDown.csv')

    d1000.rename(columns={d1000.columns[2]:'Depth'},inplace = True)
    d1001.rename(columns={d1001.columns[2]:'Depth'},inplace = True)



    a = d1000.groupby('geneId')['distanceToTSS'].apply(list).tolist()
    b = d1000.groupby('geneId')['Depth'].apply(list).tolist()
    c = d1001.groupby('geneId')['distanceToTSS'].apply(list).tolist()
    d = d1001.groupby('geneId')['Depth'].apply(list).tolist() 


    list3 = []
    list4 = []

    list1 = []
    list2 = []
    for i in range(len(a)):
        listt = zip(b[i],a[i])
        maxa,maxb = max(list(listt))
        list1.append(maxa)
        list2.append(maxb)

    for i in range(len(c)):
        listt = zip(d[i],c[i])
        maxa,maxb = max(list(listt))
        list3.append(maxa)
        list4.append(maxb)

    #aa = pd.read_csv('ar.csv')
    tt= pd.read_csv('./Module4Files/countDistalUp.csv')
    pp = pd.read_csv('./Module4Files/countDistalDown.csv')
    pp['DownDTSS']= list2
    tt['UpDTSS']= list4
    ii = pd.merge(pp,tt, on='geneId', how ='outer')
    mm = pd.merge(ii,d20,on='geneId',how='outer')
    gg = pd.merge(mm,d22, on='geneId', how = 'outer')
    #zz = gg.fillna(0)
    zz = gg
    zz.rename(columns={'geneId':'GeneId'}, inplace = True)
    #zz['Sum_Dpth_Pks_D_Up'] *= -1 
    #zz['MaxUpPeak'] *= -1
    zz['SumDisUp/sumPromDep'] = 1/(zz['SumPromDep']/zz['Sum_Dpth_Pks_D_Up'])
    zz['SumDisDown/sumPromDep'] = 1/(zz['SumPromDep']/zz['Sum_Dpth_Pks_D_Down'])
    zz['Mx_D_UpPeak/MxProm'] = 1/(zz['MaxUpPeak']/zz['MaxDepthProm'])
    zz['Mx_D_DownPeak/MxProm'] = 1/(zz['MaxDepthProm']/zz['MaxDownPeak'])
    zz['Mx_IntronPeak/MxProm'] = 1/(zz['MaxDepthProm']/zz['MaxIntrPeak'])
    zz['Sum_D_Up/Sum_Introns']= zz['Sum_Dpth_Pks_D_Up']/zz['Sum_Dpth_Pks_Introns']
    zz['Sum_D_Up/Sum_D_Down']= zz['Sum_Dpth_Pks_D_Up']/zz['Sum_Dpth_Pks_D_Down']
    zz.replace([np.inf, -np.inf], np.nan, inplace = True)
    zz['AbsUpDTSS'] = zz['UpDTSS']*(-1)



    #zz.fillna(0,inplace = True)
    #kk = zz.dropna()
    for i in range(bin):
        nn = pd.read_csv("./Module1Files/finalModified.csv")
        r = pd.read_csv('./Module1Files/Lastbin'+str(i)+'.csv')
        r.drop(['Unnamed: 0', 'Freq1'],axis = 1, inplace = True)
        toto = pd.merge(r,zz,on='GeneId', how='inner')
        pop = pd.merge(nn, toto, on= 'GeneId', how = 'inner' )
        pop.drop(['ann_x', 'ann_y','baseMean','pvalue'],axis = 1, inplace = True) 
        pop['AbsLog2FC'] = pop['log2FoldChange']
        pop.to_csv(r'./Module4Files/rat'+str(i)+'.csv', index = False)
    listtt = []
    listtt1 = []
    listtt2 = []
    listtt3 = []
    listtt4 = []
    listtt5 = []
    listtt6 = []
    listtt7 = []
    listtt8 = []
    for i in range(bin):
        a = pd.read_csv('./Module4Files/rat'+str(i)+'.csv')
        b = a['SumDisUp/sumPromDep'].tolist()
        c = a['SumDisDown/sumPromDep'].tolist()
        d = a['Mx_D_UpPeak/MxProm'].tolist()
        e = a['Mx_D_DownPeak/MxProm'].tolist()
        f = a['Mx_IntronPeak/MxProm'].tolist()
        g = a['Sum_D_Up/Sum_Introns'].tolist()
        h = a['Sum_D_Up/Sum_D_Down'].tolist()
        aa = a['AbsUpDTSS'].tolist()
        bb = a['#DistalUpPeaks'].tolist()
        listtt.append(b)
        listtt1.append(c)
        listtt2.append(d)
        listtt3.append(e)
        listtt4.append(f)
        listtt5.append(g)
        listtt6.append(h) 
        listtt7.append(aa)
        listtt8.append(bb)
        a.dropna(inplace = True)
    plt.figure(figsize=(4,6)) 
    sns.boxplot(data = listtt, showcaps = True,width = 0.6, notch = True, showfliers = False,medianprops = {'color': 'red'}, showmeans =True, meanprops = {"markerfacecolor":"red"}, color = 'blue', whiskerprops = dict(linestyle = '--', linewidth = 2 ))
    plt.title('Sum of Distal Upstream intergenic depths(per gene)/Sum of promoter of depths(per gene)')
    plt.savefig(r'./Module4Graphs/SumDUPdepths(per gene)-Sumpromoterdepths(per gene).svg')
    
    plt.figure(figsize=(4,6)) 
    sns.boxplot(data = listtt7, showcaps = True,width = 0.6, notch = True, showfliers = False,medianprops = {'color': 'red'}, showmeans =True, meanprops = {"markerfacecolor":"red"}, color = 'blue', whiskerprops = dict(linestyle = '--', linewidth = 2 ))
    plt.title('Distal Upstream Enhancers distance to TSS per gene/bin')
    plt.savefig(r'./Module4Graphs/DUPdistanceTSSgene.svg')
    
    plt.figure(figsize=(4,6)) 
    sns.boxplot(data = listtt8, showcaps = True,width = 0.6, notch = True, showfliers = False,medianprops = {'color': 'red'}, showmeans =True, meanprops = {"markerfacecolor":"red"}, color = 'blue', whiskerprops = dict(linestyle = '--', linewidth = 2 ))
    plt.title('Distal Upstream Enhancers number of peaks per gene/bin')
    plt.savefig(r'./Module4Graphs/DUpNbrpksgene.svg')
    plt.show()
    #sns.boxplot(data = listtt1, showcaps = False,width =0.4, notch = True, showfliers = False,medianprops = {'color': 'red'}, showmeans =True, meanprops = {"markerfacecolor":"red"}, palette = 'RdYlGn_r', whiskerprops = dict(linestyle = '--', linewidth = 2 ))
    #plt.title('Sum of promoter of depths(per gene) / Sum of Distal Downstream intergenic depths(per gene)')
    #plt.show()
    plt.figure(figsize=(4,6)) 
    sns.boxplot(data = listtt2, showcaps = True,width = 0.6, notch = True, showfliers = False,medianprops = {'color': 'red'}, showmeans =True, meanprops = {"markerfacecolor":"red"}, color = 'blue', whiskerprops = dict(linestyle = '--', linewidth = 2 ))
    plt.title('Highest D_Upstream_Peak(per gene)/Highest promoter peak(per gene)')
    plt.savefig(r'./Module4Graphs/HighestDUPpeakGENE-HighestPrompeakGENE.svg')
    

def flavours(bin):
    pd.set_option('display.max_rows', None)
    pd.set_option('display.width', None)
    #bin = int(input('What is your desired number of bins:'))

    lis = [1,2,3,4,5,6,7]

    for i in range(len(lis)):
        inn = lis[i]

        if inn == 1:
            for i in range(bin):   
                df = pd.read_csv('./Module4Files/rat'+str(i)+'.csv', sep =',')
                df.drop(['MaxUpPeak','MaxDownPeak', 'Sum_Dpth_Pks_Introns','Sum_Dpth_Pks_D_Down', 'MaxDepthProm',
                'MaxIntrPeak', 'SumDisUp/sumPromDep', 'SumDisDown/sumPromDep',
                'Mx_D_UpPeak/MxProm', 'Mx_D_DownPeak/MxProm', 'Mx_IntronPeak/MxProm',
                'Sum_D_Up/Sum_Introns', 'Sum_D_Up/Sum_D_Down'], axis = 1, inplace = True)
                df0 = df.fillna(0)
                df00 = df0.drop(df0[df0['#PromotersPeaks'] == 0].index)
                df00.reset_index(drop = True)
                df1 = df00.drop(df00[df00['#IntronsPeaks'] != 0].index)
                df1.reset_index(drop = True)
                df2 = df1.drop(df1[df1['#DistalUpPeaks'] != 0].index)
                df2.reset_index(drop = True)
                df3 = df2.drop(df2[df2['#DistalDownPeaks'] != 0].index)
                df3.reset_index(drop = True)
                print("Number of genes [Promoter peaks only]:",len(df3))
                df3.to_csv(r'./Module5Files/Promoteronly'+str(i)+'.csv')
            

        elif inn == 2: 
            for i in range(bin):   
                df = pd.read_csv('./Module4Files/rat'+str(i)+'.csv', sep =',')
                df.drop(['MaxUpPeak','MaxDownPeak', 'Sum_Dpth_Pks_Introns','Sum_Dpth_Pks_D_Down', 'MaxDepthProm',
                'MaxIntrPeak', 'SumDisUp/sumPromDep', 'SumDisDown/sumPromDep',
                'Mx_D_UpPeak/MxProm', 'Mx_D_DownPeak/MxProm', 'Mx_IntronPeak/MxProm',
                'Sum_D_Up/Sum_Introns', 'Sum_D_Up/Sum_D_Down'], axis = 1, inplace = True)
                df0 = df.fillna(0)
                df00 = df0.drop(df0[df0['#PromotersPeaks'] == 0].index)
                df00.reset_index(drop = True)
                df1 = df00.drop(df00[df00['#IntronsPeaks'] == 0].index)
                df1.reset_index(drop = True)
                df2 = df1.drop(df1[df1['#DistalUpPeaks'] != 0].index)
                df2.reset_index(drop = True)
                df3 = df2.drop(df2[df2['#DistalDownPeaks'] != 0].index)
                df3.reset_index(drop = True)
                print("Number of genes [Promoter + Intronic only]:",len(df3))
                df3.to_csv(r'./Module5Files/Promoter&Intronic'+str(i)+'.csv')
        elif inn == 3: 
            for i in range(bin):   
                df = pd.read_csv('./Module4Files/rat'+str(i)+'.csv', sep =',')
                df.drop(['MaxDownPeak', 'Sum_Dpth_Pks_Introns','Sum_Dpth_Pks_D_Down', 'MaxDepthProm',
                'MaxIntrPeak', 'SumDisDown/sumPromDep', 'Mx_D_DownPeak/MxProm', 'Mx_IntronPeak/MxProm',
                'Sum_D_Up/Sum_Introns', 'Sum_D_Up/Sum_D_Down'], axis = 1, inplace = True)
                df0 = df.fillna(0)
                df00 = df0.drop(df0[df0['#PromotersPeaks'] == 0].index)
                df00.reset_index(drop = True)
                df1 = df00.drop(df00[df00['#IntronsPeaks'] != 0].index)
                df1.reset_index(drop = True)
                df2 = df1.drop(df1[df1['#DistalUpPeaks'] == 0].index)
                df2.reset_index(drop = True)
                df3 = df2.drop(df2[df2['#DistalDownPeaks'] != 0].index)
                df3.reset_index(drop = True)
                print("Number of genes [Promoter + DistalUp only]:",len(df3))
                df3.to_csv(r'./Module5Files/Promoter&DistalUp'+str(i)+'.csv')
        elif inn == 4: 
            for i in range(bin):   
                df = pd.read_csv('./Module4Files/rat'+str(i)+'.csv', sep =',')
                df.drop(['MaxUpPeak','MaxDownPeak', 'Sum_Dpth_Pks_Introns','Sum_Dpth_Pks_D_Down', 'MaxDepthProm',
                'MaxIntrPeak', 'SumDisUp/sumPromDep', 'SumDisDown/sumPromDep',
                'Mx_D_UpPeak/MxProm', 'Mx_D_DownPeak/MxProm', 'Mx_IntronPeak/MxProm',
                'Sum_D_Up/Sum_Introns', 'Sum_D_Up/Sum_D_Down'], axis = 1, inplace = True)
                df0 = df.fillna(0)
                df00 = df0.drop(df0[df0['#PromotersPeaks'] == 0].index)
                df00.reset_index(drop = True)
                df1 = df00.drop(df00[df00['#IntronsPeaks'] != 0].index)
                df1.reset_index(drop = True)
                df2 = df1.drop(df1[df1['#DistalUpPeaks'] != 0].index)
                df2.reset_index(drop = True)
                df3 = df2.drop(df2[df2['#DistalDownPeaks'] == 0].index)
                df3.reset_index(drop = True)
                print("Number of genes [Promoter + DistalDown only]:",len(df3))
                df3.to_csv(r'./Module5Files/Promoter&DistalDown'+str(i)+'.csv')
        elif inn == 5: 
            for i in range(bin):   
                df = pd.read_csv('./Module4Files/rat'+str(i)+'.csv', sep =',')
                df.drop(['MaxUpPeak','MaxDownPeak', 'Sum_Dpth_Pks_Introns','Sum_Dpth_Pks_D_Down', 'MaxDepthProm',
                'MaxIntrPeak', 'SumDisUp/sumPromDep', 'SumDisDown/sumPromDep',
                'Mx_D_UpPeak/MxProm', 'Mx_D_DownPeak/MxProm', 'Mx_IntronPeak/MxProm',
                'Sum_D_Up/Sum_Introns', 'Sum_D_Up/Sum_D_Down'], axis = 1, inplace = True)
                df0 = df.fillna(0)
                df00 = df0.drop(df0[df0['#PromotersPeaks'] == 0].index)
                df00.reset_index(drop = True)
                df1 = df00.drop(df00[df00['#IntronsPeaks'] != 0].index)
                df1.reset_index(drop = True)
                df2 = df1.drop(df1[df1['#DistalUpPeaks'] == 0].index)
                df2.reset_index(drop = True)
                df3 = df2.drop(df2[df2['#DistalDownPeaks'] == 0].index)
                df3.reset_index(drop = True)
                print("Number of genes [Promoter + DistalUp + DistalDown]:",len(df3))
                df3.to_csv(r'./Module5Files/Promoter&DUP&dDOWN'+str(i)+'.csv')
        elif inn == 6: 
            for i in range(bin):   
                df = pd.read_csv('./Module4Files/rat'+str(i)+'.csv', sep =',')
                df.drop(['MaxUpPeak','MaxDownPeak', 'Sum_Dpth_Pks_Introns','Sum_Dpth_Pks_D_Down', 'MaxDepthProm',
                'MaxIntrPeak', 'SumDisUp/sumPromDep', 'SumDisDown/sumPromDep',
                'Mx_D_UpPeak/MxProm', 'Mx_D_DownPeak/MxProm', 'Mx_IntronPeak/MxProm',
                'Sum_D_Up/Sum_Introns', 'Sum_D_Up/Sum_D_Down'], axis = 1, inplace = True)
                df0 = df.fillna(0)
                df00 = df0.drop(df0[df0['#PromotersPeaks'] == 0].index)
                df00.reset_index(drop = True)
                df1 = df00.drop(df00[df00['#IntronsPeaks'] == 0].index)
                df1.reset_index(drop = True)
                df2 = df1.drop(df1[df1['#DistalUpPeaks'] == 0].index)
                df2.reset_index(drop = True)
                df3 = df2.drop(df2[df2['#DistalDownPeaks'] != 0].index)
                df3.reset_index(drop = True)
                print("Number of genes [Promoter + Intronic + DistalUp]:",len(df3))
                df3.to_csv(r'./Module5Files/Promoter&Intronic&DUP'+str(i)+'.csv')
        elif inn == 7: 
            for i in range(bin):   
                df = pd.read_csv('./Module4Files/rat'+str(i)+'.csv', sep =',')
                df.drop(['MaxUpPeak','MaxDownPeak', 'Sum_Dpth_Pks_Introns','Sum_Dpth_Pks_D_Down', 'MaxDepthProm',
                'MaxIntrPeak', 'SumDisUp/sumPromDep', 'SumDisDown/sumPromDep',
                'Mx_D_UpPeak/MxProm', 'Mx_D_DownPeak/MxProm', 'Mx_IntronPeak/MxProm',
                'Sum_D_Up/Sum_Introns', 'Sum_D_Up/Sum_D_Down'], axis = 1, inplace = True)
                df0 = df.fillna(0)
                df00 = df0.drop(df0[df0['#PromotersPeaks'] == 0].index)
                df00.reset_index(drop = True)
                df1 = df00.drop(df00[df00['#IntronsPeaks'] == 0].index)
                df1.reset_index(drop = True)
                df2 = df1.drop(df1[df1['#DistalUpPeaks'] != 0].index)
                df2.reset_index(drop = True)
                df3 = df2.drop(df2[df2['#DistalDownPeaks'] == 0].index)
                df3.reset_index(drop = True)
                print("Number of genes [Promoter + Intronic + DistalDown]:",len(df3))
                df3.to_csv(r'./Module5Files/Promoter&Intronic&DDOWN'+str(i)+'.csv')

def corre(bin):

    list1 = []
    list2 = []
    #list4 = ['0','1','2','3','4','5','6','7','8','9','10','11']
    list3 =[]
    list5 =[]
    for i in range(bin):
        d = pd.read_csv('./Module5Files/Promoter&DistalUp'+str(i)+'.csv')
        d.drop(['Unnamed: 0','#DistalUpPeaks', 'GeneId', 'DownDTSS', '#IntronsPeaks', '#DistalDownPeaks',
               '#PromotersPeaks', 'AbsUpDTSS', 'AbsLog2FC'], axis = 1, inplace = True)
        stat, pv = stats.spearmanr(d['log2FoldChange'],d['Mx_D_UpPeak/MxProm'])
        rst = round(stat,4)
        rpv = round(pv,4)
        sns.lmplot(data = d,x='Mx_D_UpPeak/MxProm',y='log2FoldChange')
        a = mpatches.Patch(color = 'White', label = 'Correlation coef r: '+str(rst))
        b = mpatches.Patch(color = 'White', label = 'pValue: '+str(rpv))
        plt.legend(handles =[a,b] ,loc = 'best', title = 'Spearman correlation test')  
        plt.savefig(r'./Module6Graphs/CorrelationBin'+str(i)+'.svg')
        
        d['UpDTSS'] *= -1
        a = d['UpDTSS'].tolist()
        list1.append(a)
        c = d['MaxUpPeak'].tolist()
        list2.append(c)
        e = d['Mx_D_UpPeak/MxProm'].tolist()
        list3.append(e)
        h = d['SumDisUp/sumPromDep'].tolist()
        list5.append(h) 
    plt.figure(figsize=(4,6)) 
    sns.boxplot(data = list1,showfliers = False,showcaps = True, width = 0.6, medianprops = {'color': 'red'}, showmeans = True, meanprops = {"markerfacecolor": "red"}, color = 'blue', whiskerprops = dict(linestyle = '--', linewidth = 2));
    plt.title("Average Distance to TSS")
    plt.savefig(r'./Module6Graphs/AvgDisttoTSS.svg')
    plt.figure(figsize=(4,6)) 
    sns.boxplot(data = list2,showfliers = False,showcaps = True, width = 0.6, medianprops = {'color': 'red'}, showmeans = True, meanprops = {"markerfacecolor": "red"}, color = 'blue', whiskerprops = dict(linestyle = '--', linewidth = 2));
    plt.title('Average Max D_up Peaks')
    plt.savefig(r'./Module6Graphs/AvgMaxDUPpks.svg')
    plt.figure(figsize=(4,6)) 
    sns.boxplot(data =list3,showfliers = False,showcaps = True, width = 0.6, medianprops = {'color': 'red'}, showmeans = True, meanprops = {"markerfacecolor": "red"}, color = 'blue', whiskerprops = dict(linestyle = '--', linewidth = 2));
    plt.title('Average ratios maxDupPeak/maxPromPeak')
    plt.savefig(r'./Module6Graphs/AvgRatiomaxDupPeak-maxPromPeak.svg')
    plt.figure(figsize=(4,6)) 
    sns.boxplot(data = list5, showfliers = False,showcaps = True, width = 0.6, medianprops = {'color': 'red'}, showmeans = True, meanprops = {"markerfacecolor": "red"}, color = 'blue', whiskerprops = dict(linestyle = '--', linewidth = 2));

    plt.title('SumDisUpDepths/SumPromDepths')
    plt.savefig(r'./Module6Graphs/SumDisUpDepths-SumPromDepths.svg')
 




cutoff(fi,cut)
log2FC()
pvalueSort()
binbar(bin)
peaksGenes(annot, bin)
transcriptL(quant, bin)
GeneL(annot, bin)
#gElement(annot, bin)
TSS(annot,bin)
lastP(annot,bin)
flavours(bin)
corre(bin)
