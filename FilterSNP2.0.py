import pandas as pd
import numpy as np

def filter_uninformative_snp(file):
   
    #read file and convert to dataframe
    #file = "(bench)CA346865.hmp.txt"
    df = pd.read_csv(file, sep="\t", low_memory=False)
    
    #separate taxa and genotype from other informations
    #rs,alleles,chrom,pos,strand,assembly,center,protLSID,assayLSID,panelLSID,QCcode
    other_informations = {'rs': df['rs'],
                          'alleles': df['alleles'],
                          'chrom': df['chrom'],
                          "pos": df['pos'],
                          "strand": df['strand'],
                          "assembly": df['assembly'],
                          "center": df['center'],
                          "protLSID": df['protLSID'],
                          "assayLSID": df['assayLSID'],
                          "panelLSID": df['panelLSID'],
                          "QCcode": df['QCcode']}
    
    #drop columns other than genotype data
    eliminate_column = ['rs','alleles','chrom','pos','strand','assembly','center','protLSID','assayLSID','panelLSID','QCcode']
    df.drop(columns=eliminate_column, inplace=True)
      
    #add column SNP ID
    df.insert(0, 'SNP', other_informations['rs'], allow_duplicates=False)
    
    #identify number of genotype within each SNP
    report=df.nunique(axis=1, dropna=True)

    #identify uninformative SNP
    homogenous_snp = []
    for i in range(len(df.index)):
        if report[i]==1:
            homogenous_snp.append(i)
                   
    #delete uninformative SNP
    df.drop(labels = homogenous_snp, axis=0, inplace=True)
    
    #return dataframe
    return df, other_informations

def filter_snp_by_accession(df, other_informations):
   
    #get accession number
    taxa = list(df.columns)
    taxa.remove('SNP')
    accession = []
    for i in range(len(taxa)):
        sample = taxa[i]
        Ca = sample.rsplit("_",1)
        accession.append(Ca[0])
    
    #Reformat the dataset
    df2 = df.T 
    df2.columns = other_informations['rs']
    df2.drop('SNP', inplace=True)
    df2['Accession'] = accession

    #remove SNP based on accession's genotype similarity
    df2 = df2.groupby('Accession').nunique()
    snp_to_drop = df2.columns[(df2 != 1).any()]
    snp_to_drop = list(snp_to_drop)
    
    #reformat the dataset
    df.drop(columns='SNP', inplace=True)
    df_info = pd.DataFrame(other_informations)
    dataset = pd.concat([df_info, df], axis=1, sort=False)
    dataset = dataset.set_index('rs')
    dataset.drop(labels=snp_to_drop, inplace=True, axis=0)
    dataset.reset_index(inplace=True)
    
    return dataset

def main():
    file = input("enter file name or file path: ")
    df, other_informations = filter_uninformative_snp(file)
    dataset = filter_snp_by_accession(df, other_informations)
    dataset.to_csv("(filtered)CA346865.hmp.txt", sep='\t',na_rep='NA', index=False)
    
if __name__ == "__main__":
    main()
