import pandas as pd
import sys
import os

### Loading files, the input is the output from interproscan analysis

file=pd.read_csv(sys.argv[1], sep='\t', header=None)
filtrato=file.drop_duplicates(subset=[0,11])
filtrato = filtrato[filtrato.loc[:,11].str.contains("-") == False]
gene=filtrato.loc[:,[0,11]]
# sys.argv[3] argument in the script IPR333,IPRetc....= personalize domains from the user. Important: the domains have to be equal to interproscan 
# domains 
input_domains=sys.argv[3].split(',')
gene_to_find=set(input_domains)

name=list(gene.loc[:,0])
uniprot=list(gene.loc[:,11])

if os.path.exists(sys.argv[2]):
    os.remove(sys.argv[2])

# Check if the domain of the orthoDB proteins are equal to user input and take only the proteins with the same domains)
# sys.argv[2]= output .txt file
with open(sys.argv[2], 'x') as f:
    count=[]
    for x in range(0,len(name)):
        if x==0 and uniprot[0] in gene_to_find:
            count.append(uniprot[0])
        # if the next protein is equal to the previous one take the uniprot domain
        elif name[x] == name[x-1]:
            count.append(uniprot[x])
        else:
        # the proteins are different so compare the set that will be generated with the set including the official protein domains
        #The order of the original list items is not important, because the == operator returns true when each set contains identical items in any order.
            count1=set(count)
            if count1 == gene_to_find:
                f.write(name[x-1])
                f.write('\n')
            else:
                pass
            # reset the count and check the next protein domain after that
            count=[]
            uni=uniprot[x]
            if uni in gene_to_find:
                count.append(uni)
            else:
                count=['1','2','3','4','5']
           
        

f.close()


