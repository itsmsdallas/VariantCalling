### Prerequisites: For this homework you will need to install bcftools
# !apt install bcftools

#Manipulating VCF files.

# Below there are the links to the chromosomal variant calls of a very, very famous person. Let's pry together!

# **DO NOT CLICK DOWNLOAD OR RUN**

# the files are very big and we just need to play with the data a little bit

# #### Here are the links

# !gdown 1IKSdvcpxJsAfZmKMxJUGs1tM4sPrvkJg     ### CHR Y
# !gdown 1nt-0mhYZzLeKc16bVnL2QWnnvIGbtiVt     ### CHR X
# !gdown 1K7EiX_3zE6vRFvKxhGfI4Q0MQGjAAXFC     ### CHR 21
# !gdown 1l9gI3onjyaxW4bO-6HxJEOb1XMwiB6DP     ### CHR 20
# !gdown 1Qa0ucdGfS4M9WrtecyVWxxhBRKIL6_Fo     ### CHR 19
# !gdown 1rAlsfFqlXx2283S-47L7xTDIvBr8ZH9B     ### CHR 18
# !gdown 1EttqAojXgjJb5xpR8PNrGQfhzrIB7kh5     ### CHR 17
# !gdown 1eCgqHeo1G-CkB5gG_7o1mRk8alQ-MNUM     ### CHR 16
# !gdown 1snKgxRigdqwF7uhHHYlhr_cB5zQkKDt-     ### CHR 15
# !gdown 1ZV8O9SwwAIgvx-GU42hkr69ZvyzfEs7Q     ### CHR 14
# !gdown 1Y7YXnw5lFm6sTwAtC7ttJpQyP6UuMzlQ     ### CHR 13
# !gdown 1DnFYEPwsoythazhU5fRQUWEdTIZ2S4T-     ### CHR 12
# !gdown 1b0vb63fQrUbm1tMKIwmYd-RLTCg12AK2     ### CHR 11
# !gdown 16Zd8BKA3q2E70DjrhPU49VB-KtqdcOPD     ### CHR 10
# !gdown 15WqxXHn4xBQGG42i1Ekj758dg3eH9gsC     ### CHR 9
# !gdown 1-TQIFClZ-sK9UxB8SVVE3_8NorV0GusB     ### CHR 8
# !gdown 1y3Bzaa-Rl6Ecx8s8YuBkOqFUM-svHWN0     ### CHR 7
# !gdown 1aF8kB6kmlAjEjaYkTahvuewrGe7lkmgG     ### CHR 6
# !gdown 1vTEUeUA3d4CuLIlGd_jEwfoXrfovQcAG     ### CHR 5
# !gdown 1y8NEsqFjan1OPFVVJYVrSVQVgOiTE3JQ     ### CHR 4
# !gdown 1XfDMzhacP0HoqbLbSkgO9LXJ2N5vtaLk     ### CHR 3
# !gdown 1xSyEMbxaJCw3BBxO00CYyLdpo_-yAjnZ     ### CHR 2
# !gdown 1bwEteqU3CBkXAp9Wkm9kdN_A9Fyd79uB     ### CHR 1


#### Pick up any chromosome that you wanna take a look at, copy the line of code to the box below and run

# HINT: Smaller chromosomes will run faster

### Download the data and unzip the data.
# !gdown 1snKgxRigdqwF7uhHHYlhr_cB5zQkKDt-     ### CHR 15



# from google.colab import drive
# drive.mount('/content/drive')

# %cd /content/drive/MyDrive/Bioinformatics in Python

# !gzip -df chr15.vcf.gz

# !more chr15.vcf

### THE VCFS IN THIS HOMEWORK HAVE ALL SITES, EVEN HOMOZYGOUS REFERENCE (GENOTYPE 0). WE WANNA LOOK AT VARIANT SITES ONLY SO WE RUN THIS:
### EXAMPLE: !bcftools view -v snps,indels chr8.vcf > Variable-SitesOnly.vcf

bcftools view -v snps,indels chr15.vcf > Variable-SitesOnly.vcf

# Using the python libraries that we used in class cyvcf2 parse out the vcf and write the output into a temporary file.

# HINT: Remember first to install cyvcf2

### Write your command to install cyvcf2
# !pip install cyvcf2

# Now, use the same code we have used before to parse the vcf file and save the output chromosome.tsv
# 
# We want to only output, CHROM, POS, REF, ALT, QUAL, and genotype. + 10 points

from cyvcf2 import VCF

vcf_reader = VCF('Variable-SitesOnly.vcf')

# Open a file in write mode to store the selected columns
with open('output_file.tsv', 'w') as file:
    # Loop through variants in the VCF file
    for variant in vcf_reader:
      genotype = variant.gt_types[0]
      data_string = f"{variant.CHROM}\t{variant.POS}\t{variant.REF}\t{','.join(variant.ALT)}\t{variant.QUAL}\t{genotype}\n"
      file.write(data_string)

!more output_file.tsv

# Import this file into a dataframe, call it anything you want. Write your code below, remember the output has fewer columns than that it did in class.



# HINT: `columns = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'GENOTYPE']`

### Write your code to export this into a dataframe + 10 points

import pandas as pd

columns = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'GENOTYPE']

df = pd.read_csv("output_file.tsv", sep = "\t", header=None, names= columns )
df

#Let's start exploring the data.

# Count the number of ocurrances of each genotype in the chromosome that you chosee to work with.
# HINT: Call the column `'GENOTYPE'` and apply the function `*.value_counts()*`

### Write your code here, execute it before saving it/sending it. + 4 points
df['GENOTYPE'].value_counts()

# Find the maximum, minimum and mean quality scores of the variant calls. Write them in independent code cells, + 2 point each.

#### Show the max. + 2 points
df['QUAL'].max()

#### Show the min + 2 points
df['QUAL'].min()

#### Show the mean + 2 points
df['QUAL'].mean()

# Subset the dataframe (name it df_subset). to include only positions in between 100000 and 1000000. Print the dataframe. + 5 points

#### Subset the dataframe to include positions in between 100000 and 1000000 + 3 points
subset_df = df[(df['POS'] > 100000) & (df['POS'] < 1000000)]
subset_df

##THIS RANGE IS TOO SMALL SEE THE NEXT CELL ##

filtered_df = df[(df['POS'] > 40000000) & (df['POS'] < 100000000)]
filtered_df



df['POS'].min()

df['POS'].max()



# Identify the number of mutations where the reference genome has a G or a T.
# HINT: Subset the reference genome using REGEX!! + 5 points

#### Write your code here to subset the original DF to show sites where the reference has a G or a T   (+ 3 points)
nuc = ['G', 'T']
G_or_T_df = df[df['REF'].isin(nuc)]
G_or_T_df

G->T mutations are very common. Find how many G->T Mutations happened in this genome.

# NOTE: Use the original DATAframe

# HINT: There are many ways to do it. I want you to do it the following way.
# First subset the original dataframe to positions where the reference has a G. Then subset this temporary df to show positions where the alt shows a T. Print this last DF. + 15 points

### Show your code here + 15 points
ref_has_G = df[df['REF'] == 'G' ]

alt_has_T = ref_has_G[ref_has_G['ALT'] == 'T']

alt_has_T

# Go on https://www.snpedia.com/, find a trait of interest, for example https://www.snpedia.com/index.php/Rs2472297 which is associated with coffee consumption. Does this person have this trait predictive snp? What can you tell me about it. NOTE. Find your own SNP! (This might take some time as you need to make sure you find them in your genome)

#### Work your code and results here + 4 points
## HERES THE SNPid rs182018947, this SNP is at position 48767448, it is associated with 	Primary autosomal recessive microcephaly
## HOMOZYGOSITY FOR C OR G at THIS POSITION IS ASSOCIATED WITH PATHOGENICITY


filtered_df[filtered_df['POS'] == '49059645']

# this person does not have this predictive snp (which makes sense)