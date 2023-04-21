##This script will convert the format filed to info fields


##Extract information for the info fields
vcffile="3929_MungeSumstats_Out.vcf.gz"
os.system('''gatk VariantsToTable -F CHROM -F POS -F REF -F ALT -ASGF INFO -ASGF FRQ -ASGF SNP  -V '''+vcffile + " -O MasterTable.tsv ")


##Read file
sparq=SparkSession.Builder().appName("Dataframe").getOrCreate()
commondf=sparq.read.csv("MasterTable.tsv",header=True,inferSchema=True,sep="\t")
new_cols=['CHROM', 'POS', 'REF', 'ALT', 'INFO', 'FRQ', 'SNP']
commondf = commondf.toDF(*new_cols)
commondf = commondf.withColumn('SNP', F.regexp_replace('SNP', r'_', ':'))

for col_name in ['INFO', 'FRQ']:
    commondf = commondf.withColumn(col_name, col(col_name).cast(StringType()))


commondf=commondf.withColumn("INFO", concat(lit("INFO="), "INFO"))
commondf=commondf.withColumn("FRQ", concat(lit("FRQ="), "FRQ"))
commondf=commondf.withColumn("SNP", concat(lit("SNP="), "SNP"))


commondf = commondf.withColumn("INFO", concat(col("SNP"), lit(";"), col("INFO"), lit("-"), col("INFO")))
commondf=commondf.select(['CHROM', 'POS', 'REF', 'ALT','INFO'])
commondf=commondf.withColumn("QUAL", lit("."))  
commondf=commondf.withColumn("FILTER", lit("."))
commondf=commondf.withColumn("ID", lit("."))
commondf=commondf.select(['CHROM','POS','ID','REF','ALT','QUAL',"FILTER",'INFO'])




os.system(f"zgrep \"##\"  {vcffile} >vcf_header.txt")
header='''##INFO=<ID=INFO,Type=Float,Description="INFO Score">
##INFO=<ID=FRQ,Type=Float,Description="FRQ of Alt allele">
##INFO=<ID=SNP,Type=String,Description="SNP ID">
'''


f = open("vcf_header.txt", "a")
f.write(header)
f.close()


commondf.coalesce(1).write.option("sep", "\t").csv('result.csv',header=True)
os.system(f"mv  result.csv/*csv Annotation_Input.tsv")
os.system("sed -i 's/CHROM/#CHROM/' Annotation_Input.tsv")


os.system("cat vcf_header.txt Annotation_Input.tsv > Annotation_Input.vcf")


os.system("bgzip Annotation_Input.vcf")
os.system("tabix -p vcf Annotation_Input.vcf.gz")

bcftools annotate -a Annotation_Input.vcf.gz -c INFO/SNP,INFO/INFO,INFO/FRQ 3931_MungeSumstats_Out.vcf.gz | bgzip -c > Annotated_output.vcf.gz

