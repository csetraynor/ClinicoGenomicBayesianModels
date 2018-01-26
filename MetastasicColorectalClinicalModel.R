#Download Colorectal Metastasic Cancer Dataset

#get data from cbioportal 

require(cgdsr)
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
test(mycgds)

listofcancerstudies = getCancerStudies(mycgds)
#View(listofcancerstudies)

#get study id
id_sutdy = getCancerStudies(mycgds)[155,1] #targeted colorectal cancer
#id_sutdy_history = getCancerStudies(mycgds)[56,1] Get another study for power prior?

#get case list
case_list = getCaseLists(mycgds, id_sutdy)[1,1] #All tumours include also NGS, and CNA
#case_list_ngs = getCaseLists(mycgds, id_sutdy)[2,2] #substract only NGS?

#get avialable genetic profiles
geneticprofile = getGeneticProfiles(mycgds, id_sutdy)

#get data slices for a specified list of genes, genetic profile and case list
#gene_data_slices = getProfileData(mycgds, c('PIK3CA', 'TP53'), geneticprofile, case_list)

#Get Clinical Data for the case list
clinical_data <-  getClinicalData(mycgds, case_list)
#documentation
#help('cgdsr')

