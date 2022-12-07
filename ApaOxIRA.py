"""
                              ApaOxIMOD
                    Apatite Oxygen Isotopes Simulation
                    
Monte-Carlo simulation of the alteration of the original oxygen isotope composition of biogenic apatites

                    JP Flandrois* and C Lecuyer**
                            2021  mk02

* Author and maintainer
** Equations and solution

How it works :
The risk of not identifying a diagenesis process by using the Normality criteria (-r) is estimated by an iterative Monte-Carlo simulation. Monte Carlo simulation is based on the creation of an array of A/W values taken from a uniform distribution between given threshold values of δ18O. This represents the diversity of the deposit situation and physical state of the biological apatite. The array is composed of subarrays with length corresponding to the number n of samples in a sample set. Typically 100,000 sets of length n are simulated. The second phase uses equation (4) to compute the δ18O in biological apatite δ18OAf at the equilibrium on each item of the whole array given the temperature T, initial δ18O in biological apatite (δ18OAi) and the δ18O in water δ18OWi. (These two steps may be repeated once, using the array issued from the previous computation (the array of δ18OAf values) as initial value of δ18OA with different values of A/W, δ18OWi and T.) The resulting δ18OAf array is then optionnaly used to simulate the analytical uncertainty by using Normal law random sampling whose mean being the individal δ18OAf and given standard deviation (default 0.05). Finaly each subarray representing a sample set of lenth n of δ18OAf values is submitted to the normality tests. Shapiro-Wilk and the Anderson-Darling are combined. Both have demonstrated a high potential in normality detection, but they are using different approaches and have different sensitivity to skew and tails. Normality is rejected if anyone test reject it at the given risk α. The number of sets in the array were normality cannot be rejected is computed. The result is expressed as a percentage of the number of sample sets. As Normality is unexpected from the construction but may arise by chance (and the shorter the sample, the higher the chance ), this is the MC estimated risk of the sample strategy with sets of length n. Lastly, a bootstrap (BT) estimator of the mean and standard deviation is computed. The whole simulation is done iteratively with increasing values for the length of the sets with given limits for the exploration and the parameters described in a file using the yaml format.

Statistical (-s) analysis of experimental data: Basic statistics (mean, variance, skew and kurtosis) are computed for a list of data (δ18O, one column only). The normality tests (Shapiro-Wilk, Anderson-Darling) are used and a joined normality test is computed that enable us to detect or reject normality. The interpretation criteria being that of the MC simulation parameters, the results may be compared to the MC simulation. MC simulation performed for the given number of samples and basic parameters give us the _risk_ of getting a normal distributed population by chance after diagenesis. Additionaly an histogram, the adjusted Normal law and the quantile-quantile plot is also returned.
"""
import argparse, pickle, os
from multiprocessing import Pool, cpu_count
#import multiprocessing
import yaml #yaml format for the parameters
#import numpy as np
from numpy import array as npArray
import matplotlib.pyplot as plt
from scipy.stats import describe,shapiro,anderson,skew,skewtest,kurtosistest,probplot
from scipy.stats import norm as Norm
from numpy import mean as npMean
from numpy import std as npStd
from numpy import random as npRandom
from numpy import linspace
from numpy import histogram_bin_edges
#import O18parameters as param
import time




def readParameters(parameterYAML):
    """
    The 8 parameters are read from a YAML file O18parameters.yaml where they are stored as a python dictionaries
    And all the possibilities of MC computing are described
    This is by far the coolest way
    Reading is done at the beginning
    The structure of the dictionary parameterdictionary is shown in the function basicParameters()
    """
    #global parameterdictionary
    yaml_file = open(parameterYAML, 'r')
    parameterdictionary = yaml.load(yaml_file, Loader=yaml.SafeLoader) #this is the safe way tp prevent exploits
    
    return parameterdictionary
    
def extractParameters(parameterdictionary,phase):
    """
    The 8 parameters are read from a YAML file like O18parameters.yaml where they are stored as a python dictionaries
    And all the possibilities of MC computing are described
    This is by far the coolest way
    We return the necessary datas
    """
    #basic cas of one diagenesis event
    T=parameterdictionary["parameterdictionary0"]["T"]
    d18OWi=parameterdictionary["parameterdictionary0"]["d18OWi"]
    d18OAi=parameterdictionary["parameterdictionary0"]["d18OAi"]
    sampleNB=parameterdictionary["parameterdictionary0"]["sampleNB"]
    alpha=parameterdictionary["parameterdictionary0"]["alpha"]
    arraysize=parameterdictionary["parameterdictionary0"]["arraysize"]
    parameterdictionary["parameterdictionary0"]["arraysize"]=(arraysize//5000)*5000 #we need a number that can be divided in parts of 5000
    normalizedArraysize=arraysize//5000 # the number of arrays is divided by batches of 5000 to be used more rapidly by pool
    parameterdictionary["parameterdictionary0"]["arraysize"]=normalizedArraysize
    WAlow=parameterdictionary["parameterdictionary0"]["WAlow"]
    WAhigh=parameterdictionary["parameterdictionary0"]["WAhigh"]
    #two events
    SecondDiagenesis=parameterdictionary["SecondDiagenesis"] #if True: chaining two diagenesis events
    T1=parameterdictionary["parameterdictionary1"]["T1"]
    d18OWi1=parameterdictionary["parameterdictionary1"]["d18OWi1"]
    WAlow1=parameterdictionary["parameterdictionary1"]["WAlow1"]
    WAhigh1=parameterdictionary["parameterdictionary1"]["WAhigh1"]
    #step increase
    steppingOn=parameterdictionary["stepwiseIncreaseNbSamples"] #if True: step increase of stepNumbers
    stepNumbers=parameterdictionary["stepping"]["stepNumbers"]
    stepAmplitude=parameterdictionary["stepping"]["stepAmplitude"]
    #lab process
    analyticalProcessSimulation=parameterdictionary["analyticalProcessSimulation"]
    sigmaLab=parameterdictionary["experimental"]["sigma"]
    
    if phase == "FormatedSummary": #This is the generator of the bottom 'summary' block for the MC simulation template
        Answer="\n\n############################# PARAMETERS OF THE  MC SIMULATION ##################################\n\n* General conditions *\n"+"Number of simulations : "+str(arraysize)+"\n"+"Admittted risk for p value :"+str(alpha)+"\n"+"\n* Primary MC simulation *"+"Temperature :"+str(T)+"\n"+"Initial δ18O Water: "+str(d18OWi)+"\n"+"Initial δ18O Apatite: "+str(d18OAi)+"\n"+"A/W ratio range: "+str(WAlow)+"-"+str(WAhigh)+"\n"+"Number of samples: "+str(sampleNB)+"\n"
        Answer+=("\nSimulation of the lab process, sigma: "+str(sigmaLab)+"\n" if analyticalProcessSimulation else "NO Simulation of the lab process\n")
        if parameterdictionary["SecondDiagenesis"]:
            Answer+="\n* Second environment MC simulation *\n"+"Temperature: "+str(T1)+"\n"+"Initial δ18O Water: "+str(d18OWi1)+"\n"+"A/W ratio range: "+str(WAlow1)+"-"+str(WAhigh1)+"\n"
        if parameterdictionary["stepwiseIncreaseNbSamples"]:
            Answer+="\n* Increase of the number of samples by steps *\n"+"number of steps: "+str(stepNumbers)+"\n"+"level of the step: "+str(stepAmplitude)+"\n" #
        return Answer
    if phase == "DiagenesisTwo" :
        return (T1,d18OWi1,WAlow1,WAhigh1,analyticalProcessSimulation,sigmaLab) #this "DiagenesisTwo" information is from inside the program, not from the yaml file
    elif phase != "iteratiV" : #no iteration over the number of samples
        if not parameterdictionary["SecondDiagenesis"] :#True only when a second diagenesis phase
            return (False,steppingOn,T,d18OWi,d18OAi,sampleNB,alpha,arraysize,WAlow,WAhigh,analyticalProcessSimulation,sigmaLab) #no second phase
        else: return (True,steppingOn,T,d18OWi,d18OAi,sampleNB,alpha,arraysize,WAlow,WAhigh,analyticalProcessSimulation,sigmaLab) #seconf phase
    else :
        return (stepNumbers,stepAmplitude) #parameters for iteration
        

def arraysFiller(T,d18OWi,d18OAi,sampleNB,alpha,arraysize,WAlow,WAhigh,analyticalProcessSimulation,sigmaLab):
    """
    An array is constructed (WAarray) by sampling in an uniform repartition law with boundaries defined by WAlow,WAhigh
    Then a new array (d18OAf_array) is built by applying equation of
    At the end the normality tests are done and the error ratio returned
    When a second diagenesis is required this is done twice but analytical process is then suppressed during the first simulation
    """
    CL=(117.4-T)/4.5 #equation 3 Lécuyer et al. (2013) - CL is for ChristopheLécuyer
    WAarray=npRandom.uniform(WAlow,WAhigh, (arraysize,sampleNB))
    d18OAf_array=computeEquation(CL,d18OWi,d18OAi,WAarray)
    if analyticalProcessSimulation: #in the case of dual diagenesis the first simulation is done alone, this is also the case of one diagenesis if the analyticalProcessSimulation is not required
        d18OAf_array=npRandom.normal(d18OAf_array,sigmaLab,(arraysize,sampleNB))
    return d18OAf_array
    

def computeEquation(CL,d18OWi,d18OAi,WA):
    """
    This is the computation corresponding to equation (4) in our paper.
    It has been separated as a function to clarified the program.
    As the p value level comparison is not relevant, only the putative normality is returned as a boolean
    """
    Eq4result=(WA*(d18OWi+CL)+d18OAi)/(WA+1.0)
    return Eq4result
    
def AndersonDarlingSimple(cell,alpha): # a renommer
    """
    The Anderson-Darling test is based on the distance between the empirical distribution function and the hypothetic Normal distribution.
    AD tests is one in the family of the quadratic class of the EDF statistic as it used squared differences.
    Note that the computing time is increased for the low numbers and high sparcity of the samples.
    As the p value level comparison is not relevant, only the putative normality is returned as a boolean
    """
    A,B,C = anderson(cell)
    return False if A > B[2] else True
    
def ShapiroWilkTestSimple(cell,alpha):
    """The Shapiro-Wilk test is based on the test statistic called W that may be interpreted as the coefficient of determination between quantiles from the Normal law and empirical quantiles obtained from the data.
    Note that the computing time is increased for the low numbers and high sparcity of the samples.
    """
    stat, psh = shapiro(cell)
    #print('SH Statistics=%.3f, p=%.3f' % (stat, psh))
    return  False if psh < alpha else True
#
def normalityTestDecisionLevel(d18OAf_array,dimension,alpha,VerboseOutput) :
    """
    Normality is unexpected from the construction but may arise by chance (and the shorter the sample, the higher the chance exists).
    The percentage of samples following the normal law by chance is also the confidence level of the sampling strategy.
    ---
    This function computes two normality tests, the Shapiro-Wilk and the Anderson-Darling tests.
    Departure from normality due to skewness or kurtosis is more taken into account with the Shapiro-Wilk test and Aderson-Darling test is sensitive to tail.
    The Shapiro-Wilk test is based on the test statistic called W that may be interpreted as the coefficient of determination between quantiles from the Normal law and empirical quantiles obtained from the data.
    The Anderson-Darling test is based on the distance between the empirical distribution function and the hypothetic Normal distribution.
    AD tests is one in the family of the quadratic class of the EDF statistic as it used squared differences.
    Both tests have demonstrated a high potential in normality detection, but as they are using different approaches it seems interesting to use both.
    To improve the speed, the Anderson-Darling test is done only if the Shapiro-Wilk test does not evoques normality.
    Note that the computing time is increased for the low numbers and high sparcity of the samples.
    ---
    The ratio of positive answers is returned.
    This ratio is expressed as a percentage.
    
    """
    NormalityCount=0
    times=[]
    i=0
    while i<len(d18OAf_array):
        cell = d18OAf_array[i]
        if ShapiroWilkTestSimple(cell,alpha) and AndersonDarlingSimple(cell,alpha):
                NormalityCount +=1
        else: pass
        i+=1
    try:
        return len(d18OAf_array),len(cell),NormalityCount # !!! results are returned as numerical value not ratio
    except:
        print ('---error---',d18OAf_array)
        stop('normalityTestDecisionLevel'+' try-except error detected') #brutal end essentially for testing

    
def MCsimulationOnePhase(T,d18OWi,d18OAi,sampleNB,alpha,arraysize,WAlow,WAhigh,VerboseOutput,analyticalProcessSimulation,sigmaLab):
    """
    the array d18OAf_array is generated by the function arraysFiller, see explainations for this function
    then the function ComputeerrorRatio launches normalityTestDecisionLevel
    this is the sheduling function doing nothing except transmitting from one side to another
    """
    d18OAf_array=arraysFiller(T,d18OWi,d18OAi,sampleNB,alpha,arraysize,WAlow,WAhigh,analyticalProcessSimulation,sigmaLab)
    flatArray=d18OAf_array.flatten()
    lend18OAf_array,lencell,NormalityCount=normalityTestDecisionLevel(d18OAf_array,arraysize,alpha,True)# verbose false
    listResultsNormalityCount=[]
    listResultsNormalityCount.append([sampleNB,NormalityCount])
#        print("\n# Results for the Quasi-bottstrap Estimated Statistics\n")
#        print ("Estimated mean ",npMean(listResuPercent),"Estimated std",npStd(listResuPercent), end='')
    return d18OAf_array,listResultsNormalityCount # this is the last array, so we may run the simulation again with other parameters

def automatic(parameterdictionary,cpuCount):
    """
    the function is basically sending the job to do to the threads via the Pool+map process
    The first phase is used to divide the required array size by slides of 5000 that will be used by the bootstrap
    Of course we may increase the size of the slides but it is then difficult to do a decent bootstrap under 250 000
    We adjust also the computing charge on all the threads by increasing the number of slides
    The idea is to have (statistically) the same computing charge on all the threads and maximize the number of slides
    Then the process is launched on all the threads
    """
    #adjusting the Nb of threads to obtain a multiple of 5000 and the number of threads to be as close of the users demand as possible (and over if needed)
    INI=parameterdictionary["parameterdictionary0"]["arraysize"]
    if parameterdictionary["parameterdictionary0"]["arraysize"] < 100000: parameterdictionary["parameterdictionary0"]["arraysize"]=100000
    if parameterdictionary["parameterdictionary0"]["arraysize"] > 5000000: parameterdictionary["parameterdictionary0"]["arraysize"]=5000000
    
    NbSlices=parameterdictionary["parameterdictionary0"]["arraysize"]//5000
    
    if parameterdictionary["parameterdictionary0"]["arraysize"]*NbSlices < INI:
        NbSlices+=1 #to solve the problem of illegal arraysize (not multiple of 5000)
    parameterdictionary["parameterdictionary0"]["arraysize"]=5000 #the threads will always work with 5000 positions
    if NbSlices%cpuCount !=0 and NbSlices > cpuCount: #but we need to work with a multiple of the cpuCount
                                                    #this is not mandatory but the objective is that all the threads are fully working up to the end
                                                    #suppressing gives a gain of anly few seconds and I prefer to fully use the threads until completion
        i=1
        while NbSlices%cpuCount !=0:
            NbSlices+=1
            i+=1
    print ("Desired size",INI,"adjusted size",parameterdictionary["parameterdictionary0"]["arraysize"]*NbSlices)
    print ("WORKING with the dimension ",parameterdictionary["parameterdictionary0"]["arraysize"]*NbSlices, "on ",cpuCount," Threads"," and ",NbSlices,"slices")
    
    pickle.dump(parameterdictionary, open("tempo.pkl", "wb")) #parameterdictionary is send to the core function used by thethreads by this way
    #unfortunately sending parameterdictionary to each process was buggy so this very odd and dirty solution
    #Problem: dont launch the program twice as "tempo.pkl" will be shared and you may mix the parameters...
    #note that we cannot use a global definition ! Some clever solution may exist but it works fine.
    
    with Pool(cpuCount) as p: #sending to the n threads (i.e cpuCount)
        results=p.map(coreFunctionForPool, list(range(0,NbSlices))) #Slices are arrays of 5000
        return results,parameterdictionary["parameterdictionary0"]["arraysize"]*NbSlices,cpuCount,NbSlices
    #all the work is done, results contain the results of the n threads that we will analyse
    
def coreFunctionForPool(a):
    """
    This is what is done by the threads individually
    Essentially a scheduler
    Note the two situation : one or two diagenesis
    The result = the Normality tests results and NOT the arrays !!!
    """
    parameterdictionary=pickle.load(open("tempo.pkl", "rb")) #read the parameters
    Chained_MC,steppingOn,T,d18OWi,d18OAi,sampleNB,alpha,arraysize,WAlow,WAhigh,analyticalProcessSimulation,sigmaLab = extractParameters(parameterdictionary,"DiagenesisOne")
    #we send "DiagenesisOne" that is a default value in extractParameters function i.e it is not explicitely described there :)
    if steppingOn:
        stepNumbers,stepAmplitude=extractParameters(parameterdictionary,"iteratiV")
    else:
        stepNumbers=0
        stepAmplitude=1 #mandatory to stop the loops
    cpt=0
    #the first MC is always done
    maxIter=sampleNB+stepNumbers
    
    resuresu=[]
    while cpt <= stepNumbers:#aka sampleNB <= maxIter:
        if Chained_MC: #two diagenesis this is the information send inside the program
            T1,d18OWi1,WAlow1,WAhigh1,analyticalProcessSimulation,sigmaLab=extractParameters(parameterdictionary,"DiagenesisTwo")
            # first pass
            d18OAf_array,resu=MCsimulationOnePhase(T,d18OWi,d18OAi,sampleNB,alpha,arraysize,WAlow,WAhigh,False,False,sigmaLab) #trick, we put analyticalProcessSimulation as False to prevent doing the lab process in the past :)
            # d18OAf_array is the resulting array after the first diagenesis and replaces INITDELTA
            #second pass ==> the output is run again with the seconds parameters
            Rd18OAi,resu =MCsimulationOnePhase(T1,d18OWi1,d18OAf_array,sampleNB,alpha,arraysize,WAlow1,WAhigh1,True,analyticalProcessSimulation,sigmaLab) #here the laboratory process is in accordance to the pilot.
            
            resuresu.append(resu)
            cpt+=1
            sampleNB+=stepAmplitude
            
        else:   # case of increase of sampleNB by stepAmplitude
                # with sampleNB <= maxIter
                # verified by another program :)
            Rd18OAi,resu=MCsimulationOnePhase(T,d18OWi,d18OAi,sampleNB,alpha,arraysize,WAlow,WAhigh,True,analyticalProcessSimulation,sigmaLab)
            cpt+=1
            sampleNB+=stepAmplitude
            resuresu.append(resu)
    return resuresu
    
def MC_run(proto,fileIs):
    """
    The entrance of the program MS simulation; dispatcher of the tasks for MC simulation
    It read the parameter and decide what to do when we need a step by step increase of the number of samples and when to launch the second diagenesis simulation
    It is also used for the publication of the results
    """
    tjob0 = time.time()
    #reading the schedules
    #global UserResu #to be used elsewere to write the result, easier to declare a global variable
    #with open(fileIs, "w") as UserResu: pass #emptying, ready to use
    if args.Normality : #MC simulation in the case of users data
        parameterdictionary=basicParameters() #a minimal set for somewhat situation
        parameterdictionary["parameterdictionary0"]["sampleNB"]=proto #change for the true nb of sample of the users data
        
    
    else: parameterdictionary=readParameters(proto) #reading for all but without global definition
    
    if parameterdictionary["threadsCount"]=="auto": #simple but efficient due to the permissive shape of the thread-count - speed curve
                                                    #70% of the threads allowed and a maximum of 16 determined from our tests (see below)
                                                    #this could be an under-estimation so it may be changed
                                                    #ideally we should have tested on the user machine by running a test seet with various options
                                                    #but it may was longer than running the program even with an approximation of the optimum core number :)
        if int(cpu_count()*0.7) !=0 :
            cpuCount=min(int(cpu_count()*0.7),16)
        else : cpuCount=1 #
    else : cpuCount= int(parameterdictionary["threadsCount"]) #user fixed
    #We could have limited the number to 10-12 as on our machine, but this may be computer specific :
    #Speed is linked to the number of involved threads, it increases rapidly until 10 threads on our computer (20 threads),
    #then the gain is very poorly growing (a growing plateau) on the test machine. When more than 10 threads are used,
    #the speed is also varying marginaly due to the automated increase of the number of sets
    #that harmonize the number of simulation computed on each thread.
    #Thus working with a threadCount in the plateau region is OK at the expense of waste of energy with marginal gain.
    
    resultsGlobal,MCarray,cpuCount,NbSlices=automatic(parameterdictionary,cpuCount) #the parallelization
    #NOW THE OUTPUT OF THE RESULTS (basic presentation)
    tjob1 = time.time()
    print ("##### RESULTS ##### Time to end of MC : {:.2f} #####".format(float(tjob1-tjob0)))
    cpt=0
    dictionnaireStatsSteps={}
    for u in resultsGlobal:
        for i in u:
            #print (cpt,i[0])
            if i[0][0] in dictionnaireStatsSteps:
                dictionnaireStatsSteps[i[0][0]].append(i[0][1])
            else :
                dictionnaireStatsSteps[i[0][0]]=[i[0][1]]
        cpt+=1
    
    alpha=parameterdictionary["parameterdictionary0"]["alpha"]
    
    for i in dictionnaireStatsSteps:
        ratioList=[]
        for u in dictionnaireStatsSteps[i]:
            ratioList.append(100*(u/(MCarray/NbSlices)))
        cpt=0

        #print ("Number of samples: {} | Estimated BT mean: {:.3}% | Estimated BT std: {:.3}".format(i,float(npMean(ratioList)),float(npStd(ratioList)),width=5))
    if args.Risk or args.xRepeatedMC:
        with open(fileIs, "w") as UserResu: pass #emptying
    with open(fileIs, "a") as UserResu:
        UserResu.write("#####################################  ApaOxIS ###################################################\n")
        UserResu.write("          MC simulation of evolution of δ18O in apatite during diagenesis\n")
        UserResu.write("                       C. Lécuyer & JP. Flandrois 2021\n")
        UserResu.write("<reference of the paper>\n<reference of the program repository>\n")
        UserResu.write("\n###################################  RESULTS  ###################################################\n(you may copy the results to any place as csv)\n\n")
        UserResu.write("Nb_Samples"+";"+"Risk_Level %"+";"+"std"+"\n")
        print("\n\n##################################  MC SIMULATION ################################################\n")
        listNbSamples=[]
        listMeanBTestimated=[]
        listStdBTestimated=[]
        for i in dictionnaireStatsSteps:
            ratioList=[]
            for u in dictionnaireStatsSteps[i]:
                ratioList.append(100*(u/(MCarray/NbSlices)))
            cpt=0
            
            #arrayOf100BT=npRandom.choice(ratioList,100)
            rpt=0
            arrayOf100BTMeans=[]
            #listOfStd=[]
            while rpt < 100:
                arrayOf100BT=npRandom.choice(ratioList,int(len(ratioList)*0.75))
                #print (float(npMean(arrayOf100BT)),float(npStd(arrayOf100BT)))
                arrayOf100BTMeans.append(float(npMean(arrayOf100BT)))
                #listOfStd.append(float(npStd(arrayOf100BT)))
                rpt+=1
            EstMean=float(npMean(arrayOf100BTMeans))
            EstStd=float(npStd(arrayOf100BTMeans))
            print ("Number of samples: {} | Estimated BT mean: {:.3} %| Estimated BT std: {:.3}".format(i,EstMean,EstStd,width=5))
            listNbSamples.append(i)
            listMeanBTestimated.append(EstMean)
            listStdBTestimated.append([EstMean+2*EstStd,EstMean-2*EstStd])
            UserResu.write("{} ; {:.3} ; {:.3} \n".format(i,EstMean,EstStd))
            #UserResu.write("{} ; {:.3} ; {:.3} \n".format(i,EstMean,EstStd))
        if args.Normality: parameterdictionary=basicParameters()
        else: parameterdictionary=readParameters(proto) #read the original file to get the requirements of the user
        
        UserResu.write("\n"+extractParameters(parameterdictionary,"FormatedSummary")+"\n")
        UserResu.write("The adapted number of simulated conditions were: "+str(MCarray)+" and computing was done on "+str(cpuCount)+" threads in parallel\n")
        tjob2 = time.time()
        UserResu.write("Time to results : {:.2f} seconds".format(tjob2-tjob0))
        UserResu.write("\n#################################################################################################\n")
        print ("##### Time to results : {:.2f} #####".format(float(tjob2-tjob0)))
        print("\n\n############################# THE RESULTS ARE READY  #############################################\n")
        
        print("You will find it HERE ==> "+fileIs)
        print("                                                                         ")
        
        if len(listNbSamples)>1 and not args.xRepeatedMC:
            plt.plot(listNbSamples,listStdBTestimated,"_k")
            plt.title("ApaOxIS Risk assessment for a samples set")
            plt.ylabel('Estimated mean and 95% Confidence interval')
            plt.xlabel('Nb of samples')
            plt.plot(listNbSamples,listMeanBTestimated,".r")
            plt.savefig(fileIs+".graphRisk.pdf")
            print("and graph output is HERE ==> "+fileIs+".graphRisk.pdf")
            #plt.show()
        print("\n\n##################################  ApaOxIS END ##################################################\n")
    return MCarray,i,EstMean,EstStd


def basicParameters():
    """
    In the case of the statistics of normality for submitted data we use some basic parameters: a globally satisfying set of conditions
    the parameters are, as in the simulation, gatered in nested dictionaries in parameterdictionary
    Here we do not need to read the pilot file
    This is a the function linked to the users data and we do not complicate the main function used by the MC simulation
    """
    parameterdictionary={'parameterdictionary0': {'arraysize': 250000, 'sampleNB': 100, 'T': 15, 'd18OWi': -8.0, 'd18OAi': 20.0, 'alpha': 0.05, 'WAlow': 0.05, 'WAhigh': 0.95}, 'analyticalProcessSimulation' : True, 'experimental': {'sigma' : 0.05}, 'stepwiseIncreaseNbSamples': False, 'stepping': {'stepNumbers': 10, 'stepAmplitude': 1}, 'SecondDiagenesis': False, 'parameterdictionary1': {'T1': 40, 'd18OWi1': -8.0, 'WAlow1': 0.1, 'WAhigh1': 0.5},'analyticalProcessSimulation':True,'experimental':{'sigma':0.05},'threadsCount':'auto'}
    
    return parameterdictionary
    
def createFolder(directory):
    """
    From BIBIengine and TEALS (jp Flandrois 2020), extremely basic and I suppose that it can be found elsewhere
    """
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
        else:
            os.system('rm -r '+directory)
            os.makedirs(directory)
    except OSError:
        print(('Error: Creating directory. ' +  directory))

def multiMCstatistics(DirIs, fileout,tjob0):
    """
    the function that treat a 30 time repetition of the MC and computes statistics and built graphs
    This function has been written to built the graphs for the paper
    30 x runs to get a minimal statistic power
    """
    with open(fileout, "w") as UserResu: pass #emptying
    with open(fileout, "a") as UserResu:
        UserResu.write("#####################################  ApaOxIS ###################################################\n")
        
        UserResu.write("                    Normality tests according of bio-apatite                                      \n")
        UserResu.write("                       C. Lécuyer & JP. Flandrois 2021\n")
        UserResu.write("<reference of the paper>\n<reference of the program repository>\n")
        UserResu.write("\n#######################  STATISTICS ON A x30 REPLICATION ########################################\n")
        UserResu.write("Nb_Samples"+";"+"Risk_Level %"+";"+"std"+"\n")
        dictionnaireData={}
        dictionnaireStd={}
        for u in os.listdir(DirIs):#
            if u.find(".DS_Store") ==-1 :
                with open(os.path.join(DirIs,u),"r") as thefile:
                    THEFILE=thefile.read()
                    STHEFILE=THEFILE.split("std\n")[1].split("\n\n\n")[0].strip()
                    #print("-----")
                    #print (STHEFILE)
                    for j in STHEFILE.strip().split('\n'):
                        #print (j)
                        laligne=j.split(';')
                        if int(laligne[0].strip()) not in dictionnaireData.keys():
                            dictionnaireData[int(laligne[0].strip())]=[float(laligne[1].strip().replace(',','.').replace('%',''))]
                            #dictionnaireStd[int(laligne[0].strip())]=[float(laligne[2].strip().replace(',','.').replace('%',''))]
                        else:
                            dictionnaireData[int(laligne[0].strip())].append(float(laligne[1].strip().replace(',','.').replace('%','')))
                            #dictionnaireStd[int(laligne[0].strip())].append(float(laligne[2].strip().replace(',','.').replace('%','')))
        Echelle=list(dictionnaireData.keys())
        Data=[]
        Std=[]
        print ("Statistics for the replications")
        for m in dictionnaireData:
            meanGroup=npMean(dictionnaireData[m])
            stdGroup=npStd(dictionnaireData[m])
            Data.append(meanGroup)
            Std.append([meanGroup-2*stdGroup,meanGroup+2*stdGroup])
            print("Number of samples: {}  Mean {:.3}  std {:.3}" .format(m,float(meanGroup),float(stdGroup)))
            
            UserResu.write("{} ; {:.3} ; {:.3} \n".format(m,float(meanGroup),float(stdGroup)))
            
            
        plt.plot(Echelle,Std,"_k")
        plt.title("Repetitions of the ApaOxIS Monte Carlo : \nRisk assessment and samples number")
        plt.ylabel("Mean and 95% Confidence interval")
        plt.xlabel("Number of samples")
        plt.plot(Echelle,Data,".k")
        plt.show()
        plt.plot(Echelle,Std,"_r")
        
        plt.plot(Echelle,Data,".k")
        plt.savefig(fileout+".Graphe.pdf")
        tjob2 = time.time()
        UserResu.write("Time to final results : {:.2f} seconds".format(tjob2-tjob0))
        UserResu.write("Individual results are in the directory : ApaOxIS_workplace")
        UserResu.write("\n#################################################################################################\n")
        print ("##### Time to final results : {:.2f} #####".format(float(tjob2-tjob0)))
        print("\n\n############################# THE RESULTS ARE READY  #############################################\n")
        
        print("You will find it HERE ==> "+fileout+"  and the graph is HERE: " +fileout+".Graphe.pdf")
        print("                                                                         ")


#### les tests ####
def skew_test(cell,alpha):#The sample size must be at least 8.
    """
    [from the Scipy documentation] Compute the sample skewness of a data set. For normally distributed data, the skewness should be about zero. For unimodal continuous distributions, a skewness value greater than zero means that there is more weight in the right tail of the distribution. The function skewtest can be used to determine if the skewness value is close enough to zero, statistically speaking.
    This function tests the null hypothesis that the skewness of the population that the sample was drawn from is the same as that of a corresponding normal distribution.
    if skewness is between -0.5 and 0.5, the distribution is approximately symmetric
    """
    stask, pask = skewtest(cell)
    
    
    if len(cell) > 8 : interpretation= " The nul hypothesis HO is that the skew is different from the normal distribution (for normally distributed data, the skewness should be about zero)\n H0 (skew from normal distrib.) rejected" if pask < alpha else " H0 (skew from a normal distrib.) cannot be rejected" #
    else: interpretation= "Sample not large enough for the skew test (less than 8)"
    
    if stask >3 or stask <-3 : notesSK=" The data are excessively skewed !!!"#The kurtosis is normalized so that it is zero for the normal distribution
    elif stask < -1.0 or stask > 1.0  : notesSK=" The data are highly skewed !!"
    elif stask < -0.5 or stask > 0.5 : notesSK=" The data are moderately skewed !"
    else : notesSK=" The data are fairly symmetrical as expected for a Normal law **"
    
    
    
    return (stask, pask ,interpretation,notesSK)
    
    """
    RULES FOR SKEWNESS/KURTOSIS
    -If the skewness is between -0.5 and 0.5, the data are fairly symmetrical
    -If the skewness is between -1 and – 0.5 or between 0.5 and 1, the data are moderately skewed
    
    -If the kurtosis is close to 0, then a normal
    distribution is often assumed (mesokurtic).
    -If the kurtosis is less than zero, then the
    distribution has light tails (platykurtic).
    -If the skewness is less than -1 or greater than 1, the data are highly skewed
    -If the kurtosis is greater than zero, then the
    distribution has heavier tails (leptokurtic).
    """
def kurtosis_test(cell,alpha):#Valid only for n>20.
    """
    [from the Scipy documentation] Compute the kurtosis (Fisher or Pearson) of a dataset.
    Kurtosis is the fourth central moment divided by the square of the variance. If Fisher’s definition is used, then 3.0 is subtracted from the result to give 0.0 for a normal distribution.
    In Fisher’s definiton, the kurtosis of the normal distribution is zero.
    Test whether a dataset has normal kurtosis.

    This function tests the null hypothesis that the kurtosis of the population from which the sample was drawn is that of the normal distribution. Valid only for n>20.
    """
    staku, pku = kurtosistest(cell)
    
    if len(cell) >20  : interpretation= " H0 (kurtosis from a normal distrib.) cannot be rejected" if pku < alpha else "H0 (normal distrib.) cannot be rejected"
    else : interpretation= "Sample not large enough for the skew test (less than 20)"
    if staku > 1 : notesKU=" The distribution is too peaked !!!"
    elif staku <-1 : notesKU=" The distribution is too flat !!!"
    elif staku >0.5 : notesKU=" The distribution is moderately peeked !!"
    elif staku <0.5 : notesKU=" The distribution is moderately flat !!"
    else : notesKU=" Shape close to that expected for a Normal law **"
    return (staku, pku,interpretation,notesKU)
    

def ShapiroTest(cell,alpha):
    """The Shapiro-Wilk test is based on the test statistic called W that may be interpreted as the coefficient of determination between quantiles from the Normal law and empirical quantiles obtained from the data. The Shapiro-Wilk test tests the null hypothesis that the data was drawn from a normal distribution.
    Note that the computing time is increased for the low numbers and high sparcity of the samples.
    """
    stat, psh = shapiro(cell)
    inter="H0 (normal distrib.) rejected" if psh < alpha else "H0 (normal distrib.) cannot be rejected"
    flag= False if psh < alpha else True
    return  (stat,psh,inter,flag)
    
def AndersonTest(cell,alpha):
    """
    The Anderson-Darling test is based on the distance between the empirical distribution function and the hypothetic Normal distribution.
    AD tests is one in the family of the quadratic class of the EDF statistic as it used squared differences.
    Note that the computing time is increased for the low numbers and high sparcity of the samples.
    As the p value level comparison is not relevant, only the putative normality is returned as a boolean
    if the returned statistic (A) is larger than these critical values then for the corresponding significance level (B[2] if α=0.05), the null hypothesis that the data come from the chosen distribution can be rejected. The returned statistic is referred to as ‘A2’ in the references.
    """
    A,B,C = anderson(cell)
    
    inter="H0 (normal distrib.) rejected" if A > B[2] else "H0 (normal distrib.) cannot be rejected"
    flag= False if A > B[2] else True
    return  (A,B[2],C[2],inter,flag)


### la fonction USER (user sample)
def UserStatistics(dataUser,fileIs):
    """
    This funtion is used for the analysis of an user set of samples results.
    1) Basically it performs computation of mean, variance and skewness/kurtosis
    And the skewness/kurtosis test with H0 == Normal law
    2) Then it applied the Shapiro-Wilk and Anderson-Darling tests for normality
    It analyse the Shapiro-Wilk and Anderson-Darling results by using the same rule as in MC
    So that the submitted sample may be classify accordingly
    3) Then it compures the estimation of the risk of mis-interpretation by using our MC approach
    
    The outputs are on the terminal (ugly as usual)
    and in files.
    The dataUser=user_sample_name file(a text file with a list of δ180 values, one per line
    file_is= the name of the output files, say xxxx
    There are the following outputs :
    xxxx.txt = textual output with a lot of informations concerning statistics, methods, explainations and interpretation but with rather primitive format
    xxxx.GraphSummary.pdf = a graphic summary and some limited outputs (statistics and interpretation)
    xxxx.histogram.pdf = the histogram of the δ18O file
    xxxx.fittingNormal.pdf = the fitting of the user data to the Normal model superposed to the histogram
    xxxx.qqplot.pdf the corresponding quantiles-quantiles graph
    
    Note that there are 3 levels for outputs : screen (extensive and ugly), text file (extensive with explainations and some formating), graphic summary (a summary of course)
    """
    #global UserResu #to be used elsewer to write the result, easier to declare a global variable
    print("\n\n#####################################  ApaOxIS ###################################################\n")
    print("\n\n###############################  NORMALITY STATISTICS ############################################\n")
    with open(fileIs, "w") as UserResu: pass #emptying if redoing the analysis do not add the results to the file
    with open(fileIs, "a") as UserResu:
        UserResu.write("#####################################  ApaOxIS ###################################################\n\n")
        
        UserResu.write("                    Normality tests according of bio-apatite contents                                     \n")
        UserResu.write("                       C. Lécuyer & JP. Flandrois 2021\n\n")
        UserResu.write("Publication: A statistical toolbox designed to help detecting the alteration of the original oxygen \nisotope composition of bioapatites C. Lécuyer & JP. Flandrois 2022\n<reference of the program repository>\n\n")
        UserResu.write("\n###################################  RESULTS  ###################################################\n\n")
        UserResu.write("for the submitted sample: "+dataUser+" \n")
        UserResu.write("\n#################################################################################################\n\n")
        
        
        with open(dataUser,"r") as originalfile:
            OF=originalfile.read().strip().replace(',','.').split('\n')
        OF2=[]
        for u in OF: OF2.append(float(u))
        print (min(OF2),max(OF2))
        #verifications
        if min(OF2) < 10.0 or max(OF2) > 25.0:
            txt = input("\n\n\n\n\nOne or more values of δ18O are outside ]10.0,25.0[ \n\nPlease validate that your file is only containing apatite δ18O values\n\nValidate the file (Y/N)\n\n\n")
            if txt.upper()=="Y" :
                print(f"Thank you for the confirmation")
            else:
                print(f"please submit your corrected sample file")
                stop("Non conformity of the user file, see you soon.")
        alpha=0.05
        OFarray=npArray(OF2)
        nbSubmitted=len(OFarray)
        
        #### STATISTICS ####
        #compute
        n, (smin, smax), sm, sv, ss, sk =describe(OFarray)#basic mean, std, skew etc. all direct in scipy
        sk,skt,skinterpretation, notesSK=skew_test(OFarray,alpha) #tests skew
        ku,kut,kuinterpretation, notesKU=kurtosis_test(OFarray,alpha) # kurtosis
        #render
        RSUinit=(f'## BASIC STATISTICS\n\n')
        print (RSUinit)
        UserResu.write(RSUinit)
        RSUbasic=(f'* Number of samples:  {n} \n * Mean: {sm:.3f} \n * Variance: {sv:.3f} \n \n ---- \n * Skew: {ss:.3f} test statistics z-score: {sk:.3f} p-value: {skt:.3f} {skinterpretation} \n    ** {notesSK}\n ---- \n * Kurtosis {ku:.3f} test: z-score: {kut:.3f} p-value: {skt:.3f} {kuinterpretation}\n    ** {notesKU}\n')
        RSU_0=(f'Number of samples:  {n} Mean: {sm:.3f} Variance: {sv:.3f}')
        RSU_1a=(f'Skew: {ss:.3f} test statistics z-score: {sk:.3f} p-value: {skt:.3f}')
        RSU_1b=(f'Kurtosis {ku:.3f} test: z-score: {kut:.3f} p-value: {skt:.3f}')
        RSU_2a=(f'{skinterpretation}')
        RSU_3a=(f'{kuinterpretation}')
        RSU_2b=(f'{notesSK}')
        RSU_3b=(f'{notesKU}')
        print (RSUbasic)
        UserResu.write(RSUbasic)
        #### NORMALITY TESTS ####
        """
        Departure from normality due to skewness or kurtosis is more taken into account with the Shapiro-Wilk test and Aderson-Darling test is sensitive to tail.
        The Shapiro-Wilk test is based on the test statistic called W that may be interpreted as the coefficient of determination between quantiles from the Normal law and empirical quantiles obtained from the data.
        The Anderson-Darling test is based on the distance between the empirical distribution function and the hypothetic Normal distribution.
        AD tests is one in the family of the quadratic class of the EDF statistic as it used squared differences.
        Both tests have demonstrated a high potential in normality detection, but as they are using different approaches it seems interesting to use both.
        """
        #compute
        A,B,C,interpretationAD,flagAD =AndersonTest(OFarray,alpha)
        statsh, psh ,explainSH,flagSH= ShapiroTest(OFarray,alpha)
        #render (note that this is done step by step and in two phases 1) building the sentence 2 printing and sending to the file)
        RSU=(f'\n## NORMALITY TESTS \n\n[The null hypothesis HO is that the observed distribution follows the Normal (gaussian) model\ni.e the sample was drawn from a normal distribution]\n\n')
        print (RSU)
        UserResu.write(RSU)
        RSU=(f'Departure from normality due to skewness or kurtosis is more taken into account with the Shapiro-Wilk (SW) test,\nthe Aderson-Darling (AD) test is much sensitive to tail.\n')
        print (RSU)
        UserResu.write(RSU)
        RSUSW=(f'                \nShapiro-Wilk\n     statistics: {statsh:.3f} p-value: {psh:.3f} {explainSH} \n ')
        print (RSUSW)
        UserResu.write(RSUSW)
        RSUAD=(f'                \nAnderson-Darling\n     statistics: {A:.3f} critical value α: {B:.3f} {interpretationAD} \n ')
        print (RSUAD)
        UserResu.write(RSUAD)
        print("\n## CONVERGENCE OF THE NORMALITY TESTS\n\n")
        UserResu.write("\n## CONVERGENCE OF THE NORMALITY TESTS\n\n")
        #flagAD flagSH flagAGPEAR flagKS
        valueFlag=0
        RSU=''
        #Succinct rules to use the convergence of the tests as in the MC process
        if flagAD and flagSH :#and flagAGPEAR:
            #excellent !
            RSU= (f'All the normality tests are convergent: the sample seems issued from a Normal distributed population\nThis is certainly clear with the graph of the data fitted to a Normal law:{fileIs}.fittingNormal.pdf and the QQ-plot representation{fileIs}.qqplot.pdf')
            RSUGLOB1= (f'Convergence of the tests: Normal distributed population on this basis')
            RSUGLOB2= (f'This is certainly clear with the two graphs upward')
        else:
            if flagAD and not flagSH:
                # bad
                RSU= (f'AD is influenced by the tails and DP is better for short-tailed distribution, more SW is not convergent: the interpretation is tricky and the sample may _not_ be issued from a Normal distributed population\nTherfore a critical analysis of the graph of the data fitted to a Normal law: {fileIs}.fittingNormal.pdf and the QQ-plot representation {fileIs}.qqplot.pdf is _mandatory_')
                RSUGLOB1= (f'Divergence of the tests: Normal distributed population rejected on this basis')
                RSUGLOB2= (f'A critical analysis of the two graph upward may explain')
            else: #very bad
                RSU= (f'The sample is almost certainly _not_ issued from a Normal distributed population\nTo be convinced, have a look to the graph of the data fitted to a Normal law: {fileIs}.fittingNormal.pdf\n and the QQ-plot representation {fileIs}.qqplot.pdf\n')
                RSUGLOB1= (f'The sample is almost certainly _not_ issued from a Normal distributed population')
                RSUGLOB2= (f'This is should be observed in the two graph upward')
        print (RSU)
        UserResu.write(RSU)
        ## risk estimation for a sample nb == user set  COMPUTING THE LINKED ERROR RISK ##
        UserResu.write("\nCOMPUTING THE CORRESPONDING RISK FOR "+str(n)+" SAMPLES ... \n\n")
        UserResu.flush()
        print("\n\n########################  COMPUTING THE LINKED ERROR RISK ########################################\n")
        ##sending to the risk analysis module (Monte-Carlo)
        MCnb,Nb_poputest,Risk_Level,EstSD=MC_run(nbSubmitted,fileIs) #a simple run
        ##rendering
        riskline1=(f'Sample Nb: {Nb_poputest}')
        riskline2=(f'With such a sample ({Nb_poputest}) of bio-apatite being submitted to standard conditions of diagenesis')
        riskline2=(f'There is a {Risk_Level:.3f}+/-{EstSD:.3f} % risk of normality of δ18O data by chance for a {Nb_poputest} sample count.')
        
        ### GRAPHS ###
        
        ##HISTOGRAM
        ## the colorful elegant solution is taken from Jeeteshgavande30
        ## https://www.geeksforgeeks.org/plotting-histogram-in-python-using-matplotlib/
        ## its interest is to colorize the highest parts so that we immediately see if it is two-peeked
        from matplotlib import colors
        from matplotlib.ticker import PercentFormatter
        # Creating histogram
        fig, axs = plt.subplots(1, 1,
                                figsize =(10,7),
                                tight_layout = True)
        #figure(figsize=(11.69,8.27)) # for landscape
        # Remove axes splines
        for s in ['top', 'bottom', 'left', 'right']:
            axs.spines[s].set_visible(False)
        # Remove x, y ticks
        axs.xaxis.set_ticks_position('none')
        axs.yaxis.set_ticks_position('none')
        # Add padding between axes and labels
        axs.xaxis.set_tick_params(pad = 5)
        axs.yaxis.set_tick_params(pad = 10)
        # Add x, y gridlines
        axs.grid(visible = True, color ='grey',
                linestyle ='-.', linewidth = 0.5,
                alpha = 0.6)
         
        # Add Text watermark
#        fig.text(0.9, 0.15, 'Jeeteshgavande30 https://www.geeksforgeeks.org/plotting-histogram-in-python-using-matplotlib/',
#                 fontsize = 12,
#                 color ='red',
#                 ha ='right',
#                 va ='bottom',
#                 alpha = 0.7)
        # Creating histogram
        
        N, bins, patches = axs.hist(OFarray)#, bins="doane")#, bins = len(OFarray/10))
        # Setting color
        fracs = ((N**(1 / 5)) / N.max())
        normal = colors.Normalize(fracs.min(), fracs.max())
        for thisfrac, thispatch in zip(fracs, patches):
            color = plt.cm.viridis(normal(thisfrac))
            thispatch.set_facecolor(color)
        # Adding extra features
        plt.xlabel("δ18O apatite")
        plt.ylabel("Population")
        #plt.legend(legend)
        plt.title('δO18 apatite from '+dataUser)
        # Show plot #testing
        #plt.show() #testing
        plt.savefig(fileIs+".histogram.pdf")
        plt.close()
        
        ############################### residuals
        from numpy import histogram
        
        y, x = histogram(OFarray, density=True)# bins="doane"
        # Milieu de chaque classe
        from numpy import roll
        x = (x + roll(x, -1))[:-1] / 2.0
        print (x,y)
        from lmfit.models import GaussianModel
        mod = GaussianModel()
        pars = mod.guess(y, x=x)
        print ("A",pars)
        
        out = mod.fit(y, pars, x=x)
        
        X = linspace(min(x),max(x), 100)
        Y = mod.eval(out.params,x=X)
        print (X,Y)
#        plt.plot(X, Y,'-r')
#        plt.show()
#        stop()
        print(out.fit_report(min_correl=0.25))
        
        plt.plot(x, y,'o',label='Observed δ18O')
        plt.plot(X, Y,'-r',label='Gaussian fit')
        dely = out.eval_uncertainty(x=x,sigma=3)#Evaluate the uncertainty of the model function. This can be used to give confidence bands for the model from the uncertainties in the best-fit parameters.
        #https://cars9.uchicago.edu/software/python/lmfit/model.html
        plt.fill_between(x, out.best_fit-dely,out.best_fit+dely, color='#BEC7BD',label='3-$\sigma$ uncertainty band')
        #plt.plot(x, out.init_fit, '--', label='initial fit')
        #plt.plot(x, out.best_fit, 'o', label='best fit')
        
        plt.legend()
        plt.savefig(fileIs+".UncertaintyfittingNormal.pdf")
        plt.close()
        print("plot residuals")
        out.plot_residuals()
        plt.savefig(fileIs+".UncertaintyResiduals.pdf")
        plt.close()
#        stop()
        print ("B",out.params)
        #print (linspace(set(x)))
        
        ###ma solution mais c'est fait au dessus automatiquement :)
#        r=mod.eval(out.params,x=x)
#        print (r)
#        residuals=y-r
#        print ('residuals:',residuals)
#        plt.show()
#        plt.close()
#        pars = mod.guess(residuals, x=x)
#        out = mod.fit(residuals, pars, x=x)
#        print('RESIDUALS\n',out.fit_report(min_correl=0.25))
#        neutralLinex = npArray([min(x),max(x)])
#        neutralLiney = npArray([0.0,0.0])
#        print (neutralLiney,neutralLinex)
#
#        plt.plot(neutralLinex,neutralLiney, linestyle = 'dotted')
#        plt.plot(x, residuals,'o')
#        #plt.plot(x, out.best_fit, '-', label='best fit')
#        #plt.legend()
#        plt.show()
        
        ############################# fin residuals#
        
        
        ### FITTING TO NORMAL DISTRIBUTION + SIMPLE HISTOGRAM
        ### Fit a normal distribution to the user-data ##
        fig, axs = plt.subplots(1, 1,
                                figsize =(10, 7),
                                tight_layout = True)
        mu, std = Norm.fit(OFarray)
        # Histogram easy way
        plt.hist(OFarray, density=True, alpha=0.6, color='g')
        # PDF
        xmin, xmax = plt.xlim()
        x = linspace(xmin, xmax, 100)
        p = Norm.pdf(x, mu, std)
        plt.plot(x, p, 'k', linewidth=2)
        plt.xlabel("δ18O apatite")
        plt.ylabel("Population")
        title = "δO18 apatite from "+dataUser+": Fitted to Normal lax (mu = %.2f,  std = %.2f)" % (mu, std)
        plt.title(title)
        plt.savefig(fileIs+".fittingNormal.pdf")
        plt.close()
        ### QUANTILE-QUANTILE PLOT
        # QQ-plot
        fig, axs = plt.subplots(1, 1,
                                figsize =(10, 7),
                                tight_layout = True)
        probplot(OFarray, dist="norm", plot=plt)
        #fig.suptitle("ApaOXIS")
        axs.set_title("Probability plot for "+dataUser)
        axs.get_lines()[0].set_marker('o')
        axs.get_lines()[0].set_markerfacecolor('g')
        axs.get_lines()[0].set_markerfacecolor('g')
        axs.get_lines()[0].set_markersize(5.0)
        axs.get_lines()[1].set_color("g")
        plt.savefig(fileIs+".qqplot.pdf")
        plt.close()
        
        ###SUMMARY GRAPH + TEXT
        #graph prob plot of global sales
        plt.figure(figsize = (10,14))
        plt.suptitle("Datas from "+dataUser)
        # Histogram easy way
        ax1 = plt.subplot(311)
        mu, std = Norm.fit(OFarray)
        #tests
        #bins =histogram_bin_edges(OFarray, bins="auto", range=None, weights=None)
        #ax1=plt.hist(OFarray, density=True, bins=bins, alpha=0.6, color='g')
        #end tests
        ax1=plt.hist(OFarray, density=True , alpha=0.6, color='g')#, bins="doane"
        # PDF model
        xmin, xmax = plt.xlim()
        x = linspace(xmin, xmax, 100)
        p = Norm.pdf(x, mu, std)
        ax1=plt.plot(x, p, 'k', linewidth=2)
        ax1=plt.xlabel("δ18O apatite")
        ax1=plt.ylabel("Population")
        title = "δO18 apatite Fitted to Normal lax (mu = %.2f, std = %.2f)" % (mu, std)
        ax1=plt.title(title)
        #graph prob plot
        ax2 = plt.subplot(312)
        res = probplot(OFarray, plot=plt)
        ax2.set_title('Probability Plot the samples')
        ###TEXT SUMMARY
        ax3 = plt.subplot(313)
        ax3.set(xlim=(0, 10), ylim=(0, 10))
        ax3.text(0, 1,"SUMMARY",horizontalalignment='left',verticalalignment='center',fontweight='bold', fontvariant='small-caps',transform = ax3.transAxes)
        ax3.text(0.01, 0.93,"BASIC STATISTICS",horizontalalignment='left',verticalalignment='center',transform = ax3.transAxes)
        ax3.text(0, 0.87,RSU_0,horizontalalignment='left',verticalalignment='center',transform = ax3.transAxes)
        ax3.text(0, 0.82,RSU_1a,horizontalalignment='left',verticalalignment='center',transform = ax3.transAxes)
        ax3.text(0.05, 0.78,RSU_2a,horizontalalignment='left',verticalalignment='center',transform = ax3.transAxes)
        ax3.text(0.05, 0.73,RSU_2b,horizontalalignment='left',verticalalignment='center',transform = ax3.transAxes)
        ax3.text(0, 0.68,RSU_1b,horizontalalignment='left',verticalalignment='center',transform = ax3.transAxes)
        ax3.text(0.05, 0.62,RSU_3a,horizontalalignment='left',verticalalignment='center',transform = ax3.transAxes)
        ax3.text(0.05, 0.58,RSU_3b,horizontalalignment='left',verticalalignment='center',transform = ax3.transAxes)
        ax3.text(0.01, 0.5,"NORMALITY TESTS",horizontalalignment='left',verticalalignment='center',transform = ax3.transAxes)
        ax3.text(0.02, 0.45,RSUSW.replace('                \n','').replace('\n',''),horizontalalignment='left',verticalalignment='center',transform = ax3.transAxes)
        ax3.text(0.02, 0.40,RSUAD.replace('                \n','').replace('\n',''),horizontalalignment='left',verticalalignment='center',transform = ax3.transAxes)
        ax3.text(0.04, 0.35,RSUGLOB1,horizontalalignment='left',verticalalignment='center',fontweight='semibold',transform = ax3.transAxes)
        ax3.text(0.01, 0.25,"RISK ESTIMATION (Monte-Carlo diagenesis in-silico "+str(MCnb)+" iterations), sample count: "+str(Nb_poputest),horizontalalignment='left',verticalalignment='center',transform = ax3.transAxes)
        ax3.text(0, 0.15,"Mean and Std.Dev. estimated by the bootstrap method",horizontalalignment='left',verticalalignment='center',transform = ax3.transAxes)
        ax3.text(0, 0.10,riskline2,horizontalalignment='left',verticalalignment='center',fontweight='semibold',transform = ax3.transAxes)
        ax3.text(0, 0.00,"Computed with ApaOsIS C. Lécuyer & JP. Flandrois 2022",horizontalalignment='left',verticalalignment='center',transform = ax3.transAxes)
        #ax3.text(0, 0.01,"A statistical toolbox designed to help detecting the diagenesis...",horizontalalignment='left',verticalalignment='center',transform = ax3.transAxes)
        ax3.axis('off')
        plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)
        plt.savefig(fileIs+".GraphSummary.pdf")
        print("\n\n########################  THE RESULTS AND GRAPHS ARE READY #######################################\n")
        
        print("You will find it HERE ==> the results as a text file:"+fileIs+"\nthe histogram:"+fileIs+".histogram.pdf"+"\nthe data fitted to a Normal law:"+fileIs+".fittingNormal.pdf"+"\nthe QQ-plot representation"+fileIs+".qqplot.pdf")
        print("                                                                         ")
        print("\n\n##################################  ApaOxIS END ##################################################\n")
        
def stop(valeur='\n-------------------\nPROGRAM ABORTED'):
    #a very primitive error-catching during the debugging
    print('\n-------------------\nPROGRAM ABORTED: \n'+valeur)
    exit(1)
    
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(prog="ApaOxIS",description='''
    APAtite OXygen Isotope Simulation (a toolbox)
    ApaOxIS has been designed to
    1) help to improve the sampling strategy when working with δ18O in biological apatites, it uses Monte-Carlo simulation of the alteration of the original oxygen isotope composition of biogenic apatites (-m) to compute an estimation of the risk to observe a normally distributed sample by chance
    2) run the normality test of δO18 apatite data (-s) by using the same method as in (1) and the risk to observe a normally distributed sample by chance for the number of samples.
    ''',
    epilog="""
    place our references here (paper, site...)
    """)

    #print("\n\n#####################################  ApaOxIS ###################################################\n")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-r", "--Risk", action='store_true', help="Monte-Carlo Simulation of sampling scheme (sample size) -output as error rate for a given sampling scheme-")
    group.add_argument("-x", "--xRepeatedMC", action='store_true', help="30x Repeated Monte-Carlo Simulation Simulation of sampling scheme (sample size) final graph on the collected results")
    group.add_argument("-s", "--Normality", action='store_true', help="Statistics for Normality testing on users data and computing the error rate using MC with the sample size corresponding to the data")
#    group.add_argument("-m", "--Model", action='store_true', help="Monte-Carlo Simulation -output as histograms of frequencies of δ18Oa- ")
    parser.add_argument("i", type=str, help="the parameter file (yaml file, WARNING: structure is different between -m / others); or the data file (one column plain text)")
    parser.add_argument("o", type=str, help="output results file")
    args = parser.parse_args()
    #print (args)
    if args.Risk:
        print("\n\n#####################################  ApaOxIS ###################################################\n")
        MC_run(args.i,args.o)
    if args.xRepeatedMC:
        print("\n\n#####################################  ApaOxIS ###################################################\n")
        tjob0 = time.time()
        #  for the publication statistics only. Usualy not used
        #creation of a result directory
        createFolder("ApaOxIS_workplace")
        WAPfile=os.path.join("ApaOxIS_workplace","tempo")
        for i in range(0,30): #for the publication statistics
            MC_run(args.i,WAPfile+str(i))
        multiMCstatistics("ApaOxIS_workplace", args.o,tjob0)
    if  args.Normality:
        UserStatistics(args.i,args.o)
