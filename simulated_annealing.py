from metropolis_hastings import *
from init_and_save import *


def simulated_annealing(Tmax,Tmin,delta,K,T,microstateX,rotate_matrices,matrix_l,kb,alfa):
    
    T = Tmax
    simulation = []
    
    while T>=Tmin:
      
        print "Temp:",T
        
        aM, number_of_contacts, energy = Metropolis_Hastings(K,T,microstateX,rotate_matrices,delta,matrix_l)
        microstateX = aM[-1] #last microstate from previous step will be the first of current step

        for j in range(len(aM)):
            simulation.append( (T, aM[j], number_of_contacts[j], energy[j]) ) #temp, matrix, number_of_contacts, energy

        T-=alfa

    return simulation


  
if __name__ == '__main__':

    ## polymer based on HP model (hydrophobic-polar protein folding model)
    #l="PHPPHPPHHPPHHPPHPPHP"
    l = "HPPPHHPPHPHHHHHH"
    #l="HPPPHHPPHPHHPHHH"


    K = 200         #number of steps of the Metropolis-Hastings algorithm
    kb = 1          #Boltzman constant 
    Tmax = 1        #initial temperature
    Tmin = 0.15     #final temperature
    delta = 1       #constant used during calculating energy for each microstate
    alfa = 0.05     #temperature will decrease by alfa

    start_microstate, rotation_matrices, matrix_polymer = initialization(l)

    simulation = simulated_annealing(Tmax,Tmin,delta,K,Tmax,start_microstate,rotation_matrices,matrix_polymer,kb,alfa)
    simulation2pdb(simulation, matrix_polymer, "output/trajectory_sa.pdb")


    kb = 1          # Boltzman constant 
    Tmax = 1        # max temperature
    Tmin = 0.15     # min temperature
    delta = 1       # constant used during calculating energy for each microstate