from metropolis_hastings_algorithm import *  
  
  
def replicas_temp(Tmax, Tmin, how_many_replicas):
    deltaT = (Tmax - Tmin) /(how_many_replicas-1)
    return [Tmin + (deltaT * i) for i in xrange(how_many_replicas)] # temperature for each replica
  
     
    
def replica_numbers2exchange(total_number_of_replicas):
    # randomly chosen numbers of replicas to exchange with its neighbour
    a1 = random.randint(0,total_number_of_replicas-2) #minus 2, because I index replicas from 0 and without last one
    return a1, a1+1


def prob_replica_exchange(Ti,Tj,Ei,Ej,kb):
    return (1.0/(kb*Tj) - 1.0/(kb*Ti))*(Ei - Ej)



def replica_exchange(length_of_simulation,number_of_replicas,every_x_step_exchange,T,delta,K,microstateX,rotation_matrices,matrix_polymer,kb):
    """
    The replica exchange algorithm implementation
    """
    
    g = range(0, length_of_simulation, every_x_step_exchange) # each x I allow for the occurrence of exchange
    
    
    simulation = []
    while length_of_simulation: #Monte Carlo step
        print "step of simulation", length_of_simulation
        replicas = []
        
        for t in xrange(len(T)):
            if simulation==[]:
                aM, number_of_contacts, energia = Metropolis_Hastings_algorithm(K,T[t],microstateX,rotation_matrices,delta,matrix_polymer)
            else:
                aM, number_of_contacts, energia = Metropolis_Hastings_algorithm(K,T[t],simulation[-1][t][1],rotation_matrices,delta,matrix_polymer)
           

            microstate,number_of_contacts,energia = aM[-1],number_of_contacts[-1],energia[-1]
            replicas.append((T[t], microstate, number_of_contacts, energia))
              
        
        if length_of_simulation in g:
	    #replica exchange
            
            nr_a1,nr_a2 = replica_numbers2exchange(number_of_replicas) #indexes of replicas that I want to exchange
            a1,a2 = replicas[nr_a1],replicas[nr_a2]
            
            random_number = math.log(random.random(1)[0]) #array([ 0.25290701]); from uniform distribution

            Ti = a1[0]
            Tj = a2[0]
            Ei = a1[2]
            Ej = a2[2]

            p = prob_replica_exchange(Ti,Tj,Ei,Ej,kb) #logarithmized
            if p < random_number: #both logarithmized
                print "prob. of replica exchange=",p,"random number=",random_number
                print "2 replicas has exchanged ( %f , %f )" % (float(replicas[nr_a1][0]),float(replicas[nr_a2][0]))
                replicas[nr_a1] = a2
                replicas[nr_a2] = a1
            else:
                print "prob. of replica exchange=",p, "random number=",random_number
                print "we don't exchange replicas!"

        simulation.append(replicas)    
        length_of_simulation-=1
        
    return simulation
 