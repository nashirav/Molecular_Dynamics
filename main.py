from metropolis_hastings_algorithm import *  
from replica_exchange import *


def rotation_matrix(alfa):
    return matrix([[round(math.cos(alfa)),     round(-math.sin(alfa))],
		  [round(math.sin(alfa)),     round(math.cos(alfa))]])


def initialization(l):

  matrix_polymer = array([i for i in l.upper()])

  #initialization
  #a = array((1,0)*len(l)).reshape((len(l),2))
  #but there is no need to take two first coordinates of amino acids:
  start_microstate = array((1,0)*(len(l)-2)).reshape((len(l)-2,2))

  m90 = rotation_matrix(math.pi / 2.0)
  m180 = rotation_matrix(math.pi)
  m270 = rotation_matrix( (3.0/2.0) * math.pi )

  return start_microstate, [m180,m90,m270], matrix_polymer

    
def matrix2coord(m):
    return m.cumsum(axis=0)


def simulation2pdb(simulation, hp_polymer, output_path):
    
    f = open(output_path, "w")
    aa_type = ["lysine" if x=="P" else "alanine" for x in hp_polymer.tostring() ]
    
    for m in xrange(len(simulation)):
      temp, matrix, number_of_contacts, energy = simulation[m][0], simulation[m][1], simulation[m][2], simulation[m][3]
      f.write("MODEL " + str(m + 1) + "\n")
      f.write("COMMENT\t" + str(number_of_contacts))
      
      coords = matrix2coord(matrix)
      for i in xrange(coords.shape[0]):
	c = coords[i].tolist()[0]
	f.write("ATOM      0  CA  " + aa_type[i] + "   0  A     " + str(float(c[0])) + "   " + str(float(c[1])) + "   0.000")
	
      f.write("ENDMDL")
    f.close()





if __name__ == '__main__':
  
    
    print "start replica exchange"
    
    #polymer based on HP model (hydrophobic-polar protein folding model)
    #l="PHPPHPPHHPPHHPPHPPHP"
    #l="HPPPHHPPHPHHPHHH"
    l = "HPPPHHPPHPHHHHHH"

    kb = 1          # Boltzman constant 
    Tmax = 1        # max temperature
    Tmin = 0.15     # min temperature
    delta = 1       # constant used during calculating energy for each microstate

    K = 10          # number of steps of the Metropolis-Hastings algorithm
    symulation_length = 100
    number_of_replicas = 5
    every_x_step_exchange = 10

    T = replicas_temp(Tmax, Tmin, number_of_replicas)

    start_microstate, rotation_matrices, matrix_polymer = initialization(l)
    simulation = replica_exchange(symulation_length,number_of_replicas,every_x_step_exchange,T,delta,K,start_microstate,rotation_matrices,matrix_polymer,kb)    
        
    replica_Tmin = [i[0] for i in simulation]
    simulation2pdb(replica_Tmin, matrix_polymer, "output/trajectory.pdb")

  