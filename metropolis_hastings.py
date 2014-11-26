"""

Metropolis-Hastings algorithm implementation

"""
from numpy import *


def generate_neighbors(macierz_przek,a):
    """
    for the rotate matrix generating all possible next microstates
    """
    s=[]
    ma = a * macierz_przek
    for i in xrange(a.shape[0]):
        s.append(vstack((a[:i,:],
                        ma[i:,:])))
    return s

def neighbors(a,rotate_matrices):
    """
    generating all possible next microstates for microstate a
    """
    l=[]
    for macierz_przek in rotate_matrices:
        s = generate_neighbors(macierz_przek,a)
        l +=s       
    return l


def absolute_coordinates(a): 
    """
    calculation of absolute coordinates in 2D
    i.e. adding two first rows to current matrix and cumsum on it 
    """
    a = vstack((array([[0,0],[1,0]]),a)) #back to the original matrix
    return a.cumsum(axis=0) #add coordinates "cumulative" around the x-axis


def unique_rows(c):
    return array(list(set(map(tuple, c.tolist()))))

def is_shape_allowed(a):
    c = absolute_coordinates(a)
    return unique_rows(c).shape == c.shape #True if shape is allowed, i.e. the polymer does not cross
    
    
def check_allowed_neighbors(l):
    return [neighbor for neighbor in l if is_shape_allowed(neighbor)]
        
    
def number_of_contacts(a,matrix_l):
    d = 0
    #left = array([-1,0]) #there is no need
    right = array([1,0])
    #down = array([0,-1]) #there is no need
    up = array([0,1])

    c = absolute_coordinates(a)
    
    m_logic = matrix_l == "H" #True when hydrophobic, false when polar 'row' in matrix
    m_hydro = c[m_logic,:]
    for w in [right,up]:
        m_moved = m_hydro+w
        m_wspolna = unique_rows(vstack((m_hydro,m_moved))) 
        d += (m_hydro.shape[0] + m_moved.shape[0]) - m_wspolna.shape[0]
    return d


def h_next2eachother(matrix_l):
    n=0
    for i in xrange(len(matrix_l)-1):
        if matrix_l[i]=="H" and matrix_l[i+1]=="H":
            n+=1
    return n
    

def how_many_contacts_opposite(matrix_l,a):
    return number_of_contacts(a,matrix_l) - h_next2eachother(matrix_l)


  
def count_energy(microstate,delta,matrix_l):
    return -1 * delta * how_many_contacts_opposite(matrix_l,microstate)
    
       
def create_new_microstate(microstate,Lmacierzy):
    neighbors1 = neighbors(microstate,Lmacierzy)
    neighbors_allowed = check_allowed_neighbors(neighbors1)
    m = random.random_integers(0,len(neighbors_allowed)-1)
    return len(neighbors_allowed),neighbors_allowed[m]
    
    
    
    
    
def Metropolis_Hastings(K,T,microstateX,rotate_matrices,delta,matrix_l):
    
    """ 
        X-old microstate, Y-new microstate
    """
    
    A = [] #all microstates
    n_of_contacts = []
    energy = []
    
    
    while K:
        
        energyX = count_energy(microstateX,delta,matrix_l)   
        n_of_contacts_X = how_many_contacts_opposite(matrix_l,microstateX)
        
        n_neighborsX,microstateY = create_new_microstate(microstateX,rotate_matrices)

        energyY = count_energy(microstateY,delta,matrix_l)
        n_of_contacts_Y = how_many_contacts_opposite(matrix_l,microstateY)
        
        n_neighborsY,microstateZ = create_new_microstate(microstateY,rotate_matrices) #microstate Z is not needed to further calculations
        
        prop_accept_microstateY = (-energyY/T) - math.log(n_neighborsY) - (-energyX/T) + math.log(n_neighborsX)
        random_number = math.log(random.random(1)[0]) #e.g.array([ 0.25290701])

        if random_number < prop_accept_microstateY: # both numbers are logarithmized
            #creating new state
            A.append(microstateY)
            n_of_contacts.append(n_of_contacts_Y)
            energy.append(energyY)
            microstateX = microstateY
        else:
            #old state stays
            A.append(microstateX)
            n_of_contacts.append(n_of_contacts_X)
            energy.append(energyX)
        K-=1
    return A,n_of_contacts,energy
              