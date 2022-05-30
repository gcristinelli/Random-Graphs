import numpy as np
import math
import random 
import matplotlib.pyplot as plt


n=500
r=5
p=1
nruns=50

def initial_adjacency_matrix(n,r):     
    adj_mat = np.zeros([n,n])
    for i in range(0,n):
        #left
        if(i-r>=0):
            adj_mat[i][i-r:i] = 1
        else:
            diff = r-i
            adj_mat[i][0:i] = 1
            adj_mat[i][n-diff:n+1] = 1
        #right
        if(i+r<n):
            adj_mat[i][i+1:i+r+1] = 1 #+1 to avoid self loop and up to sym value
        else:
            diff = i+r-n
            adj_mat[i][i+1:n+1] = 1
            adj_mat[i][0:diff+1] = 1

    return adj_mat
    
def rewired_adjacency_matrix(n,r,p):
    adj_mat=initial_adjacency_matrix(n,r)
    for i in range(0,n):
        #print("computing row # %d/%d"%(i+1,n))
        for j in range(0,r):  #for each link to vertex i
            if (random.random()<p): #attempt a rewire
                #performing the rewire
                #    - Choose which of the connected edge to rewire->deleted_edge
                #    - Choose were to rewire it among the available positions->candidates
                #    - Perform the connection/delete old connection/update mirror adjmat

                #choose which edge to remove: [+periodic boundary conditions]
                deleted_edge = i+1+j
                if deleted_edge>n-1:
                    deleted_edge = deleted_edge-n
                
                #chose available position:
                candidates = np.where(adj_mat[i]==0)
                candidates=np.delete(candidates[0], np.where(candidates[0]==i)) #no self loop
                new_edge = random.choice(candidates)

                #print("candidates list = ",candidates)
                #print("new edge chosen = ",new_edge)
                
                #create new wire
                adj_mat[i][new_edge]=1
                adj_mat[new_edge][i]=1

                #deleate old wir	e
                adj_mat[i][deleted_edge]=0
                adj_mat[deleted_edge][i]=0
                
    return adj_mat
    
def degree(n,r,p,nruns):
    dens=np.zeros(n)
    deg=np.zeros([nruns,4*r])
    for k in range(0,nruns):
        print("computing degree distribution of run=%d/%d"%(k+1,nruns))
        adj_mat=rewired_adjacency_matrix(n,r,p)
        for i in range(0,n):
            dens[i]=np.sum(adj_mat[i]) 
        for i in range(0,4*r):
            deg[k][i]=np.count_nonzero(dens == i)/n
    return deg

def plot_degree(avg,std, compare=True):
    plt.figure(1,figsize=[10,10])
    ax = plt.gca()
    ax.set_xticks(np.arange(0, 4*r, 1))
    plt.grid(b=True, color='gray', ls=':')
    ax.set_ylim([-.01,0.35])
    #ax.xaxis.tick_top()
    T_plot = np.linspace(0,4*r-1,4*r)
    plt.errorbar(T_plot,avg,yerr=std,fmt='.',ls=':', color='blue',ecolor='green', elinewidth=1, capsize=2)
    if compare==True:
    	P=np.zeros(4*r)
    	for m in range(r,2*r):
    	    for n in range (0,m-r+1):
    	        ex=1/math.exp(p*r)
    	        fact=1/(math.factorial(n)*math.factorial(r-n)*math.factorial(m-r-n))
    	        P[m]=P[m]+math.factorial(r)*(1-p)**n*p**(r-n)*(r*p)**(m-r-n)*ex*fact
    	for m in range(2*r,4*r):
    	    for n in range (0,r+1):
    	        ex=1/math.exp(p*r)
    	        fact=1/(math.factorial(n)*math.factorial(r-n)*math.factorial(m-r-n))
    	        P[m]=P[m]+math.factorial(r)*(1-p)**n*p**(r-n)*(r*p)**(m-r-n)*ex*fact
    	plt.plot(T_plot,P,marker='.',ls='-.', color='red')
    plt.savefig("ComparedDegreeDistribution_n=%d_r=%d_p=%d_nruns=%d.png" % (n,r,p*100,nruns), delimiter=",",dpi=250)
    plt.show()

def plot_adjacency_matrix(n,r,p,adj_mat):
    plt.figure(3,figsize=[10,10])
    ax = plt.gca()
    ax.xaxis.tick_top()
    plt.imshow(adj_mat, cmap='binary', interpolation='nearest')
    plt.savefig("adjacency_n=%d_r=%d_p=%d.png" % (n,r,p*100), delimiter=",",dpi=250)

    	
   
DD=degree(n,r,p,nruns)
DD_avg=np.average(DD,axis=0)
DD_std=np.std(DD,axis=0)
plot_degree(DD_avg,DD_std)
#np.savetxt("n=%d_r=%d_p=%d_adjacency_matrix.csv" % (n,r,p), adj_mat, delimiter=",")
#plot_adjacency_matrix(n,r,p,adj_mat)



def circle(n):
    #Returns x,y coordinates of points on an unit circle with spacing 2π/n
    x,y = [],[]
    step_angle = -2*math.pi/n
    for i in range(0,n):
        x.insert(len(x),math.cos(i*step_angle))
        y.insert(len(y),math.sin(i*step_angle))
    return x,y
    
def plot_graph(n,p,r,adj_mat):
    plt.figure(2,figsize=[10,10])
    plt.axis('off')
    ax = plt.gca()
    coords = circle(n) 
    labels_nodes = np.zeros(n)
    for i in range(0,n):
        connections = np.where(adj_mat[i]==1)
        print("plotting node # %d/%d"%(i+1,n))
        for k in range(0,len(connections[0])):
            #print("    wiring edge # %d/%d"%(k+1,len(connections[0])))
            ax.plot([coords[0][i],coords[0][connections[0][k]]],[coords[1][i],coords[1][connections[0][k]]],linewidth=.5,color='green')
    ax.plot(coords[0],coords[1],color='black',ls='none',marker='.',markersize=5)
    plt.savefig("graph_n=%d_r=%d_p=%d.png" % (n,r,p*100), delimiter=",",dpi=250)
    #plt.show()
    
#plot_graph(n,p,r,adj_mat)
    
    
    
    
    
