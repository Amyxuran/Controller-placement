import numpy as np
import matplotlib.pyplot as plt
import random
import networkx as nx
from decimal import *
print(getcontext())
Context(prec=28, rounding=ROUND_HALF_EVEN, Emin=-999999, Emax=999999, capitals=1, clamp=0, flags=[], traps=[InvalidOperation, DivisionByZero, Overflow])
getcontext().prec = 50

def FindPath(network,n):
    nodes=np.array(range(n))
    minPathwidth=n
    minPathDecomposition=0
    for node in nodes:
        maxVertexSep=0
        pathDecomposition=[]
        pathDecomposition.append(node)
        remainNodes=nodes[nodes!=node]
        
        while len(remainNodes)!=0:
            a=np.nonzero(network[:,pathDecomposition])
            Neighborhood=list(set(a[0]).intersection(set(remainNodes)))
            random.shuffle(Neighborhood)
            tempVertexSep=n
            
            for neighbor in Neighborhood:
                tempPD=pathDecomposition.copy()
                tempPD.append(neighbor)
                tempRN=remainNodes[remainNodes!=neighbor]
                VertexSep=np.count_nonzero(sum(network[tempPD][:,tempRN]))
                
                if VertexSep < tempVertexSep:
                    tempVertexSep=VertexSep
                    chosen=neighbor
            
            pathDecomposition.append(chosen)
            remainNodes=remainNodes[remainNodes!=chosen]
            
            if tempVertexSep > maxVertexSep:
                maxVertexSep=tempVertexSep

        if maxVertexSep < minPathwidth:
            minPathwidth=maxVertexSep
            minPathDecomposition=pathDecomposition

    # print('The minimum pathwidth is ',minPathwidth)
    # print('The decomposition path is ',minPathDecomposition) # with node numbering start from 0
    
    return minPathwidth,minPathDecomposition
# find the decomposition series (activate node/deactivate node/activate link)
def FindSeries(network,Path):
    path=np.add(Path,1) # with node numbering start from 1
    for i in path:
        int(i) 
    series=[]
    DelNodes=[]
    i=0
    for node in path:
        i=i+1
        # activate node
        series.append(node)
        # activate egde
        a=np.nonzero(network[:,node-1])
        Neighbors=np.add(a[0],1) # all neighbors of node
        for Neighbor in Neighbors:
            if Neighbor in path[0:i]:
                series.append([Neighbor,node])

                # deactivate node
                NeighborColumn=network[:,Neighbor-1]
                NeighborColumn=NeighborColumn[path[i:]-1]  # if neighbor is connected to the remain part
                neighborSepSet=sum(NeighborColumn)
                if neighborSepSet==0 and (Neighbor not in DelNodes):
                    series.append(-Neighbor)
                    DelNodes.append(Neighbor)

        # deactivate node
        nodeColumn=network[:,node-1]
        nodeColumn=nodeColumn[path[i:]-1]
        nodeSepSet=sum(nodeColumn)
        if nodeSepSet==0 and (node not in DelNodes):
            series.append(-node)
            DelNodes.append(node)
            
    # print('The decomposition series is',series)  # with node numbering start from 1
    return series
def add_poly(L1,L2):
    R=[]
    if len(L1)>len(L2):
        L1,L2=L2,L1
    i=0
    while i<len(L1):
        R.append(L1[i]+L2[i])
        i+=1
    R=R+L2[len(L1):len(L2)]
    return R

def multiply_poly(L1,L2):
    if len(L1)>len(L2):
        L1,L2=L2,L1
    zero=[];R=[]
    for i in L1:
        T=zero[:]
        for j in L2:
            T.append(i*j)
        R=add_poly(R,T)
        zero=zero+[0]
    return R
def Decomposition(network,n,series,Link_to_add):
    # Reliability from decomposition method
    # Polynomial is stored as [a_0,a_1,...a_n]
    # classes are stored as following example:
    # node class1 class2
    #  1     0      0
    #  2     0      0
    #  3     1      1
    #  4     2      1
    #  5     0      0
    #  6     0      0
    # Meaning: 3 and 4 are activated node, two ways to group them  (3/4) or (34)

    R=[]
    R.append([1,0]) # initial reliability (cannot use [1] at here)
    X=np.zeros((n,2))
    X[:,0]=range(1,n+1) # Create classes storage array
    X[series[0]-1,1]=1 # First step and activate first node

    for s in range(1,len(series)):  # Start from second step in series
        step=series[s]
        if isinstance(step, (int, np.integer)):    # If it is int number, thus is node process
            if s==len(series)-1:
                # print(R)
                z=1
            else:
                if step > 0:
                    # activation of node
                    # print('This step is activation of node',step)
                    for Class in range(1,X.shape[1]): 
                        New_class=len(set(X[:,Class]))
                        X[step-1,Class]=New_class
                    # print(X,R)
                else:
                    # deactivation of node
                    # print('This step is deactivation of node',step)
                    ColumnDel=[]  # record which class make node disconnected
                    for Class in range(1,X.shape[1]):
                        ClassNum_del=X[abs(step)-1,Class]
                        X[abs(step)-1,Class]=0 # deactivate node
                        if ClassNum_del not in X[:,Class]:
                            # this deactivated node is disconnected
                            ColumnDel.append(Class)
                        else:
                            # renumbering
                            ClassNum=len(set(X[:,Class]))
                            b=range(1,ClassNum)
                            X[:,Class]
                            j=0
                            for i in range(len(X[:,Class])):
                                if X[i,Class]>0:
                                    X[:,Class][X[:,Class]==X[i,Class]]=-b[j]
                                    j=j+1
                            X[:,Class]=abs(X[:,Class])
                            
                    X = np.delete(X, ColumnDel, axis=1)
                    RowDel=[i-1 for i in ColumnDel]
                    RowDel.sort(reverse=True)
                    for i in RowDel:
                        del R[i]
                    
                    
                    # Merge the class with same partition
                    for Class in range(1,X.shape[1]):
                        ColumnMerge=[]
                        ColumnMerge.append(Class)
                        for remainClass in range(Class+1,X.shape[1]):
                            if (X[:,Class]==X[:,remainClass]).all():
                                ColumnMerge.append(remainClass)
                        if len(ColumnMerge)>=2:
                            X = np.delete(X, ColumnMerge[1:], axis=1)
                            RowMerge=[i-1 for i in ColumnMerge]
                            R_merge=[]
                            for i in RowMerge:
                                R_merge.append(R[i])
                            R_new=R_merge[0].copy()
                            for i in range(1,len(R_merge)):
                                R_new=add_poly(R_new,R_merge[i])
                            a=RowMerge[1:]
                            a.sort(reverse=True)
                            for i in a:
                                del R[i]
                            R[RowMerge[0]]=R_new
                    # print(X,R)
        else:
            # activation of edge
            # print('This step is activation of edge',step)
            node1,node2=step   # neighbor, node
            ClassNum_old=X.shape[1]
            SpecialLink=0
            for i in range(len(Link_to_add)):
                if set((node1,node2))==set(Link_to_add[i]):
                    SpecialLink=1
            if SpecialLink==1:
                for Class in range(1,ClassNum_old): 
                    Class_change=X[:,Class].copy()
                    
                    a=Class_change[node1-1] 
                    b=Class_change[node2-1] 
                    if a<b:
                        Class_change[Class_change==b]=Class_change[node1-1]
                    else:
                        Class_change[Class_change==a]=Class_change[node2-1]
                    a=Class_change.copy()
                    a=np.insert(a, 0, values=0, axis=0)
                    ClassNum=len(set(a))
                    
                    for number in range(1,ClassNum+1):
                        if number not in set(a):
                            for j in range(n):
                                if Class_change[j]>number:
                                    Class_change[j]=Class_change[j]-1    
                            
                    X[:,Class]=Class_change
            else:
                R_new=[]
                for Class in range(1,ClassNum_old): 
                    # Get new X (insert new class in case edge is connected)
                    Class_add=X[:,Class*2-1].copy()
                    a=Class_add[node1-1] 
                    b=Class_add[node2-1] 
                    if a<b:
                        Class_add[Class_add==b]=Class_add[node1-1]
                    else:
                        Class_add[Class_add==a]=Class_add[node2-1]
                    
                    # Make group nummbering continuous here, ex: change 0,1,3,1,0,0 to 0,1,2,1,0,0
                    a=Class_add.copy()
                    a=np.insert(a, 0, values=0, axis=0)
                    ClassNum=len(set(a))

                    for number in range(1,ClassNum+1):
                        if number not in set(a):
                            for j in range(n):
                                if Class_add[j]>number:
                                    Class_add[j]=Class_add[j]-1
                    X=np.insert(X, Class*2, values=Class_add, axis=1)
                    
                    # Get new R
                    R_new.append(multiply_poly(R[Class-1],[1,-1]))  # mutiply to 1-p
                    R_new.append(multiply_poly(R[Class-1],[0,1]))   # mutiply to p
                R=R_new
                
            
            # Merge the class with same partition
            for Class in range(1,X.shape[1]):
                ColumnMerge=[]
                ColumnMerge.append(Class)
                for remainClass in range(Class+1,X.shape[1]):
                    if (X[:,Class]==X[:,remainClass]).all():
                        ColumnMerge.append(remainClass)
                if len(ColumnMerge)>=2:
                    X = np.delete(X, ColumnMerge[1:], axis=1)
                    RowMerge=[i-1 for i in ColumnMerge]
                    R_merge=[]
                    for i in RowMerge:
                        R_merge.append(R[i])
                    R_new=R_merge[0].copy()
                    for i in range(1,len(R_merge)):
                        R_new=add_poly(R_new,R_merge[i])
                    a=RowMerge[1:]
                    a.sort(reverse=True)
                    for i in a:
                        del R[i]
                    R[RowMerge[0]]=R_new
            # print(X,R)
    return R
def Link_between_sensor(Sensors,G):
    NumSensor=len(Sensors)
    Link_to_add=[]
    
    check_G=G.subgraph(Sensors)
    listCC = [len(c) for c in sorted(nx.connected_components(check_G), key=len, reverse=True)]
    List=sorted(nx.connected_components(check_G)) # components with node
    minD_node=np.zeros(len(listCC))
    for i in range(len(listCC)):
        component=G.subgraph(List[i])
        degree=np.array(G.degree(List[i]))
        degree=degree[np.lexsort(degree.T)] 
        minD_node[i]=degree[0,0]
        if len(component.edges()) !=0:
            for i in component.edges():
                Link_to_add.append(i)
    for i in range(len(minD_node)-1):
        Link_to_add.append((minD_node[i],minD_node[i+1]))

    Link_to_add=np.array(Link_to_add)
    return Link_to_add

def R_poly(G,Sensors):
    g=G.copy()
    n=G.number_of_nodes()
    Link_to_add=Link_between_sensor(Sensors,G)
    for j in Link_to_add:
        g.add_edge(j[0],j[1])
        
    Adj=nx.adjacency_matrix(g)
    A=Adj.todense()
    network=A.copy()
    PathWidth,Path=FindPath(network,n)
    series=FindSeries(network,Path)
    R=Decomposition(network,n,series,Link_to_add)
    R=R[0] # The all - terminal reliability when placed two sensor
    return R

def main():
    # import network at here 
    s='Real1.txt'
    lines=[]
    # with open('/home/ranxu/TopologyZooNetworks/'+s, 'r') as f:
    with open('D:/TUD/Code/python/Thesis/TopologyZooNetworks/'+s, 'r') as f:
        for line in f.readlines():
            line = line.replace('\n','').replace('\t',' ')
            lines.append(line)

    edges=[]
    for i in range(len(lines)):
        a=list(map(int,lines[i].split()))
        edges.append(a)
    edges=np.array(edges)
    G=nx.Graph()
    for edge in edges:
        G.add_edge(edge[0], edge[1])
    G=nx.convert_node_labels_to_integers(G,first_label=1) # Cannot delete this line !!!
    nx.draw(G, pos=nx.kamada_kawai_layout(G),node_size=300,with_labels = True)
    plt.show()


    Controllers=[1,2]
    reliability=R_poly(G,Controllers)
    print(reliability)

if __name__ == "__main__":
    main()