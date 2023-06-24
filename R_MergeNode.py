import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
import networkx as nx
from decimal import *
from timeit import default_timer as timer
print(getcontext())
Context(prec=28, rounding=ROUND_HALF_EVEN, Emin=-999999, Emax=999999, capitals=1, clamp=0, flags=[], traps=[InvalidOperation, DivisionByZero, Overflow])
getcontext().prec = 50

def FindPath2(network,n):
    nodes=np.array(range(n))
    minAN=n
    minANtotal=n*n
    minPathDecomposition=0
    for node in nodes:
        maxAN=0
        ANtotal=0
        pathDecomposition=[]
        pathDecomposition.append(node) # node is the first node to activate
        remainNodes=nodes[nodes!=node]
        
        while len(remainNodes)!=0:
            a=np.nonzero(network[:,pathDecomposition])
            Neighborhood=list(set(a[0]).intersection(set(remainNodes)))
            random.shuffle(Neighborhood)
            tempAN=n
            
            for neighbor in Neighborhood:
                tempPD=pathDecomposition.copy()
                tempPD.append(neighbor)
                tempRN=remainNodes[remainNodes!=neighbor]
                ActivateNode=len(tempPD)
                for i in tempPD:
                    if network[i,tempRN].any()==0:
                        ActivateNode=ActivateNode-1
                
                if ActivateNode < tempAN:
                    tempAN=ActivateNode
                    chosen=neighbor
            
            if tempAN==len(tempPD):
                connection=0
                for neighbor in Neighborhood:
                    tempConnection=sum(network[pathDecomposition,neighbor])
                    if tempConnection>connection:
                        connection=tempConnection
                        chosen=neighbor
            pathDecomposition.append(chosen)
            remainNodes=remainNodes[remainNodes!=chosen]
            ANtotal=ANtotal+tempAN
            
            if tempAN > maxAN:
                maxAN=tempAN

        if maxAN <= minAN and ANtotal<=minANtotal:
            minAN=maxAN
            minANtotal=ANtotal
            minPathDecomposition=pathDecomposition

    # print('The minimum pathwidth is ',minPathwidth)
    # print('The decomposition path is ',minPathDecomposition) # with node numbering start from 0
    
    return minAN,minPathDecomposition

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

def Decomposition(series,network,Adj_dic,omega):
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
    n=len(network)
    R=[]
    R.append(Decimal(1)) # initial reliability
    X=np.zeros((n,2))
    X[:,0]=range(1,n+1) # Create classes storage array
    X[series[0]-1,1]=1 # First step and activate first node

    for s in range(1,len(series)):  # Start from second step in series
        step=series[s]
        # print(s,'/',len(series),X.shape,len(np.nonzero(X[:,-1])[0]))
        if isinstance(step, (int, np.integer)):    # If it is int number, thus is node process
            if s==len(series)-1:
                R[0]=Decimal(R[0]*omega)
                # print(R[0])
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
                            R_new=R_merge[0]
                            for i in range(1,len(R_merge)):
                                R_new=Decimal(R_new+R_merge[i])
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
            link_p=Adj_dic[str([min(node1-1,node2-1),max(node1-1,node2-1)])]
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
                R_new.append(Decimal(R[Class-1]*(1-link_p)))  # mutiply to 1-p
                R_new.append(Decimal(R[Class-1]*link_p))  # mutiply to p
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
                    R_new=R_merge[0]
                    for i in range(1,len(R_merge)):
                        R_new=Decimal((R_new+R_merge[i]))
                    a=RowMerge[1:]
                    a.sort(reverse=True)
                    for i in a:
                        del R[i]
                    R[RowMerge[0]]=R_new
            # print(X,R)
    return R[0]

def trim(A,Adj_dic,omega):
    change=1
    degree=np.sum(A,1)
    n=len(degree)
    while change==1:
        change=0
        for node in range(n):
            if degree[node]==1:  # degree=1
                neighbors=np.nonzero(A[node])[0]
                neighbor=int(neighbors)
                neighborColumn=A[neighbor]
                if np.sum(neighborColumn)>1: # if neighbor's degree>=2
                    Min=min(node,neighbor)
                    Max=max(node,neighbor)
                    omega=Decimal(omega*Adj_dic[str([Min,Max])])
                    del Adj_dic[str([Min,Max])]
                    degree[node]=0
                    degree[neighbor]=degree[neighbor]-1
                    change=1
                    A[node,neighbor]=0
                    A[neighbor,node]=0
                            
            elif degree[node]==2:
                neighbors=np.nonzero(A[node])[0]
                oldValue=0
                k1=str([min(node,neighbors[0]),max(node,neighbors[0])])
                value1=Adj_dic[k1]
                k2=str([min(node,neighbors[1]),max(node,neighbors[1])])
                value2=Adj_dic[k2]
                k3=str([min(neighbors[0],neighbors[1]),max(neighbors[0],neighbors[1])])
                if k3 in Adj_dic:
                    oldValue=Adj_dic[k3]
                
                omegaTerm=Decimal(1-(1-value1)*(1-value2))
                omega=Decimal(omega*omegaTerm)
                newValue=Decimal(value1*value2/omegaTerm)
                
                if oldValue==0:
                    Adj_dic[str([min(neighbors[0],neighbors[1]),max(neighbors[0],neighbors[1])])]=newValue
                    A[neighbors[0],neighbors[1]]=1
                    A[neighbors[1],neighbors[0]]=1
                else:
                    adjustedValue=Decimal(1-(1-newValue)*(1-oldValue))
                    Adj_dic[k3]=adjustedValue
                    degree[neighbors[0]]=degree[neighbors[0]]-1
                    degree[neighbors[1]]=degree[neighbors[1]]-1
                
                del Adj_dic[k1]
                del Adj_dic[k2]
                degree[node]=0
                A[node,neighbors[0]]=0
                A[neighbors[0],node]=0
                A[node,neighbors[1]]=0
                A[neighbors[1],node]=0
                change=1
    # G_new=nx.Graph()
    # for key,value in Adj_dic.items():
    #     key=key.replace('[', '')
    #     key=key.replace(']', '')
    #     key=list(map(int,key.split(',')))
    #     G_new.add_edge(key[0],key[1])
    #     print(key, value)
    # nx.draw(G_new,with_labels=True,node_size=200)
    # plt.show()
    return A,Adj_dic,omega

def reduction(A,Adj_dic,omega):
    A,Adj_dic,omega=trim(A,Adj_dic,omega)
    # left edges
    new_edges=[]
    L=[]
    for key,value in Adj_dic.items():
        key=key.replace('[', '')
        key=key.replace(']', '')
        new_edges.append(list(map(int,key.split(','))))
        L+=list(map(int,key.split(',')))
    NodeList=list(set(L))
    # number the left node
    node_dic={}
    k=0
    for i in NodeList:
        node_dic[i]=k
        k=k+1
    # creat new Adj
    Adj_new={}
    n=len(NodeList)
    A=np.zeros((n,n))
    for edges in new_edges:
        node1=min(node_dic[edges[0]],node_dic[edges[1]])
        node2=max(node_dic[edges[0]],node_dic[edges[1]])
        Adj_new[str([node1,node2])]=Adj_dic[str(edges)]
        A[node1,node2]=1
        A[node2,node1]=1
    # G_new=nx.Graph()
    # for edge in new_edges:
    #     G_new.add_edge(node_dic[edge[0]],node_dic[edge[1]])
    # nx.draw(G_new,with_labels=True,node_size=200)
    # plt.show()
    # print(new_edges)
    # for key,value in Adj_new.items():
    #     print(key, value)
    return A,Adj_new,omega
 
def MergeSensor(G,SensorPlacement,p):
    n=G.number_of_nodes()
    Adj=nx.adjacency_matrix(G)
    A=Adj.todense()
    Adj_dic={}
    # initialize the link probability
    for i in range(n):
        for j in range(i+1,n,1):
            if A[i,j]==1:
                Adj_dic[str([i,j])]=p
    FirstNode=SensorPlacement[0]-1
    MergeNode=SensorPlacement[1:]
    for node in MergeNode:
        node=node-1
        for i in range(n):
            if A[min(node,i),max(node,i)]==1:
                if i==FirstNode:
                    del Adj_dic[str([min(node,i),max(node,i)])]
                    A[node,i]=0
                    A[i,node]=0
                    continue
                elif A[min(FirstNode,i),max(FirstNode,i)]==1:
                    Adj_dic[str([min(FirstNode,i),max(FirstNode,i)])]=1-(1-p)*(1-Adj_dic[str([min(FirstNode,i),max(FirstNode,i)])])
                else:
                    Adj_dic[str([min(FirstNode,i),max(FirstNode,i)])]=p
                    A[FirstNode,i]=1
                    A[i,FirstNode]=1
                del Adj_dic[str([min(node,i),max(node,i)])]
                A[node,i]=0
                A[i,node]=0
    new_edges=[]
    L=[]
    for key,value in Adj_dic.items():
        key=key.replace('[', '')
        key=key.replace(']', '')
        new_edges.append(list(map(int,key.split(','))))
        L+=list(map(int,key.split(',')))
    NodeList=list(set(L))
    # number the left node
    node_dic={}
    k=0
    for i in NodeList:
        node_dic[i]=k
        k=k+1
    # creat new Adj
    Adj_new={}
    n=len(NodeList)
    A=np.zeros((n,n))
    for edges in new_edges:
        node1=min(node_dic[edges[0]],node_dic[edges[1]])
        node2=max(node_dic[edges[0]],node_dic[edges[1]])
        Adj_new[str([node1,node2])]=Adj_dic[str(edges)]
        A[node1,node2]=1
        A[node2,node1]=1
    # G_new=nx.Graph()
    # for edge in new_edges:
    #     G_new.add_edge(node_dic[edge[0]],node_dic[edge[1]])
    # G_new=nx.convert_node_labels_to_integers(G_new,first_label=1)
    # nx.draw(G_new,with_labels=True,node_size=200)
    # plt.show()
    # for key,value in Adj_new.items():
    #     print(key, value)
    return A,Adj_new
 
def R(G,SensorPlacement,probability):
    omega=Decimal(1)
    p=Decimal(probability)
    A,Adj_dic=MergeSensor(G,SensorPlacement,p)
    A,Adj_dic,omega=reduction(A,Adj_dic,omega)  
    network=A.copy()
    n=len(network)
    PathWidth,Path=FindPath2(network,n)
    # print('PathWidth: ',PathWidth)
    series=FindSeries(network,Path)
    R=Decomposition(series,network,Adj_dic,omega)
    return R

def OS3Egraph():
    g=nx.Graph()
    g.add_edge(9,33)
    g.add_edge(33,32)
    g.add_edge(33,29)
    g.add_edge(21,33)
    g.add_edge(21,27)
    g.add_edge(27,4)
    g.add_edge(1,32)
    g.add_edge(1,29)
    g.add_edge(1,15)
    g.add_edge(29,15)
    g.add_edge(15,23)
    g.add_edge(6,29)
    g.add_edge(6,14)
    g.add_edge(14,5)
    g.add_edge(5,23)
    g.add_edge(6,20)
    g.add_edge(4,20)
    g.add_edge(7,20)
    g.add_edge(7,25)
    g.add_edge(5,25)
    g.add_edge(25,31)
    g.add_edge(4,11)
    g.add_edge(8,11)
    g.add_edge(2,8)
    g.add_edge(2,17)
    g.add_edge(17,31)
    g.add_edge(13,25)
    g.add_edge(18,13)
    g.add_edge(18,19)
    g.add_edge(16,18)
    g.add_edge(2,16)
    g.add_edge(4,30)
    g.add_edge(3,16)
    g.add_edge(30,24)
    g.add_edge(24,26)
    g.add_edge(26,28)
    g.add_edge(28,22)
    g.add_edge(10,22)
    g.add_edge(10,34)
    g.add_edge(12,34)
    g.add_edge(12,30)
    g.add_edge(3,10)
    # g=nx.convert_node_labels_to_integers(g,first_label=1)
    
    return g
    


G=OS3Egraph()
G=nx.convert_node_labels_to_integers(G,first_label=1)
# nx.draw_kamada_kawai(G,with_labels=True,node_size=300,font_size=15,font_color='white')
# plt.show()

probability=0.99
n=G.number_of_nodes()
sensor=[1,2]
reliability=R(G,sensor,probability)
print(reliability)
