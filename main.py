#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 21:09:37 2019

@author: beatrizurdag
"""

import sys 
import re 
import numpy as np
import time
       
class Node(object):

    def __init__(self, data):
        pchr = re.compile("(\d*[pq])")
        pchr_n = re.compile("(\d*)[pq]")
        pchr_a = re.compile("([pq])")
        
        self.id = data[0]
        self.loci = data[1]
        self.pos = data[2]
        self.chr = pchr.findall(self.loci)[0]
        self.chr_n = int(pchr_n.findall(self.loci)[0])
        self.chr_arm = pchr_a.findall(self.loci)[0]
        self.next = None
    
    def __str__(self):
        return ("%s %s %s" %(self.id,self.loci,self.pos))
    
    def __eq__(self, other):
        return (self.id == other.id)
    
    def __lt__(self,other):
       if (self.chr_n != other.chr_n):
            return  (self.chr_n < other.chr_n)
       else:
            return (self.chr_arm < other.chr_arm)
    
    def __le__(self,other):
       if (self.chr_n != other.chr_n):
            return  (self.chr_n < other.chr_n)
       else:
            return (self.chr_arm <= other.chr_arm)
        
class LinkedList(object):
    def __init__(self):
        self.start_node = None
        
    def insert_at_start(self, data):
        """Creates and node and inserts it at the beginning of the list"""
        
        new_node = Node(data)
        new_node.next = self.start_node
        self.start_node= new_node
        
    def traverse(self):
        if self.start_node is None:
            print("Linked List is empty")
            return
        else:
            n = self.start_node
            while n is not None:
                print(n, " ")
                n = n.next
            
    def insert_node_before_item(self,x, new_node):
        """Inserts a node before the element with id == x"""
        
        if self.start_node is None:
            print("List has no element")
            return

        if x == self.start_node.id:
            new_node.next = self.start_node
            self.start_node = new_node
            return

        n = self.start_node
        while n.next is not None:
            if n.next.id == x:
                break
            n = n.next
        if n.next is None:
            print("item not in the list")
        else:
            new_node.next = n.next
            n.next = new_node
            
    def delete_element_by_value(self, x):
        """Deletes a node with id == x"""
        
        if self.start_node is None:
            print("The list has no element to delete")
            return
    
        # Deleting first node 
        if self.start_node.id == x:
            self.start_node = self.start_node.next
            return
    
        n = self.start_node
        while n.next is not None:
            if n.next.id == x:
                break
            n = n.next
    
        if n.next is None:
            print("item not found in the list")
        else:
            n.next = n.next.next
            
    def delete_node(self, to_delete):
        """Deletes a given node"""
        
        if self.start_node is None:
            print("The list has no element to delete")
            return
    
        # Deleting first node 
        if self.start_node == to_delete:
            self.start_node = self.start_node.next
            return
    
        n = self.start_node
        while n.next is not None:
            if n.next == to_delete:
                break
            n = n.next
    
        if n.next is None:
            print("item not found in the list")
        else:
            n.next = n.next.next
    
    def insertion_sort(self):
        """Sorts the linked list using the insertion sort algorithm"""
        
        if self.start_node is None:
            print("Linked List is empty")
            return
        else:
            selected = self.start_node.next
            while selected is not None: 
                cnode = self.start_node
                while ((selected != cnode) and (cnode is not None)):
                    if (selected < cnode):
                        new_selected = selected.next
                        self.delete_node(selected)
                        self.insert_node_before_item(cnode.id,selected)
                        break
                    elif (selected == cnode):
                        break
                    else:
                        cnode = cnode.next
                        new_selected = selected.next  
                selected = new_selected
        return linked_list
    
    def count_close_sequences(self,k,output_filename):
        """Counts the number of sequences pairs that are closer than a given
        threshold for each chromosome arm. Stores this information in a file.
        
        Arguments:
        k - distance threshold
        output_filename - name of the output file"""
        
        if self.start_node is None:
                print("Linked List is empty")
                return
        else:
            fdo = open(output_filename,"w")
            fdo.write("Location \t Number\n")
            fdo.close()
            fdo = open(output_filename,"a")
            n1 = self.start_node
            c_chr1 = n1.chr
            c_count = 0
            while n1 is not None:
                if n1.chr != c_chr1:
                    # Write line for the previous chr and count and restart the process
                    fdo.write("%s\t%s\n" %(c_chr1,c_count))
                    c_chr1 = n1.chr
                    c_count = 0
                
                n2 = n1.next
                found_chr = False
                while n2 is not None:
                    c_chr2 = n2.chr
                    if (c_chr1 == c_chr2):
                        found_chr = True
                        if (euclidean_distance(n1.pos,n2.pos) <= k):
                            c_count += 1
                    else:
                        if(found_chr == True):
                            break
                    n2 = n2.next
                    
                        
                
                n1 = n1.next
            fdo.write("%s\t%s\n" %(c_chr1,c_count))
            fdo.close()


def euclidean_distance(u1,u2):
    """Returns the euclidean distance between 2 bidimensional points"""
    return np.sqrt(np.square(u1[0]-u2[0]) + np.square(u1[1]-u2[1]))


def create_linkedlist_from_file(filename):
    """Creates a Linked List object from a sequence file"""
    
    linked_list = LinkedList()
    fd = open(filename,"r")
    px = re.compile("(\d*\.*\d*),")
    py = re.compile(",(\d*\.*\d*)")
    for line in fd:
        lpieces = line.split()
        # Convert position vector in a tuple
        x = float(px.findall(lpieces[2])[0])
        y = float(py.findall(lpieces[2])[0])
        cdata = (lpieces[0],lpieces[1],(x,y))
        linked_list.insert_at_start(cdata)
    fd.close()
#    linked_list.traverse()
    return linked_list


if __name__ == "__main__":
    start = time.time()
    
    # Parsing arguments and initializing variables
    filename = sys.argv[1]
    k = float(sys.argv[2])
    output_filename = "result4.txt"
    
    # Creating a linked list from the input file
    linked_list = create_linkedlist_from_file(filename) # cost n because I have used insert_at_start
    # Sorting
    linked_list = linked_list.insertion_sort()
    # Counting close sequences in the same arm and writing in output file
    linked_list.count_close_sequences(k,output_filename)
    
    print ("The script took %s s" %(time.time()- start))
       