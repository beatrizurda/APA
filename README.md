# APA
2nd APA assignment

FILE 1: 0.003159046173095703 s 
FILE 2: 0.6112880706787109 s 
FILE 3: 80.56356716156006 s

The running of files 4 and 5 takes a lot of time. 

Considerations I have assumed to improve efficiency:
- The generation of the LinkedList has a O(n) cost because I avoid the 2nd loop by always inserting the nodes in the first position. 
- The sorting of the sequences by the chromosome information reduces the number of iterations in the loops needed to compute the number of pairs that are close within the arms. In principle this cost would be O(n^2) but given that the matrix is symmetric, we only need a triangle minus the principal diagonal. Also, for each node, once the nodes belonging to the same chromosome arm have been reached, the function skips the remaining nodes and continues with the next ones. 

Taking this into account, I consider that a more efficient implementation could imply the following modification of the insertion sort: 
Merging the insertion and deletion of the nodes that have to be moved into a single function called 'move_node_before_item_x'. By doing this, we could use the fact that you always insert earlier in the linked list and delete later on. Therefore, there should be a way of merging these two operations into a single loop, reducing the computational cost of this process (insert and delete) from O(n^2) to O(n).
