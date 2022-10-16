"""
Author: Charlotte Versavel
Date:   June 2022
Last Edit: Oct 2022

                             cluster_class.py

Purpose: a class to represent a cluster of proteins

TODO: code documentation / fxn contracts
TODO: figure out how to use (oct 2022)

"""

import pandas as pd 
import numpy as np
from collections import defaultdict

from sklearn import cluster

class AllClusters:

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * MEMBER VARIABLES * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    clusters = defaultdict(lambda: []) # a dict of relation {cluster_num : list_of_proteins_in_cluster}
    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * INITIALIZERS * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    def __init__(self, csv_filename: str = "", protein_to_cluster_dict: dict() ={}) -> None:
        """  
        Parameters: csv_filename is the name of a csv file containing several 
                    clusters of proteins   
                    protein_to_cluster_dict is a dictionary with the form { protein : cluster_num }
        Purpose:    to populate several single clusters with data from a CSV 
                    file, or from a dictionary
        Returns:    n/a
        """
        self.clusters.clear()

        if csv_filename != "":
            try:
                with open(csv_filename, "r") as data:

                    for item in data:
                        list_of_proteins = item.strip().split("\t")

                        cluster_number = int(list_of_proteins.pop(0))
                        other_number = list_of_proteins.pop(0)

                        self.clusters[cluster_number] = list_of_proteins

            except FileNotFoundError:
                print(f"ERROR! file: {csv_filename} not found.")
        
        elif protein_to_cluster_dict: # dictionary not empty
            for protein in protein_to_cluster_dict.keys():
                self.add_protein_to_cluster(protein, int(protein_to_cluster_dict[protein]))
        
        else: # no filename or dictionary passed in
            print(f"ERROR! please specify a [csv_filename] or a [protein_to_cluster_dict] to initialize the clusters.")
            



    def __repr__(self): 
        """             
        Purpose:    TODO
        Returns:    TODO
        """
        return f"AllClusters has {len(self.clusters)} clusters (use the print_all method to see them)"


    def add_protein_to_cluster(self, protein:str, cluster_num:int) -> None:
        """             
        Parameters: 
            -   protein is the protein to add to a specified cluster
            -   cluster_num is the num of the cluster to add a protein to
        Purpose:    to add a protein to a cluster
        Returns:    n/a
        """
        self.clusters[cluster_num].append(protein)
        # print(f"appended cluster {cluster_num}: {self.clusters[cluster_num]}")

    
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    * * * * * * * * * * * * * * * GETTERS * * * * * * * * * * * * * * * * *  
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

    def get_cluster_proteins(self, cluster_number: int) -> list:
        """             
        Parameters: cluster_number is the number of the cluster to get
        Purpose:    to get the list of proteins from a cluster
        Returns:    the list of proteins in the cluster
        """

        return self.clusters[cluster_number]

    def get_num_clusters(self) -> int:
        """
        Purpose:    to determine the number of clusters
        Returns:    
        """
        return len(self.clusters)

    def get_all_cluster_labels(self) -> list():
        """
        TODO
        """
        print("get_all_cluster_labels fxn has not been tested")
        return self.clusters.keys()

    def get_all_clusters(self) -> dict():
        """
        TODO
        """
        return dict(self.clusters)


    def print_all(self) -> None:
        """             
        Purpose:    to print all the clusters in the dictionary
        Returns:    n/a
        """
        print(self.clusters.keys())
        
        for cluster_num in self.clusters.keys():
            print(f"Cluster {cluster_num}: {self.get_cluster_proteins(cluster_num)}")
    

    # def print_querylist_of_clusters_to_file(self, clusters_to_print: list(), query_filepath: str = "querylist.txt"):
    #     """
    #     clusters_to_print -> specify a list of which clusters to print
    #     TODO
    #     """
    #     output_file = open(query_filepath, 'w')
    #     for cluster_num in clusters_to_print:
    #         for protein in self.get_cluster_proteins(cluster_num):
    #             output_file.write(f"{protein}\tcluster_{cluster_num}\n")
    #     output_file.close()



    def print_to_file(self, filename: str = "output.txt", print_as_dict=False) -> None:
        """
        TODO
        """
        print(f"ERROR! function unfinished: print_to_file")
        output_file = open(filename, 'w')
        if print_as_dict:
            pass
        pass


